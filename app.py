from flask import Flask, render_template, request
import requests
import os
import io
import base64
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pandas as pd
from alphagenome.data import genome
from alphagenome.data.transcript import Transcript
from alphagenome.data import transcript as transcript_utils
from alphagenome.data import gene_annotation
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
BASE_URL = "https://github.com/Pamela-Hao/fortunatecow/releases/download/v1.0/"
LOCAL_DIR = "/app/gencode_split"

def get_chr_file(chrom):
    chrom = chrom.strip()
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    os.makedirs(LOCAL_DIR, exist_ok=True)
    local_path = os.path.join(LOCAL_DIR, f"{chrom}.feather")
    if not os.path.exists(local_path):
        url = f"{BASE_URL}/{chrom}.feather"
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
        with open(local_path, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1024*1024):
                if chunk:
                    fh.write(chunk)
    return local_path

app = Flask(__name__)


# You should store your AlphaGenome API key as an environment variable for safety
ALPHAGENOME_API_KEY = os.getenv("ALPHAGENOME_API_KEY")
model = dna_client.create(ALPHAGENOME_API_KEY)

@app.route("/", methods=["GET", "POST"])
def index():
    plot_url = None
    prediction = None
    error = None

    if request.method == "POST":
        chrom = request.form.get("chromosome", "").strip()
        pos = request.form.get("position")
        ref = request.form.get("ref")
        alt = request.form.get("alt")
        output_type = request.form.get("output_type") 
        window_size = request.form.get("window_size", type=int) or 2**15
        prediction_center = request.form.get("prediction_center", type=int) 
        ontology_terms_str = request.form.get("ontology_terms")
        ontology_terms = [term.strip() for term in ontology_terms_str.split(",")] if ontology_terms_str else ["UBERON:0000955"]
        ref_color = request.form.get("ref_color")
        figure_ref_color = ref_color if ref_color else "dimgrey"
        alt_color = request.form.get("alt_color")
        figure_alt_color = alt_color if alt_color else "red"
        if not chrom or not pos or not ref or not alt or not output_type:
            error = "Please fill in all fields."
        else:
            gtf_region = None
            pos = int(pos)
            center_pos = prediction_center if prediction_center else pos
            start = center_pos - window_size // 2
            end = center_pos + window_size // 2
            figure_interval = genome.Interval(chrom, int(start), int(end))
            gtf_region = None
            chr_file = get_chr_file(chrom)
            dataset = ds.dataset(chr_file, format="feather")
            for batch in dataset.to_batches():
                mask = pc.and_(
                pc.less_equal(batch["Start"], end),
                pc.greater_equal(batch["End"], start)
                )
                filtered = batch.filter(mask)
                if filtered.num_rows > 0:
                    gtf_region = filtered.to_pandas()
                    break
            del batch, filtered, dataset
            import gc; gc.collect()
                # Define variant and interval
            variant_obj = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
                )
            interval_obj = variant_obj.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
                # Convert output types to AlphaGenome OutputType objects

            columns = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'gene_id', 'gene_type', 'gene_name', 'level', 'tag', 'transcript_id', 'transcript_type', 'transcript_name', 'transcript_support_level', 'havana_transcript', 'exon_number', 'exon_id', 'hgnc_id', 'havana_gene', 'ont', 'protein_id', 'ccdsid', 'artif_dupl']
            requested_outputs = [getattr(dna_client.OutputType, output_type)]
                # Predict
            outputs = model.predict_variant(
                interval=interval_obj,
                variant=variant_obj,
                requested_outputs=requested_outputs,
                ontology_terms = ontology_terms
                )
            plot_items = []
            longest_transcripts = None
            if gtf_region is not None and not gtf_region.empty:
                gtf_transcripts = gene_annotation.filter_protein_coding(gtf_region)
                if not gtf_transcripts.empty and gtf_transcripts['Feature'].str.contains('exon').any():
                    gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
                    transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)
                    longest_transcripts = transcript_extractor.extract(figure_interval)
                    plot_items.append(plot_components.TranscriptAnnotation(longest_transcripts))
                    exons_exist = True
                else:
                    exons_exist = False
                    gene_metadata = gtf_transcripts[['gene_name', 'Start', 'End', 'Strand']].drop_duplicates(subset = 'gene_name')

            ref_track = getattr(outputs.reference, output_type.lower())  # make sure casing matches attribute
            alt_track = getattr(outputs.alternate, output_type.lower())
            

            plot_items.append(
            plot_components.OverlaidTracks(
            tdata={'REF': ref_track, 'ALT': alt_track},
            colors={'REF': figure_ref_color, 'ALT': figure_alt_color}
            )
            )
            fig = plot_components.plot(
                plot_items,
                interval=figure_interval,
                annotations=[plot_components.VariantAnnotation([variant_obj], alpha=0.8)],
            )
            if not exons_exist and gtf_region is not None and not gtf_region.empty:
                ax = fig.axes[0]  # assume first axis contains the track plot
                x_min, x_max = ax.get_xlim()
                y_min, y_max = ax.get_ylim()
                y_offset = 0.05 * (y_max - y_min)
                y_arrow = y_max + y_offset
                for _, row in gene_metadata.iterrows():
                    start = max(row['Start'], x_min)
                    end = min(row['End'], x_max)
                    if start >= end:
                        continue
                    y_frac = 1.0       # y position in axes fraction
                    text_offset = -5
                    if row['Strand'] == '+':
                        ax.annotate('',
                                    xy=(end, 1.0), xytext=(start, 1.0),
                                    arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5),
                                    xycoords=('data', 'axes fraction'))  # y in fraction of axis
                        ax.text(start, y_frac, row['gene_name'],
                                ha='right', va='bottom', fontsize=8,
                                transform=ax.get_xaxis_transform())
                    else:
                        ax.annotate('',
                                    xy=(start, 1.0), xytext=(end, 1.0),
                                    arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5),
                                    xycoords=('data', 'axes fraction'))
                        ax.text(start, y_frac, row['gene_name'],
                                ha='right', va='bottom', fontsize=8,
                                transform=ax.get_xaxis_transform())
            print("Figure created")
            # Save the figure to a PNG in memory
            img = io.BytesIO()
            fig.savefig(img, format='png', bbox_inches='tight')
            img.seek(0)
            plot_url = base64.b64encode(img.getvalue()).decode()  # convert to base64 string
    return render_template("index.html", plot_url=plot_url, error=error)

                

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)))
    #app.run(debug=True)

