from flask import Flask, render_template, request
import requests
import os
import io
import base64
import pyarrow.feather as feather
import pyarrow.ipc as ipc
import pyarrow.compute as pc
import pyarrow as pa
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import pandas as pd
GTF_URL = "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather"
LOCAL_GTF = "gencode.v46.annotation.gtf.feather"
if not os.path.exists(LOCAL_GTF):
    print("Downloading GTF Feather file...")
    resp = requests.get(GTF_URL, stream=True)
    with open(LOCAL_GTF, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1024*1024):
            f.write(chunk)
    print("Download complete.")

app = Flask(__name__)


# You should store your AlphaGenome API key as an environment variable for safety
ALPHAGENOME_API_KEY = os.getenv("ALPHAGENOME_API_KEY")
model = dna_client.create(ALPHAGENOME_API_KEY)
def query_gtf_region(chrom, start, end):
    """
    Yields Pandas DataFrames for rows matching the chromosome and interval.
    Reads Feather file in batches to avoid loading everything into memory.
    """
    with pa.memory_map(LOCAL_GTF, "r") as source:
        reader = ipc.RecordBatchFileReader(source)
        for i in range(reader.num_record_batches):
            batch = reader.get_batch(i)
            table = pa.Table.from_batches([batch])
            mask = pc.and_(
                pc.equal(table['Chromosome'], chrom),
                pc.less_equal(table['Start'], end),
                pc.greater_equal(table['End'], start)
            )
            region = table.filter(mask)
            if len(region) > 0:
                yield region.to_pandas()


@app.route("/", methods=["GET", "POST"])
def index():
    plot_url = None
    prediction = None
    error = None

    if request.method == "POST":
        chrom = request.form.get("chromosome")
        pos = request.form.get("position")
        ref = request.form.get("ref")
        alt = request.form.get("alt")
        output_type = request.form.get("output_type") 
        window_size = request.form.get("window_size", type=int) or 2**15
        ontology_terms_str = request.form.get("ontology_terms")
        ontology_terms = [term.strip() for term in ontology_terms_str.split(",")] if ontology_terms_str else ["UBERON:0000955"]

        if not chrom or not pos or not ref or not alt or not output_type:
            error = "Please fill in all fields."
        else:
            pos = int(pos)
            start = pos - window_size // 2
            end = pos + window_size // 2
            region_dfs = list(query_gtf_region(chrom, start, end))
            if not region_dfs:
                error = "No transcripts found in the selected interval."
                return render_template("index.html", error=error)

            gtf_region = pd.concat(region_dfs, ignore_index=True)
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
            gtf_transcripts = gene_annotation.filter_protein_coding(gtf_region)
            gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
            transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)
            longest_transcripts = transcript_extractor.extract(interval_obj)

            ref_track = getattr(outputs.reference, output_type.lower())  # make sure casing matches attribute
            alt_track = getattr(outputs.alternate, output_type.lower())
            interval = variant_obj.reference_interval.resize(window_size)
            fig = plot_components.plot(
                [
                    plot_components.TranscriptAnnotation(longest_transcripts),
                    plot_components.OverlaidTracks(
                        tdata={
                            'REF': ref_track,
                            'ALT': alt_track
                        },
                        colors={'REF': 'dimgrey', 'ALT': 'red'},
                    ),
                ],
                interval=interval,
                annotations=[plot_components.VariantAnnotation([variant_obj], alpha=0.8)],
            )
            print("Figure created")
            # Save the figure to a PNG in memory
            img = io.BytesIO()
            fig.savefig(img, format='png', bbox_inches='tight')
            img.seek(0)
            plot_url = base64.b64encode(img.getvalue()).decode()  # convert to base64 string

    return render_template("index.html", plot_url=plot_url, error=error)

                

if __name__ == "__main__":
    #app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)))
    app.run(debug=True)

