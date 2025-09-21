# ================================
# Stage 1: Builder
# ================================
FROM python:3.11-slim AS builder
WORKDIR /app

# Install build tools only
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        python3-dev \
        git \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

# Install system-wide (no --user), so cleanup works
RUN pip install --no-cache-dir -r requirements.txt \
    && rm -rf /root/.cache/pip \
    # Strip out large test directories to save space
    && find /usr/local/lib/python3.11/site-packages/pyarrow -name "tests" -type d -exec rm -rf {} + \
    && find /usr/local/lib/python3.11/site-packages/numpy -name "tests" -type d -exec rm -rf {} + \
    && find /usr/local/lib/python3.11/site-packages/pandas -name "tests" -type d -exec rm -rf {} + \
    && find /usr/local/lib/python3.11/site-packages -name "__pycache__" -type d -exec rm -rf {} +

# ================================
# Stage 2: Runtime
# ================================
FROM python:3.11-slim
WORKDIR /app

# Copy installed packages from builder
COPY --from=builder /usr/local /usr/local

# Copy your app code
COPY . .

# Environment
ENV PORT=8080
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Run with gunicorn
CMD ["gunicorn", "-b", "0.0.0.0:8080", "app:app", "--workers=1", "--timeout=120"]
