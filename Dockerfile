# Use ARG for Python version flexibility
ARG PYTHON_VERSION=3.11.9
FROM python:${PYTHON_VERSION}-slim

# Metadata
LABEL fly_launch_runtime="flask"

# Set working directory
WORKDIR /app

# Install system dependencies required for building wheels
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        wget \
        git \
        gcc \
        g++ \
        libssl-dev \
        libffi-dev \
        && rm -rf /var/lib/apt/lists/*

# Upgrade pip for latest wheel support
RUN pip install --upgrade pip setuptools wheel

# Copy requirements first for caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Expose port
EXPOSE 8080

# Run app using Gunicorn
CMD ["gunicorn", "-b", "0.0.0.0:8080", "app:app", "--workers", "1", "--timeout", "120"]
