# Use NVIDIA CUDA base image with Python
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    BOLTZ_CACHE="/opt/boltz_cache"

# Install system dependencies
RUN apt-get update && apt-get install -y build-essential
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    git \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a symbolic link for python
RUN ln -sf /usr/bin/python3 /usr/bin/python

# Create boltz cache directory
RUN mkdir -p ${BOLTZ_CACHE}

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Install Boltz in a fresh environment
RUN pip install --no-cache-dir boltz

# Create a working directory
WORKDIR /workspace

# Copy example input files (if needed)
# COPY examples/ /workspace/examples/

# Copy the launching script run.sh 
COPY run.sh /workspace/run.sh

ENTRYPOINT ["/bin/sh"]
