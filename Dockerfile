FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libhdf5-dev \
        tabix \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

# gdown downloads the Nextflow --premade resources (PREMADE_DOWNLOAD), so it
# is installed here rather than from PyPI on the compute node at run time.
RUN pip install --no-cache-dir . gdown==6.1.0

# Links the GHCR package to its source repository.
LABEL org.opencontainers.image.source="https://github.com/AlexandrovLab/DupCaller" \
      org.opencontainers.image.description="DupCaller: somatic mutation calling from barcoded ecNGS (duplex sequencing) data" \
      org.opencontainers.image.licenses="BSD-2-Clause"
