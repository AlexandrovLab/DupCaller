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

RUN pip install --no-cache-dir .
