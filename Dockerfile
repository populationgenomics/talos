FROM python:3.11-slim-bullseye AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        gnupg \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        openjdk-11-jdk-headless \
        procps \
        wget \
        zip \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM base AS bcftools_compiler

ARG BCFTOOLS_VERSION=1.23.1

# AS OF 11.0.0, Talos is building BCFtools from a private fork. This fork contains a single change - csq applies annotations
# to both coding and non-coding genes in the event of overlapping genes. By default BCFtools skips non-coding gene annotation
# if a coding transcript consequence was detected, but in practice this is masking clinically relevant non-coding gene variation
# in cases where the non-coding gene overlaps with a non-clinically relevant coding gene.
# The change made is to always search for non-coding variation, but to mask non-coding consequences in otherwise coding transcripts
# if a coding change was already detected.
RUN apt-get update && apt-get install --no-install-recommends -y \
    autoconf \
    gcc \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    make \
    zlib1g-dev && \
    wget https://github.com/samtools/htslib/releases/download/${BCFTOOLS_VERSION}/htslib-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf htslib-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd htslib-${BCFTOOLS_VERSION} && \
    ./configure --enable-libcurl && \
    make && \
    make DESTDIR=/bcftools_install install && \
    cd .. && \
    git clone https://github.com/populationgenomics/bcftools.git && \
    cd bcftools && \
    autoheader && \
    autoconf && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install

FROM base AS talos

COPY --from=bcftools_compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=bcftools_compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/
COPY --from=bcftools_compiler /bcftools_install/usr/local/lib/ /usr/local/lib/
RUN ldconfig

ARG ECHTVAR_VERSION=v0.2.2
ARG VERSION=11.0.3

RUN wget -q -O /bin/echtvar "https://github.com/brentp/echtvar/releases/download/${ECHTVAR_VERSION}/echtvar" && \
    chmod +x /bin/echtvar

COPY --from=ghcr.io/astral-sh/uv:0.9.26 /uv /uvx /bin/

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

WORKDIR /talos

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --all-extras --frozen --no-install-project --no-dev

# Add in the additional requirements that are most likely to change.
COPY LICENSE pyproject.toml uv.lock README.md ./
COPY src src/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --all-extras --frozen --no-dev

# Place executables in the environment at the front of the path
ENV PATH="/talos/.venv/bin:$PATH"

COPY echtvar echtvar/
