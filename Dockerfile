FROM ghcr.io/astral-sh/uv:python3.10-bookworm-slim AS base

ENV DEBIAN_FRONTEND=noninteractive

ENV VERSION=8.3.4

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        gnupg \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        openjdk-17-jdk-headless \
        procps \
        wget \
        zip \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM base AS bcftools_compiler

ARG BCFTOOLS_VERSION=1.22

RUN apt-get update && apt-get install --no-install-recommends -y \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        zlib1g-dev && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install

FROM base AS base_bcftools

COPY --from=bcftools_compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=bcftools_compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

FROM base_bcftools AS base_bcftools_echtvar

ARG ECHTVAR_VERSION=v0.2.2

ADD "https://github.com/brentp/echtvar/releases/download/${ECHTVAR_VERSION}/echtvar" /bin/echtvar

RUN chmod +x /bin/echtvar

FROM base_bcftools_echtvar AS talos

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

WORKDIR /talos

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project --no-dev

# Add in the additional requirements that are most likely to change.
COPY LICENSE pyproject.toml uv.lock README.md .
COPY src src/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install ".[cpg]"

# Place executables in the environment at the front of the path
ENV PATH="/talos/.venv/bin:$PATH"
