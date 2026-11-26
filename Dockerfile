# The single Talos image - the Python package, bcftools/htslib, echtvar, and SVAFotate.
#
# SVAFotate used to be built as a separate image because its dependency stack could not be
# reconciled with Hail's. Hail is no longer a dependency, but then the dependency constraint
# was python < 3.12
# As a result SVAFotate has been forked and rewritten, we install that version.
#
# build with "docker build -t talos:13.0.0 ."

FROM python:3.12-slim-trixie AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        procps \
        wget \
        zip \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM base AS compiler

ARG BCFTOOLS_VERSION=1.24

RUN apt-get update && apt-get install --no-install-recommends -y \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install && \
    cd htslib-${BCFTOOLS_VERSION} && \
    make && \
    make DESTDIR=/bcftools_install install

FROM base AS talos

COPY --from=compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

ARG ECHTVAR_VERSION=v0.2.2
ENV VERSION=13.0.0

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

# Place executables in the environment at the front of the path
ENV PATH="/talos/.venv/bin:$PATH"

# Add in the additional requirements that are most likely to change.
COPY LICENSE pyproject.toml uv.lock README.md ./
COPY src src/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --all-extras --frozen --no-dev

# both halves of the image have to survive each other's install - the Talos package and its cyvcf2,
# and the SVAFotate CLI. Cheap, and the failure mode it catches is otherwise a runtime one
RUN python -c "import talos, cyvcf2" && \
    bcftools --version && \
    tabix --version && \
    echtvar --version && \
    svafotate annotate -h

COPY echtvar echtvar/
