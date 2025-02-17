FROM python:3.10-bullseye AS base
COPY --from=ghcr.io/astral-sh/uv:0.5.5 /uv /uvx /bin/

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    apt clean

# install nextflow
ADD https://get.nextflow.io nextflow
RUN chmod +x nextflow && \
    mv nextflow /usr/bin && \
    nextflow self-update

WORKDIR /talos

FROM base AS talos_gcloud

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
ADD https://sdk.cloud.google.com install_glcoud.sh
RUN bash install_glcoud.sh --disable-prompts --install-dir=/opt && \
    rm install_glcoud.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

# Add in the additional requirements that are most likely to change.
COPY README.md pyproject.toml uv.lock ./

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project

COPY src src/

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --extra cpg

# Place executables in the environment at the front of the path
ENV PATH="/talos/.venv/bin:$PATH"

FROM base AS talos_none

RUN echo "Skipping cloud dependency installation"

# Add in the additional requirements that are most likely to change.
COPY README.md pyproject.toml uv.lock ./

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project

COPY src src/

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --extra cpg

# Place executables in the environment at the front of the path
ENV PATH="/talos/.venv/bin:$PATH"
