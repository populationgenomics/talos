ARG PY_VER=${PY_VER:-3.10}
ARG UV_VER=${UV_VER:-0.6.4}

FROM ghcr.io/astral-sh/uv:${UV_VER}-python${PY_VER}-bookworm-slim AS basic

RUN apt update && apt install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-17-jre-headless \
        wget \
        zip && \
    apt clean

FROM base AS talos_gcloud

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
ADD https://sdk.cloud.google.com install_glcoud.sh
RUN bash install_glcoud.sh --disable-prompts --install-dir=/opt && \
    rm install_glcoud.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

# Add in the additional requirements that are most likely to change.
COPY requirements*.txt README.md setup.py ./
COPY src src/
RUN pip install --upgrade pip && pip install --no-cache-dir .[cpg]

FROM base AS talos_none

RUN echo "Skipping cloud dependency installation"

# Add in the additional requirements that are most likely to change.
COPY requirements*.txt README.md setup.py ./
RUN pip install -r requirements.txt
COPY src src/
RUN pip install --no-cache-dir .
