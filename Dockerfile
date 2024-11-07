FROM python:3.11-bullseye AS base

# gcloud (default) and none are availabe choices
ARG cloud="gcloud"

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

FROM base AS talos_gcloud

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
RUN curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

# Add in the additional requirements that are most likely to change.
COPY requirements*.txt README.md setup.py ./
COPY src src/
RUN pip install --upgrade pip && pip install .[cpg]


FROM base AS talos_none

RUN echo "Skipping cloud dependency installation"

# Add in the additional requirements that are most likely to change.
COPY requirements*.txt README.md setup.py ./
COPY src src/
RUN pip install .[cpg]
