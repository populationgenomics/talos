FROM python:3.13-bullseye

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
RUN curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

# install nextflow
ADD https://get.nextflow.io nextflow
RUN chmod +x nextflow && \
    mv nextflow /usr/bin && \
    nextflow self-update

# # then finally, install Talos
COPY requirements*.txt README.md setup.py .
COPY src src/
RUN pip install --upgrade pip \
    && pip install .[cpg]
