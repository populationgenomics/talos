FROM python:3.10-bullseye

RUN apt update && apt install -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    # Google Cloud SDK: use the script-based installation, as the Debian package is outdated.
    curl https://sdk.cloud.google.com > install.sh && \
    bash install.sh --disable-prompts --install-dir=/opt && \
    rm install.sh

ENV PATH=$PATH:/opt/google-cloud-sdk/bin

COPY requirements*.txt .

RUN python3 -m pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY src/talos talos/
RUN pip install .[cpg]
