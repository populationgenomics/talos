FROM python:3.10-bullseye

ARG HAIL_SHA=${HAIL_SHA:-8f6797b033d2e102575c40166cf0c977e91f834e}

RUN apt update && apt install -y \
        apt-transport-https \
        bash \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        g++ \
        gcc \
        git \
        gnupg \
        liblapack3 \
        libopenblas-base \
        make \
        openjdk-11-jdk-headless \
        rsync \
        software-properties-common \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# Install Hail from the CPG fork.
RUN git clone https://github.com/populationgenomics/hail.git && \
    cd hail && \
    git checkout $HAIL_SHA && \
    cd hail && \
    # Install locally, avoiding the need for a pip package.
    make install && \
    cd ../.. && \
    rm -rf hail

COPY . /
RUN pip install .
