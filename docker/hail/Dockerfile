# Use a Python image with uv pre-installed
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

RUN apt update && apt install -y --no-install-recommends \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        gnupg \
        openjdk-17-jdk-headless \
        wget \
        zip && \
    apt clean

# Install the project into `/app`
WORKDIR /talos

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

ADD uv.lock .

# Install dependencies using the requirements file
RUN uv pip sync --system uv.lock

WORKDIR /talos

COPY *.py /talos/

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []
