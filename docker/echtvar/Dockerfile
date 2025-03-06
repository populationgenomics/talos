FROM debian:bookworm-slim

ARG VERSION=${VERSION:-v0.2.1}

ADD "https://github.com/brentp/echtvar/releases/download/${VERSION}/echtvar" /bin/echtvar

RUN chmod +x /bin/echtvar

COPY config.json /config.json

ENV ECHTVAR_CONFIG="/config.json"

# estimated total size of this image is ~35MB
