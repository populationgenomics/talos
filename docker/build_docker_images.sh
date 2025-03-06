#!/usr/bin/env bash

# Build the docker images for the annotation workflow

set -euo pipefail

# build the bcftools docker image
docker build -t bcftools:1.21 bcftools

# build the echtvar docker image
docker build -t echtvar:v0.2.1 echtvar

# build the docker image for the hail stage of the workflow - hail install + scripts
docker build -t hail:0.2.133 hail
