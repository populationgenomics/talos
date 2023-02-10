#!/usr/bin/env bash

gsutil -m -u broad-rgp cp -r \
    gs://cpg-broad-rgp-main/mt/3355c9263be6d4b6e13c88b95fb0e3bc1bc99d_1559-broad-rgp.mt \
    gs://cpg-broad-rgp-main-release/mt/3355c9263be6d4b6e13c88b95fb0e3bc1bc99d_1559-broad-rgp.mt
