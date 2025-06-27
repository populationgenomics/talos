# Mock Data

This repository contains a synthetic data generator, which can be used to create mock data for testing purposes. This is useful for unit tests, category design, and testing inheritance models (MOI).

This Mock data generator contains two elements:

- [the data model](../src/talos/data_model.py), containing the various classes and transformations required to create mock data
- [an implementation script](../nextflow/inputs/generate_test_data.py), which can be used to create mock data in a Hail Table or MatrixTable format, and export it to a VCF file.

The data model was initially designed to create Talos-compatible data, emulating a highly specific Hail-centric VEP annotation pipeline. It can still do that, but a more typical usage pattern would be generating a generic VCF with minimal annotation. [This script](../nextflow/inputs/generate_test_data.py) is an example of generating the minimal VCF, and is used to generate the dummy data used in this repository's demonstration workflow.

The data model provides granular control over samples in the VCF, genotypes and variant phase for each, and can be used to insert arbitrary INFO and ENTRY fields. The data is generated with Hail, and exports VCFs with fully populated headers, correct for the expected reference genome.
