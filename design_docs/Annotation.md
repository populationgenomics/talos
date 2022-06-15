# Annotation

Prior to running the filtering and categorisation process, we need the variant data to be annotated; population
frequencies, clinvar entries, and a range of predicted consequences. To do this we use a combination of VEP and
annotation directly from Hail Tables containing annotations.

This process is split into a number of steps:

1. If the variant (singleton, or joint call) data is provided in MatrixTable format, we revert this back to a VCF
2. We use GATK to identify intervals across the genome, so the VCF can be evenly split and annotation tasks parallelised
3. For each VCF interval, we spawn a separate job using CLI VEP to annotate the variant data, writing JSON output
4. Each interval's annotations are parsed using Hail, based on the VEP schema, creating a Hail Table of annotations
5. Once all annotation jobs are complete, interval tables are joined into one, containing all consequence annotations
6. The original VCF is loaded into a Hail MatrixTable, and in a single loop we annotate each variant with
   - population frequencies
   - variant consequences
   - clinvar records
7. The final MatrixTable, containing all variant data and all annotations, is written as an output directory

## Inputs

Annotation is currently built around the resource bundle structure defined in [hail-elasticsearch](https://github.com/broadinstitute/hail-elasticsearch-pipelines).
This allows us to leverage the Broad's resource packaging tools, whilst allowing us the flexibility to update and alter
data sources as & when we see fit.

The annotation resource structure is currently set rigidly, using the following structure:

```text
reference_location/
├─ seqr/
│  ├─ clinvar/
│  │  ├─ clinvar_hail_table
│  ├─ combined_reference/
│  │  ├─ reference_hail_table
VEP_annotation/
├─ VEP_annotation_hail_table
```

*VEP_annotation* here is the output location of a previous step (#3 in the list above), where the CLI VEP
intervals are combined into a single object.

*Reference_Location* here is set from the global runtime configuration (see [cpg-utils.config](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/config.py)),
on the path:

```python
from cpg_utils.config import get_config

_reference_root = get_config()['workflow']['reference_prefix']
```

The annotation step takes the VEP annotations in HailTable form, and combines that data with the
pre-generated [clinvar](https://github.com/broadinstitute/hail-elasticsearch-pipelines/blob/master/hail_scripts/utils/clinvar.py)
and [combined](https://github.com/broadinstitute/hail-elasticsearch-pipelines/blob/master/download_and_create_reference_datasets/v02/hail_scripts/write_combined_reference_data_ht.py)
annotations (in silico tool scores, population frequencies), to form the fully annotated joint-call MatrixTable
