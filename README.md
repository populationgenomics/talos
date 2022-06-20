# Genome-wide Variant Prioritisation (MVP)

Code Repository for the Rare Disease variant prioritisation pipeline

Currently being ported over from [a prior repository](https://github.com/populationgenomics/rare-disease/tree/initial_content)

review [design docs](design_docs) for component descriptions

## Flow

The current implementation consists of the following steps:

### [Annotation](reanalysis/annotation.py)

Summarised here [Annotation stage](design_docs/Annotation.md)

This stage can take a MT as input, but a pre-processing step will revert to VCF prior to running annotation

#### Inputs

- String: path to a VCF
- (Implicit): requires annotation resources set up in expected format (waiting for more documentation)

#### Process

1. We calculate a number of intervals across the genome, so that VCF annotation can be split and parallelised
2. For each interval, we use VEP to annotate the variant data, writing output in JSON format
3. The interval annotation data is parsed using Hail, based on a set VEP schema, creating a Hail Table of annotations
4. Once all annotation jobs are complete, the temporary tables are joined into one large table, containing all the consequence annotations
5. The original VCF is loaded into a Hail MatrixTable, and in a single loop we annotate each variant with population frequencies, variant consequences, and clinvar records

#### Outputs

The resultant MT containing all annotations, written as a `.mt` directory

### [PanelApp Query](reanalysis/query_panelapp.py)

Summarised here [PanelApp stage](design_docs/PanelApp_interaction.md)

#### Inputs

- String: a panelapp panel ID
- Optional(String): path to gene list file
- Optional(String): prior panel version

#### Process

1. Query PanelApp for the latest version of the defined panel
2. Store all Green gene data (ENSG, Symbol, MOI)
3. Optional (if prior data is provided) annotate each current gene entity with 'New' if it has been added since that
original data

#### Outputs

 - JSON: PanelApp data for use downstream

### [Labelling in Hail](reanalysis/hail_filter_and_label.py)

Summarised here [Labelling stage](design_docs/Hail_Filter_and_Label.md)

#### Inputs

- MatrixTable: Annotated Variant Data
- String: path to configuration file (JSON)
- String: Path to PanelApp data (JSON)


#### Process

1. Removal of all variants not on PanelApp 'Green' genes
2. Filtering of the Annotated data to remove common or low quality variants
3. Labelling of individual variants with categories of interest
4. Removal of all un-labelled variants
5. Construction of CSQ strings (see below)
6. Search for compound-heterozygous pairs within the remaining variants

#### Outputs

- VCF: labelled variants
- JSON: compound-het pairs

### [VCF Re-heading/CSQ String](reanalysis/interpretation_runner.py)

During the Hail labelling stage, we [create a compound CSQ string](https://github.com/populationgenomics/automated-interpretation-pipeline/blob/main/reanalysis/hail_filter_and_label.py#L478-L571)
from a number of separate annotations. This format mirrors the way that VEP annotates VCF files, and creates one or more
compound (|-delimited) strings. During this process we take instruction from the config file to decide which values to
retain, and the order to squash them together in.

When Hail writes the variant data out as a VCF file, it creates a generic entry in the header, stating that the CSQ
field is a String, but provides no way to decode it. So that the VCF can be parsed without re-referencing the config
file, we alter the default Header line to contain the full CSQ compound string. This was originally used to enable
compound-het calculation by Slivar, and is retained as it could possibly be useful for other downstream users/purposes.

[relevant Slivar issue](https://github.com/brentp/slivar/issues/120)


### [MOI Checking](reanalysis/validate_categories.py)

##### Inputs

- String: path to runtime configuration/settings (JSON)
- VCF: labelled variants
- String: path to compound-heterozygous pairs (JSON)
- String: path to PanelApp data (JSON)
- String: path to PED file for cohort

#### Process

1. iterate through the VCF, pulling out the labels & annotations for each variant
2. for each variant in turn:
   1. look up the relevant gene in the PanelApp data to find the MOI to use
   2. for each participant which has this variant, see if the MOI fits with their family structure
   3. If this search succeeds record the variant,  participant, & passing MOI

#### Outputs

- JSON: a list of all confirmation records
