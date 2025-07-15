# Configuration File

The configuration file is a TOML file, containing headed sections for each principal step of Talos.  This config file should be created in full prior to running Talos stages, a template is present at [example_config.toml](../src/talos/example_config.toml). This page goes into more detail about how we arrived at specific default settings, and where possible, describes a reasonable range of alternatives.

## Stage agnostic
  * Config Section: `categories`
  * Description: Maps the minimised category label to a longer description, presented in the final report

  * Config Section: `result_history`
  * Description: If provided `result_history` indicates a directory, local or cloud, which contains incremental result history files. With each run the most recent 'history' will be loaded in, and used to annotate the current result set with the earliest date a result was seen/classified. If a novel classification is made, or evidence changes, the date annotated into the report will be shifted to the current run date, assisting in the filtering/prioritisation process.

## Stages: `DownloadPanelApp`, `UnifiedPanelAppParser`
  * Config Section: `GeneratePanelData`
  * Description: Controls which PanelApp deployment to use, which default Panel to use as the base for analysis, which additional genes/panels to include, and which genes to blacklist from the analysis

| Field                 | Purpose                                                                 | Default                                | Alternative                                    |
|-----------------------|-------------------------------------------------------------------------|----------------------------------------|------------------------------------------------|
| `panelapp`            | PanelApp server to query                                                | https://panelapp-aus.org/api/v1/panels | https://panelapp.genomicsengland.co.uk/        |
| `default_panel`       | Base panel for the analysis, applied to app participants                | 137, the PanelApp Aus Mendeliome       |                                                |
| `require_pheno_match` | Genes to remove from analysis unless phenotype matched to a participant | FLG, GJB2, F2, F5, HFE                 |                                                |
| `forbidden_genes`     | Genes to remove from analysis                                           | N/A                                    | Late onset/incidental findings                 |
| `forced_panels`       | Panels to apply to all cohort members                                   | N/A                                    | Phenotype specific panels                      |
| `within_x_months`     | Integer, genes added to panels within X months are prioritised          | 6                                      |                                                |
| `manual_overrides`    | A manually defined panel: gene symbols, MOI                             | N/A                                    | Recently discovered genes, not yet in PanelApp |

## Stage: `RunHailFiltering`
  * Config Section: `RunHailFiltering`
  * Description: Controls filtering thresholds and de novo variant detection parameters during the Hail Stage

| Field                     | Purpose                                                                                                                                                                                                                                    | Default                                                                       | Alternative                                                                                        |
|---------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| `csq_string`              | list of named fields to concatenate into the VCF export                                                                                                                                                                                    | see example config                                                            | reduce/extend based on the available annotations                                                   | 
| `critical_csq`            | High impact Consequence terms                                                                                                                                                                                                              | frameshift, splice_acceptor, splice_donor, start_lost, stop_gained, stop_lost | if using VEP annotations, use altered VEP equivalents (e.g. "frameshift" -> "frameshift_variant")  |
| `additional_csq`          | During de novo variant filtering, we reduce the variant set to critical consequences + this list                                                                                                                                           | missense                                                                      |                                                                                                    |
| `af_semi_rare`            | A coarse population frequency filter - variants more frequent than this are removed from consideration. Later in the process more rigorous filters are applied, so this is a 'mild' filter to retain plausible compound-het candidates     | 0.01                                                                          |                                                                                                    |
| `callset_af_sv_recessive` | A coarse population and callset frequency filter - if an SV variant is more common than this either within the callset, or in the gnomAD annotations, remove                                                                               | 0.03                                                                          |                                                                                                    |

  * Config Section: `de_novo`
  * Description: Controls filtering thresholds and de novo variant detection parameters during the Hail Stage

| Field                    | Purpose                                                                                                                                                         | Default | Alternative |
|--------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|-------------|
| `min_child_ab`           | Minimum proportion of alt-supporting reads relative to overall read depth for the proband                                                                       | 0.20    |             |
| `min_depth`              | Minimum read depth for calls. In order of preference/availability this is from `DP`, total of `AD`, or if both are absent, a dummy value which will always pass | 5       |             |
| `max_depth`              | Maximum read depth for calls. In order of preference/availability this is from `DP`, total of `AD`, or if both are absent, a dummy value which will always pass | 1000    |             |
| `min_gq`                 | Minumum genotype quality score for all variants considered                                                                                                      | 25      | 40          |
| `min_alt_depth`          | If all other conditions succeed, at least this many Alt observations to confirm a de Novo call                                                                  | 5       |             |

## Stage: `ValidateMOI`
  * Config Section: `ValidateMOI`
  * Description: Controls filtering parameters during per-family MOI testing

This part of the analysis defines 4 separate Filter models, and the one selected for each variant depends on the attributes and MOI pattern being tested:
1. `ClinVarDominant` - Applied when considering Dominant inheritance, if the variant has a P/LP ClinVar rating. Strict, but less strict than non-ClinVar P/LP Dominant

| Field                             | Purpose                                                        | Default |
|-----------------------------------|----------------------------------------------------------------|---------|
| `min_callset_ac_to_filter`        | Min. instances in the callset to apply intra-callset AF filter | 10      |
| `clinvar_dominant_gnomad_max_af`  | Max gnomAD 4.1 Combined AF                                     | 0.00005 |
| `clinvar_dominant_callset_max_af` | Max intra-callset AF                                           | 0.05    |

2. `ClinVar` - Applied to ClinVar P/LP variants when considering non-Dominant inheritance patterns.

| Field                      | Purpose                                                        | Default |
|----------------------------|----------------------------------------------------------------|---------|
| `min_callset_ac_to_filter` | Min. instances in the callset to apply intra-callset AF filter | 10      |
| `clinvar_gnomad_max_af`    | Max gnomAD 4.1 Combined AF                                     | 0.05    |
| `clinvar_callset_max_af`   | Max intra-callset AF                                           | 0.05    |

3. `Dominant` - Applied to all non-ClinVar Dominant MOIs

| Field                             | Purpose                                                        | Default |
|-----------------------------------|----------------------------------------------------------------|---------|
| `min_callset_ac_to_filter`        | Min. instances in the callset to apply intra-callset AF filter | 10      |
| `dominant_callset_max_ac`         | Max intra-callset AC                                           | 10      |
| `dominant_callset_max_af`         | Max intra-callset AF                                           | 0.01    |
| `dominant_callset_sv_max_af`      | Max intra-callset AC (Structural Variants)                     | 0.01    |
| `dominant_gnomad_max_af`          | Max gnomAD 4.1 Combined AF                                     | 0.00001 |
| `dominant_gnomad_max_ac`          | Max gnomAD 4.1 Combined AC                                     | 10      |
| `dominant_gnomad_max_homozygotes` | Max gnomAD 4.1 Combined Homozygotes                            | 0       |
| `dominant_gnomad_sv_max_af`       | Max gnomAD 4.1 Combined AF (Structural Variants)               | 0       |

4. `Global` - Applied to non-Dominant MOIs when the variant is not in ClinVar

| Field                      | Purpose                                                        | Default |
|----------------------------|----------------------------------------------------------------|---------|
| `min_callset_ac_to_filter` | Min. instances in the callset to apply intra-callset AF filter | 10      |
| `callset_max_af`           | Max intra-callset AF                                           | 0.01    |
| `callset_sv_max_af`        | Max intra-callset AF (Structural Variants)                     | 0.01    |
| `gnomad_max_af`            | Max gnomAD 4.1 Combined AF                                     | 0.01    |
| `gnomad_max_homozygotes`   | Max gnomAD 4.1 Combined Homozygotes                            | 5       |
| `gnomad_max_hemizygotes`   | Max gnomAD 4.1 Combined Hemizygotes                            | 5       |
| `gnomad_sv_max_af`         | Max gnomAD 4.1 Combined AF (Structural Variants)               | 5       |

Further filtering is done on all variant types and MOI, controlled by other parameters in the `ValidateMOI` block:

| Field            | Purpose                   | Default |
|------------------|---------------------------|---------|
| `min_alt_depth`  | Min. alt-supporting reads | 5       |
| `minimum_depth`  | Min. total read depth     | 10      |

And finally, there are some further meta-parameters which control the variants being considered. These parameters are used to facilitate inclusion of less confident annotations (e.g. _in silico_ predictions) whilst preventing noise from these tools dominating the overall results:

| Field                     | Purpose                                                                                                                                                      | Default          |
|---------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
| `ignore_categories`       | Categories which are ignored for this analysis. Variants with only these cats. assigned will never reach the report                                          | [exomiser, svdb] |
| `exomiser_rank_threshold` | If Exomiser is used, limit results to the top n hits                                                                                                         | 2                |
| `phenotype_match`         | List of categories, variants with only these assigned will be removed from the report unless a phenotype-match is detected                                   | [6]              |
| `support_categories`      | List of categories, variants with only these assigned will be removed from the report unless detected with a comp-het partner having a 'full' value category | [6]              |

## Stage: `HPOFlagging`
  * Config Section: `HPOFlagging`
  * Description: Controls the phenotypic matching between genes and HPO terms by Semsimian, used to draw attention to variants in the report.

| Field            | Purpose                                                                                                | Default |
|------------------|--------------------------------------------------------------------------------------------------------|---------|
| `semantic_match` | If True, run a semantic match between gene and proband HPO terms. If False, do a pure set intersection | True    |
| `min_similarity` | If `semantic_match`, this is the minumum required similarity when comparing termsets using Semsimian   | 14.0    |

## Stage: `CreateTalosHTML`
  * Config Section: `CreateTalosHTML`
  * Description: Controls how URLs linking out to external resources are generated and embedded into the report
  * Hyperlinks: The `CreateTalosHTML.hyperlinks` block can be absent, in which case no URL links are generated. If present, a `LinkEngine` will be generated to create all external resource links for the report

| Field                         | Purpose                                                                                                                                                                                                                                                                                                                |
|-------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `external_labels`             | Optional, path to an external labels JSON file. This is a nested dictionary of `sample_id[chr-pos-ref-alt]: [label, label2]`. If provided, the labels will be annotated on the report, useful in recording known false positives, or prior diagnoses etc. A log of all provided-but-not-in-results labels is generated |
| `hyperlinks.template`         | Mandatory if `hyperlinks` section is present. A String template containing an instance of `{sample}` for interpolating with the chosen Sample ID. This should define a template link to a per-family or per-proband resource                                                                                           |
| `hyperlinks.variant_template` | Similar to `template`, but containing both `{variant}` and `{sample}`. This is used to generate a link to a per-variant resource, embedding both Sample ID and `chr-pos-ref-alt` (seqr/gnomAD-style)                                                                                                                   |
| `hyperlinks.lookup`           | Optional, a JSON file mapping the VCF sample IDs to external IDs. If absent, the VCF IDs are embedded into any URLs generated                                                                                                                                                                                          |
| `hyperlinks.external`         | Boolean, if False (default) the VCF ID will be looked up in `lookup` (if present), or embedded directly into Template strings. If True, the sample `ext_id` (present in the phenopacket data) will be used as an index in the `lookup` dict or embedded into Strings if lookup is absent.                              |

This is designed to be as flexible as possible, but may be daunting if you're starting from scratch. Worked examples:

### Links out using the exact VCF ID
- VCF ID: `SAM1`
- Target URLs: `https://seqr.org.au/project/SAM1`
  - Template: "https://seqr.org.au/project/{sample}"
- Target per-variant URLs: `https://seqr.org.au/project/variants/1-12345-A-C/family/SAM1`
  - Per-variant template: "https://seqr.org.au/project/variants/{variant}/family/{sample}"

### Links out using a looked-up ID
- VCF ID: `SAM1`
- Link ID: `ExtSam1`
- Lookup JSON
```json
{
    "SAM1": "ExtSam1",
    "SAM2": "ExtSam2"
}
```
- Target URLs: `https://seqr.org.au/project/ExtSam1`
  - Template: "https://seqr.org.au/project/{sample}"
- Target per-variant URLs: `https://seqr.org.au/project/variants/1-12345-A-C/family/ExtSam1`
  - Per-variant template: "https://seqr.org.au/project/variants/{variant}/family/{sample}"

### Links out using an external ID
- VCF ID: `SAM1`
- Sample External ID in phenopackets/pedigree: `ExtSam1`
- `CreateTalosHTML.hyperlinks.external` = True
- Target URLs: `https://seqr.org.au/project/ExtSam1`
  - Template: "https://seqr.org.au/project/{sample}"
- Target per-variant URLs: `https://seqr.org.au/project/variants/1-12345-A-C/family/ExtSam1`
  - Per-variant template: "https://seqr.org.au/project/variants/{variant}/family/{sample}"
