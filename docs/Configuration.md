# Configuration File

Talos is controlled by a **single TOML** configuration file. A fully commented template lives at
[`src/talos/example_config.toml`](../src/talos/example_config.toml).
This document explains the configuration options, the default values, and sensible alternatives.

## Global (stage-agnostic) settings

| Config section | Purpose |
|----------------|---------|
| **`[categories]`** | Maps each short category code to a human-readable label shown in the HTML report. |
| **`[result_history]`** | Local or cloud path – a directory path that stores prior Talos result-history files (one JSON per run). If defined, at the start of each run Talos loads the most recent result file in that directory. Variants that already appeared in an earlier run retain their original date; variants seen for the first time in the current run get today’s date. Downstream you can simply filter on the current run-date column to review only the newly prioritised variants.

---

## Stages: `DownloadPanelApp`, `UnifiedPanelAppParser`
  * Config Section: `GeneratePanelData`
  * Description: Talos uses PanelApp as its source of disease-associated genes. Optionally, it can also use it as a source of disease-specific gene associations used for the analysis of disease-specific cohorts.
    These settings control which PanelApp deployment is used, which default Panel to use as the base for analysis, which additional genes/panels to include, and which genes to blacklist from the analysis.

| Field                 | Purpose                                                                                                                                                                                                                                                                                                                                                                                        | Default                                                        | Alternative                                                              |
|-----------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------|--------------------------------------------------------------------------|
| `panelapp`            | PanelApp server to query                                                                                                                                                                                                                                                                                                                                                                       | https://panelapp-aus.org/api/v1/panels                         | It is highly recommended to use PanelApp Aus due to enhanced annotations |
| `default_panel`       | Base panel for the analysis, applied to all participants                                                                                                                                                                                                                                                                                                                                       | "137" (The [Mendeliome](https://panelapp-aus.org/panels/137/)) |                                                                          |
| `require_pheno_match` | List of gene symbols to remove from analysis unless phenotype matched to a participant. This is useful to exclude genes that frequently result in noise.                                                                                                                                                                                                                                       | FLG, GJB2, F2, F5, HFE                                         |                                                                          |
| `forbidden_genes`     | List of genes that will be excluded from all results.                                                                                                                                                                                                                                                                                                                                          | N/A                                                            | Late onset/incidental findings                                           |
| `forced_panels`       | List of PanelApp panel IDs to apply to all cohort members. Typically used during the analysis of disease-specific cohorts. For example, if analysing a cohort consisting of patients who all have a suspected mitochondrial disorder, you could set this to `203` to highlight any genes associated with the Aus PanelApp [Mitochondrial disease panel](https://panelapp-aus.org/panels/203/). | N/A                                                            | Phenotype specific panels                                                |
| `within_x_months`     | Integer, genes added to panels within X months are flagged in the report and used to determine “new” genes for the "ClinVar Recent Gene" category                                                                                                                                                                                                                                                     | 24                                                             |                                                                          |
| `manual_overrides`    | A manually defined panel: gene symbols, MOI. Can be used to generate an artificial panel, e.g. for recently discovered genes which are not yet in PanelApp                                                                                                                                                                                                                                     | N/A                                                            | Recently discovered genes, not yet in PanelApp                           |

---

## Stage: `RunHailFiltering`
  * Config Section: `RunHailFiltering`
  * Description: Controls filtering thresholds and de novo variant detection parameters during the Hail Stage

| Field                     | Purpose                                                                                                                                                                                                                                    | Default                                                                       | Alternative                                                                                        |
|---------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| `csq_string`              | list of named fields to concatenate into the VCF export                                                                                                                                                                                    | see example config                                                            | reduce/extend based on the available annotations                                                   |
| `critical_csq`            | High impact Consequence terms                                                                                                                                                                                                              | frameshift, splice_acceptor, splice_donor, start_lost, stop_gained, stop_lost | if using VEP annotations, use altered VEP equivalents (e.g. "frameshift" -> "frameshift_variant")  |
| `additional_csq`          | During de novo variant filtering, we reduce the variant set to critical consequences + this list                                                                                                                                           | missense                                                                      |                                                                                                    |
| `af_semi_rare`            | A coarse global population frequency filter - variants more frequent than this are hard filtered during initial processing. Later in the process more rigorous filters are applied, so this is a 'mild' filter to retain plausible compound-het candidates     | 0.01                                                                          |                                                                                                    |
| `callset_af_sv_recessive` | A coarse population and callset frequency filter - if an SV variant is more common than this either within the callset, or in the gnomAD annotations, remove                                                                               | 0.03                                                                          |                                                                                                    |

  * Config Section: `RunHailFiltering.de_novo`
  * Description: Controls filtering thresholds and de novo variant detection parameters during the Hail Stage

| Field                | Purpose                                                                                                                                                         | Default | Alternative |
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|-------------|
| `min_child_ab`       | Minimum proportion of alt-supporting reads relative to overall read depth for the proband                                                                       | 0.20    |             |
| `min_depth`          | Minimum read depth for calls. In order of preference/availability this is from `DP`, total of `AD`, or if both are absent, a dummy value which will always pass | 5       |             |
| `max_depth`          | Maximum read depth for calls. In order of preference/availability this is from `DP`, total of `AD`, or if both are absent, a dummy value which will always pass | 1000    |             |
| `min_alt_depth`      | If all other conditions succeed, at least this many Alt observations to confirm a de Novo call                                                                  | 5       |             |
| `min_proband_gq`     | Minumum genotype quality score for all variants in probands/affected participants (detected via Pedigree)                                                       | 25      | 40          |
| `min_all_sample_gq`  | Minimum GQ applied to all samples                                                                                                                               | 15      |             |

---

## Stage: `ValidateMOI`
  * Config Section: `ValidateMOI`
  * Description: Controls all filtering parameters for mode-of-inheritance (MOI) validation during per-family analysis.

Talos applies four distinct MOI filter models, each tailored to different variant types and inheritance patterns. For each variant, Talos selects the most appropriate model based on both the variant’s attributes and the MOI being evaluated.

**Parameters shared by all models:**

| Field                             | Purpose                                                        | Default |
|-----------------------------------|----------------------------------------------------------------|---------|
| `min_callset_ac_to_filter`        | Minimum allele count (AC) in the callset before intra-callset allele frequency (AF) filtering is applied. Intra-callset AF filters are highly effective for removing callset-specific artefacts, but may be too stringent for small cohorts. This parameter prevents the exclusion of rare variants in smaller callsets by only applying AF filtering when a variant reaches the specified AC threshold.  | 10      |
| `min_alt_depth`  | Min. alt-supporting reads | 5       |
| `minimum_depth`  | Min. total read depth     | 10      |

**The Four MOI Filter models are:**

1. `ClinVarDominant` - Applied when considering Dominant inheritance, if the variant has a P/LP ClinVar rating. Strict, but less strict than non-ClinVar P/LP Dominant

| Field                             | Purpose                                                        | Default |
|-----------------------------------|----------------------------------------------------------------|---------|
| `clinvar_dominant_gnomad_max_af`  | Max gnomAD 4.1 Combined AF                                     | 0.00005 |
| `clinvar_dominant_callset_max_af` | Max intra-callset AF                                           | 0.05    |

2. `ClinVar` - Applied to ClinVar P/LP variants when considering non-Dominant inheritance patterns.

| Field                      | Purpose                                                        | Default |
|----------------------------|----------------------------------------------------------------|---------|
| `clinvar_gnomad_max_af`    | Max gnomAD 4.1 Combined AF                                     | 0.05    |
| `clinvar_callset_max_af`   | Max intra-callset AF                                           | 0.05    |

3. `Dominant` - Applied to all non-ClinVar Dominant MOIs

| Field                             | Purpose                                                        | Default |
|-----------------------------------|----------------------------------------------------------------|---------|
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
| `callset_max_af`           | Max intra-callset AF                                           | 0.01    |
| `callset_sv_max_af`        | Max intra-callset AF (Structural Variants)                     | 0.01    |
| `gnomad_max_af`            | Max gnomAD 4.1 Combined AF                                     | 0.01    |
| `gnomad_max_homozygotes`   | Max gnomAD 4.1 Combined Homozygotes                            | 5       |
| `gnomad_max_hemizygotes`   | Max gnomAD 4.1 Combined Hemizygotes                            | 5       |
| `gnomad_sv_max_af`         | Max gnomAD 4.1 Combined AF (Structural Variants)               | 5       |

And finally, there are some further meta-parameters which control the variants being considered. These parameters are used to facilitate inclusion of less confident annotations (e.g. _in silico_ predictions) whilst preventing noise from these tools dominating the overall results:

| Field                     | Purpose                                                                                                                                                      | Default          |
|---------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
| `ignore_categories`       | Categories which are ignored for this analysis. Variants with only these cats. assigned will never reach the report                                          | [exomiser, svdb] |
| `exomiser_rank_threshold` | If Exomiser is used, limit results to the top n hits                                                                                                         | 2                |
| `phenotype_match`         | List of categories, variants with only these assigned will be removed from the report unless a phenotype-match is detected                                   | [6]              |
| `support_categories`      | List of categories, variants with only these assigned will be removed from the report unless detected with a comp-het partner having a 'full' value category | [alphamissense, clinvar0star] |

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
| `hyperlinks.template`         | Mandatory if `hyperlinks` section is present. A String template containing an instance of `{sample}` for interpolating with the chosen Sample ID. This should define a template link to a per-family or per-proband resource                                                                                           |
| `hyperlinks.variant_template` | Similar to `template`, but containing both `{variant}` and `{sample}`. This is used to generate a link to a per-variant resource, embedding both Sample ID and `chr-pos-ref-alt` (seqr/gnomAD-style)                                                                                                                   |

### Hyperlinks

To generate Hyperlinks to Seqr, this section should contain the `hyperlinks` block with relevant templates, alongside a file mapping the VCF sample IDs to the IDs used in the target system. This should be added to the [NextFlow config](https://github.com/populationgenomics/talos/blob/main/docs/NextflowConfiguration.md#talosconfig) as `params.seqr_lookup`, and can be a TSV, CSV, or JSON file. This is designed to be as flexible as possible, but may be daunting if you're starting from scratch. Worked examples:

### Links out using the exact VCF ID

- Scenario: I want to generate links using a template, where the sample ID from the VCF is used in the hyperlink
- VCF ID: `SAM1`
- Target URLs: `https://seqr.org.au/project/SAM1`
  - Template: "https://seqr.org.au/project/{sample}"
- Target per-variant URLs: `https://seqr.org.au/project/variants/1-12345-A-C/family/SAM1`
  - Variant template: "https://seqr.org.au/project/variants/{variant}/family/{sample}"

### Links out using an alternate ID
- Scenario: I want to generate links using a template, but the sample ID from the VCF is not suitable for the target URL, so I want to map from the VCF sample ID to a URL-appropriate ID in a TSV/CSV/JSON file
- VCF ID: `SAM1`
- Link ID: `ExtSam1`
- Lookup JSON (file identified by `params.seqr_lookup` in NextFlow config)
```json
{
    "SAM1": "ExtSam1",
    "SAM2": "ExtSam2"
}
```
- Target URLs: `https://seqr.org.au/project/ExtSam1`
  - Template: "https://seqr.org.au/project/{sample}"
- Target per-variant URLs: `https://seqr.org.au/project/variants/1-12345-A-C/family/ExtSam1`
  - Variant template: "https://seqr.org.au/project/variants/{variant}/family/{sample}"
