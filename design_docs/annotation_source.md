# Annotation Fields

AIP uses VEP annotations when filtering and classifying variants. Much of this is done during the Hail runtime, with
some values used downstream when finalising inheritance pattern checks. This document breaks down the annotations used
at each point.

## Filtering

Variants with an AC/AN value above a given threshold are removed from consideration (Allele Count, Allele Number). The
threshold used is typically 1%
Variants with an assigned FILTER value (e.g. excessHet or VQSRTrancheXXx-YYy) are removed
Variants with a gnomAD (genomes or exomes) maf above the threshold are removed. Does not use sub-population AFs.
Variants not annotated against a GREEN gene (PanelApp gene list) are removed (set intersect with the geneIds annotation)
Transcript Consequences which are not either protein coding or on a MANE Select transcript are removed
Variants have to be well normalised (i.e. exactly 2 alleles present, 1 ref 1 alt, and alt != `*`)

Most of these filters can be overridden for ClinVar Pathogenic variants -

| Variable                                | Type          |
|-----------------------------------------|---------------|
| AC                                      | Int           |
| AN                                      | Int           |
| AF threshold                            | Float         |
| gnomad_exomes.AF                        | Float         |
| gnomad_genomes.AF                       | Float         |
| filters                                 | Array[String] |
| alleles                                 | Array[String] |
| geneIds                                 | String        |  # once exploded per-geneId, this is initially an Array
| vep.transcript_consequences.gene_id     | String        |
| vep.transcript_consequences.biotype     | String        |
| vep.transcript_consequences.mane_select | String        |

## ClinVar

The ClinVar filtering makes use of two fields - ClinVar significance (P/LP/VUS/U/LB/B), and number of 'gold stars'. At
runtime these values may be overwritten using the same datatype from a private source.

| Variable                      | Type   |
|-------------------------------|--------|
| clinvar.clinical_significance | String |
| clinvar.gold_stars            | Int    |

1. Variants are removed where the Significance contains Benign and Stars > 0
2. Variants are assigned the `clinvar_aip` flag where Significance contains `Pathogenic` and not `Conflicting`
3. Variants are assigned the `clinvar_aip_strong` flag where `clinvar_aip` is assigned and Stars > 0

Note: assignment of `clinvar_aip` is used to override some quality filters (e.g. these variants will be retained even if
they are common within the callset, common within gnomAD, FILTER'd by the variant caller or VQSR)

## Category 1 - ClinVar Pathogenic

Assignment of this Category indicates that a variant is seen in ClinVar with a pathogenic rating. The previously computed
field `clinvar_aip_strong` is used.

| Variable           | Type |
|--------------------|------|
| clinvar_aip_strong | Int  |

## Category 2 - New Disease Gene Association

Assignment of this category indicates that a variant has at least a moderate consequence, but the gene-disease assc. is
recent. This is done using a set of Gene IDs which were previously found to be new associations ([see here](PanelApp_interaction.md#application-of-new)).
Alongside the new disease gene association, a number of possible annotations are used to affirm relevance.

| Variable                                      | Type        |
|-----------------------------------------------|-------------|
| new_genes                                     | Set[String] |
| vep.transcript_consequences.consequence_terms | Set[String] |
| clinvar_aip                                   | Int         |
| CADD.PHRED                                    | Float       |
| dbnsfp.REVEL_score                            | Float       |

## Category 3 - High Impact variant

Combines severe predicted transcript consequence and LOFTEE to affirm if NMD is predicted as a result. If LOFTEE ranks
as unlikely to go through NMD this can be rescued using a ClinVar Pathogenic rating

| Variable                                      | Type        |
|-----------------------------------------------|-------------|
| vep.transcript_consequences.consequence_terms | Set[String] |
| clinvar_aip                                   | Int         |
| lof                                           | String      |

## Category 4 - De Novo Variants

Combines moderate variant impact with de novo inheritance checks, with numerous associated quality checks built in to the
hail de novo method (see [this](Hail_Filter_and_Label.md#category-4--de-novo-)). Note - prior to running the de novo
search the variant table is filtered down to at least moderate consequence (Missense or worse) OR High SpliceAI score.
The latter score is assigned in the next category (Category5); 5 is assigned before 4 is tested.

The Pedigree is also ingested for this category.

| Variable                                      | Type        |
|-----------------------------------------------|-------------|
| vep.transcript_consequences.consequence_terms | Set[String] |
| categoryboolean5                              | Int         |
| GQ                                            | Int         |
| GT                                            | hail.Call   | # this is a specific representation used in Hail, though any valid representation would work - 0/1/2 for HomRef/Het/HomAlt
| PL                                            | Array[Int]  |


## Category 5 - SpliceAI predictions

Solely focused on the annotations provided by SpliceAI, tested against a threshold (taken from config)

| Variable              | Type        |
|-----------------------|-------------|
| splice_ai.delta_score | Float       |
| threshold             | Float       |


## Support - In Silico Consequence Predictors

Aims to use a number of _in silico_ prediction scores to reach a form or consensus. Relevant thresholds for each
individual tool are taken from config.

| Variable                                   | Type  |
|--------------------------------------------|-------|
| CADD.PHRED                                 | Float |
| dbnsfp.REVEL_score                         | Float |
| dbnsfp.MutationTaster_pred                 | Float |
| vep.transcript_consequences.polyphen_score | Float |
| vep.transcript_consequences.sift_score     | Float |


## Further filtering - during MOI phase

After applying the Category labels, there are some situations where further annotation fields are used. Specifically
these are during the application of inheritance pattern checks.

Examples: when processing variants we use a stricter AF check if the variant is Dominant compared to Recessive.
If the variant is in gnomAD and the number of population HOM variants is above a specified threshold we would skip it. If
the variant is in a male on chrX we would skip the variant if there are an above-threshold number of Hemi instances seen
in gnomAD.

| Variable            | Type  |
|---------------------|-------|
| gnomad_exomes.AF    | Float |
| gnomad_exomes.AN    | Float |
| gnomad_exomes.AC    | Float |
| gnomad_exomes.Hom   | Float |
| gnomad_exomes.Hemi  | Float |
| gnomad_genomes.AF   | Float |
| gnomad_genomes.AN   | Float |
| gnomad_genomes.AC   | Float |
| gnomad_genomes.Hom  | Float |
| gnomad_genomes.Hemi | Float |
| DP                  | Float | # total ref + alt depths at the site should be above a threshold (config)
