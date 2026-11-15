# Talos

Variant prioritisation for rare disease. Talos annotates cohorts' variant callsets once, then reanalyses
them on a regular cadence against current ClinVar and PanelApp data, reporting candidate diagnoses per
family.

## Language

### Runs and outputs

**Cohort**:
The named unit of input that Talos analyses: one row of the input TSV, with its callset, pedigree and
config. Outside CPG the name has no stricter meaning. At CPG, it's a Metamist cohort.

**Annotation products**:
A cohort's annotated shards, shard manifest and annotated SV VCF. They're produced once per callset and
reused by every later reanalysis cycle.
_Avoid_: annotations (the prepared reference data are also annotations)

**Pre-annotated SV VCF**:
An SV callset whose population frequencies and gene consequences were already added upstream (e.g. by
GATK-SV), so Talos's own SV annotation is skipped for it.
_Avoid_: annotated SV VCF (that's Talos's own SV annotation product)

**Analysis results**:
The outputs of one reanalysis cycle for one cohort: the results JSON, HTML report, labelled VCFs and
panel data. Each cycle's results are kept separately, by date.

**Reanalysis cycle**:
One run of Talos over a cohort's existing annotation products, using that month's prepared data. A
cycle's results JSON becomes the next cycle's history.
_Avoid_: rerun

**History**:
A cohort's results JSON from an earlier reanalysis cycle. It's used to tell new findings from ones that
have already been reported.

**Next-cycle input TSV**:
A proposed input TSV for the following reanalysis cycle. It's this run's input, with each cohort's history
pointed at the results this run produced.

### Inputs and evidence

**Callset**:
The variant calls for all of a cohort's samples, as supplied to Talos. A cohort can have up to four: small
variants (required), structural variants, mitochondrial variants and short tandem repeats.

**Shard**:
One genomic slice of a cohort's small-variant callset. Small variants are annotated and labelled a shard at
a time; every other callset is handled whole.
_Avoid_: chunk, fragment, interval

**Pedigree**:
The family structure, sex and affected status of every participant in a cohort, optionally with their
phenotype terms.

**Participant**:
One individual in the pedigree. Talos only reports variants for affected participants.
_Avoid_: sample (when you mean the person rather than their data), patient

**Prepared data**:
The reference resources Talos derives from external sources ahead of any cohort: ClinVar decisions, PanelApp
content, transcript and population-frequency data. ClinVar and PanelApp are refreshed monthly.
_Avoid_: annotations

**ClinVar decision**:
A single pathogenicity call for a variant, re-derived from ClinVar's individual submissions by
ClinvArbitration rather than taken from ClinVar's own summary.

### Genes and panels

**Panel**:
A PanelApp gene list for a disease area. Panels are the only source of the genes, modes of inheritance and
phenotype tags that Talos works with.

**Default panel**:
The panel applied to every participant: the base set of established Mendelian disease genes.
_Avoid_: base panel, Mendeliome

**Forced panel**:
A panel applied to every participant in a cohort by configuration, regardless of their phenotypes.

**Matched panel**:
A panel applied to one participant because their phenotype terms match the panel's.

**Green gene**:
A gene with PanelApp's highest confidence rating on at least one applied panel. By default Talos
only analyses Green genes.

**New gene**:
A gene that became Green within a configured recency window, so that findings in it are recent knowledge.
_Avoid_: recent gene (except in the category name ClinVar Recent Gene), novel gene (Talos doesn't do gene
discovery)

**Panel data**:
A cohort's gene and panel context for one reanalysis cycle: which panels apply to each participant, and each
gene's mode of inheritance, panels and new-gene status.

### Variants and categories

**Category**:
A named reason that a variant is of interest, such as ClinVar P/LP, High Impact or De Novo. A variant with
no category is never considered for reporting.
_Avoid_: logic module, label, flag, class

**Support category**:
A category that is too weak to report alone. A variant with only support categories is reported solely as
the second hit to a variant with a full category.
_Avoid_: supporting module, secondary category

**Labelled VCF**:
A callset reduced to the variants with at least one category, each carrying its categories.
_Avoid_: filtered VCF, categorised VCF

**Mode of inheritance (MOI)**:
The inheritance pattern that PanelApp records for a gene. It decides which genotype and family
configurations count as a possible diagnosis.

**MOI test**:
The check that a categorised variant's genotypes across a family are consistent with one mode of
inheritance, at an acceptable population frequency.
_Avoid_: MOI filter, inheritance check

**Compound het**:
Two different heterozygous variants in the same gene in one participant, not shown to sit on the same
parental copy.

**Support variant**:
The partner of a reported variant in a compound het.
_Avoid_: second hit (as a noun for the variant), partner variant

**Reported variant**:
One variant in one participant that has a category and passed an MOI test. The same variant in two
participants is two reported variants.
_Avoid_: candidate, hit, result

**Phenotype match**:
Evidence that a reported variant's gene fits the participant's phenotype terms, either through a matched
or forced panel, or through similarity between the participant's and the gene's phenotype terms.

**Flag**:
A caution attached to a reported variant about the quality or ambiguity of the evidence. It never changes
whether the variant is reported.
_Avoid_: warning, label

### Reanalysis dates

**First tagged**:
The date a reported variant first appeared for a participant, or the date its gene became Green if that
was later.
_Avoid_: first seen

**Evidence last updated**:
The most recent date that the case for a reported variant got stronger: a new category, a higher ClinVar
star rating, a first phenotype match, or its gene becoming Green.

### CPG

**Metamist**:
CPG's sample-metadata service, and the system of record for datasets, cohorts and sequencing groups.

**Dataset**:
A Metamist project, named without any `-test` suffix. It owns the buckets that its cohorts' outputs are
written to.
_Avoid_: project (except when you specifically mean the Metamist project name, which may carry `-test`)

**Metamist cohort**:
A curated collection of sequencing groups within one dataset, identified by a `COH` ID. It decides which
samples are analysed. Its dataset decides where the outputs go.

**Sequencing group (SG)**:
All of the sequencing data for one sample, of one sequencing type and technology.
_Avoid_: sample (a sample can have several SGs)

**Access level**:
Whether a run works in a dataset's production (`full`, `standard`) or test (`test`) namespace. Every run
has exactly one.
_Avoid_: namespace, environment

**Dataset bucket**:
`gs://cpg-<dataset>-main`, or `gs://cpg-<dataset>-test` at the test access level. Talos outputs are laid
out under it by a fixed convention.

**Sequencing type**:
`genome` (WGS) or `exome` (WES). Every run handles exactly one.
_Avoid_: WGS/WES in parameters and paths
