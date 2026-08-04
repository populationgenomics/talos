# SV Annotation Module

A Nextflow workflow which annotates a joint-called Structural Variant VCF with:

- **gene consequences**, using [GATK `SVAnnotate`](https://gatk.broadinstitute.org/hc/en-us/articles/13832752531611-SVAnnotate) against the MANE GTF
- **population allele frequencies**, using [SVAFotate](https://github.com/fakedrtom/SVAFotate) against gnomAD-SV

This is a re-implementation of the same core steps we (CPG) use internally. Our internal usage centres around the
GATK-SV workflow, and our CPG-Flow wrapped implementation of it. The terminal stage of this workflow is [Annotation](https://github.com/populationgenomics/cpg-flow-gatk-sv/blob/main/src/cpg_flow_gatk_sv/multisample_workflow.py#L701),
which is done using GATK's SvAnnotate tool for consequence prediction, and a complex interval-overlap-resolution step
to match variants to gnomAD frequencies.

Instead of re-implementing the exact process here, I've split the annotation into two phases:

  - Consequence: handled using SVAnnotate, and exact replica of the GATK-SV process
  - Pop.Freq: handled using [SVAFotate](https://github.com/fakedrtom/SVAFotate)

These two steps, and pre-processing of relevant input files, are engaged only if an SV file is included in the input
TSV file, with the same core conceit as small variants and Mito data - a single joint-called VCF should contain the whole
group of Samples being processed, which should also match the Pedigree and Small-variant data.

A separate sub-workflow, `SV_ANNOTATION` has been created to handle these steps. `SV_ANNOTATION` publishes an annotated
VCF per cohort, `RunHailFilteringSv` filters and labels it with `CategoryBooleanSV1`, and `ValidateMOI` folds the result
into the report. **The code is complete and has been run against the real reference data.**
What remains is [open questions about real callsets](#whats-left), not missing modules.

Existed before this work: `talos.run_hail_filtering_sv`, but wired only into the CPG-internal `cpg-flow` path,
consuming a VCF annotated by a separate upstream GATK-SV run. This module makes the open-source SV path
self-contained.

## What's left

Ordered by how much they could cost.

### 1. A full Nextflow run of the SV path

Every process has been driven by hand against the real reference data, and `CreateSequenceDictionary` has been
run under Nextflow, but `nextflow run . --input_tsv nextflow/inputs/test_sv.tsv` has never been executed start
to finish. `talos.nf`'s new joins — the inner join onto `ch_sv_annotated`, and the `remainder: true` re-attach
of the `NO_SV` sentinel — are therefore unexercised by Nextflow itself.

### 2. Breakend coverage

BNDs never receive `PREDICTED_LOF`, so they are dropped before frequency filtering is ever consulted — see
[Breakends are dropped, not flooded](#breakends-are-dropped-not-flooded). This is a coverage gap, not a volume
risk.

**Needs deciding:** whether BND coverage is wanted, and if so what a "LoF breakend" means.
`--max-breakend-as-cnv-length` is the lever, but it changes which variants reach the report. Measure BND volume
in a real callset first — this may be a non-issue.

### 3. Overlap threshold and matching fidelity

`sv_overlap_fraction = 0.5` is a convention, not a tuned value, and SVAFotate cannot express per-SVTYPE or
per-size-band thresholds — see [Matching rules](#matching-rules). Small SVs are where reciprocal overlap is least stable.

**Needs deciding:** whether 0.5 holds across the size range. Requires a real callset to answer.

### 4. Deferred and low-priority

- **Non-coding annotations are inert.** `PREDICTED_NONCODING_BREAKPOINT` and `PREDICTED_NONCODING_SPAN` are
  produced and published, but nothing reads them. `DNase` is 98.9% of the BED, so it is too common to filter on
- **`MANE_Plus_Clinical` transcripts.** The unfiltered GTF works, but 72 genes carry two transcripts and `GNAS`
  carries three — over SVAnnotate's documented limit. If a gene annotation ever looks wrong, filtering to
  `tag "MANE_Select"` is the first thing to try
- **SVAFotate pinning.** No releases or tags upstream, so `docker/SVAFotate_Dockerfile` pins commit
  `30b5004a0f4d26959c6b9a82f165651585293626`, last verified 2026-07-30. Version bumps are a manual SHA change
  with no changelog
- **Uncovered by tests:** a DEL abutting a multi-megabase DEL, and any variant off chr1

## Traps

Things that look wrong and are not, or look fine and are not. Each has cost real debugging time.

### `Best_gnomAD_ID` is populated even when nothing matched

`-a best` reports the highest-ranked *candidate* overlap regardless of whether it cleared `-f`. A naive rename
would stamp a spurious gnomAD SVID onto nearly every variant, so **the rename is gated on `gnomAD_Count > 0`**.

This is not an edge case: two of the eight variants in the test fixture hit it. A 200 bp DEL returned
`Max_AF=0;gnomAD_Count=0` while still carrying `Best_gnomAD_ID` at OFP 0.003, and the PADI6 deletion did the
same at OFP 0.13.

### Use `Max_AF`, not `Best_gnomAD_AF`

`Max_AF` is the maximum across all qualifying matches — the conservative choice for rare disease filtering.
Under-filtering a common variant is a worse outcome than losing a rare one.

### `--max-breakend-as-cnv-length` cannot be passed its documented default

`gatk SVAnnotate --help` reports `Default value: -1`, but passing `-1` is rejected with `minimum allowed
value 0`. Omit the flag rather than passing its documented default.

## How it works

### Module structure

```text
nextflow/modules/annotation/
  AnnotateSvWithGatk/main.nf          # SVAnnotate -> gene consequences
  AnnotateSvWithSvafotate/main.nf     # SVAFotate -> population AFs
  RenameSvAfFields/main.nf            # Max_AF -> gnomad_v4.1_sv_AF
nextflow/modules/prep/
  CreateSequenceDictionary/main.nf    # ref.fa -> ref.dict
nextflow/modules/talos/
  RunHailFilteringSv/main.nf          # labels CategoryBooleanSV1
nextflow/sv_annotation.nf             # SV_ANNOTATION workflow
nextflow/assets/NO_SV                 # sentinel for cohorts with no SV data
src/talos/annotation_scripts/
  rename_sv_af_fields.py              # invoked by RenameSvAfFields
```

The chain is `AnnotateSvWithGatk` -> `AnnotateSvWithSvafotate` -> `RenameSvAfFields`, mirroring how
`nextflow/annotation.nf` chains its SNV steps. `RunHailFilteringSv` lives in `talos.nf`, not here, because it
needs the PanelApp JSON, MANE JSON and pedigree that workflow already assembles.

This is a separate workflow rather than a branch inside `ANNOTATION` because that workflow's normalise,
region-filter, split and echtvar steps are all SNV-specific, SV VCFs are small enough that `vcf_split_n`
sharding is unnecessary, and the SV input is a single joint-called VCF with no equivalent of the
`shards`/`ss_vcf_dir`/`vcf` branching.

`CreateSequenceDictionary` lives under `prep/` but is invoked from `sv_annotation.nf`, so the SNV-only path
never pulls the GATK image.

The rename is a Python script rather than `bcftools annotate --rename-annots` because it reads
`RunHailFilteringSv.gnomad_population` from the Talos config — the same key `run_hail_filtering_sv` reads, so
the written and expected field names cannot drift — and because the `gnomAD_Count` gate is conditional logic a
flat rename cannot express.

### Input and gating

The SV VCF is supplied per cohort as an **`sv` column in the input TSV**, following the same sentinel
convention as `mito`: a cohort with no SV data gets `nextflow/assets/NO_SV`.

**The input VCF must be indexed.** `sv_annotation.nf` picks up `{sv}.tbi` with `checkIfExists: true`, so an
un-indexed input fails while channels are being built, before any process runs. The `complete` branch requires
the same of the previously published VCF.

`main.nf` reads the TSV header eagerly and only calls `SV_ANNOTATION` when an `sv` column exists, because the
workflow body exits if the MANE GTF or SVAFotate BED are missing — an SNV-only user must never be made to
download 478 MB of SV reference data.

Cohorts are then branched three ways:

| Branch | Condition | Behaviour |
|--------|-----------|-----------|
| `sentinel` | the SV path is `NO_SV` | dropped, no annotation attempted |
| `complete` | `{outputDir}/{cohort}_outputs/{cohort}_sv_annotated.vcf.bgz` exists | the existing VCF is emitted as-is |
| `pending` | otherwise | run the full three-step chain |

The gate reads back exactly where the previous run published, so both expensive steps are skipped on a re-run.
`CreateSequenceDictionary` is driven from the `pending` branch rather than `ch_ref_genome` directly — `first()`
over an empty channel emits nothing, so with no pending work it never runs.

### Output contract

`rearrange_annotations()` reads each of these **unconditionally**, so a missing field raises rather than
degrading:

| INFO field | Notes |
|------------|-------|
| `PREDICTED_LOF` | Array of gene symbols. Rows with an empty array are filtered out outright |
| `gnomad_v4.1_sv_AF`, `gnomad_v4.1_sv_SVID` | Prefix from `config_retrieve(['RunHailFilteringSv', 'gnomad_population'], 'gnomad_v4.1')` |
| `SVTYPE`, `SVLEN`, `END` | Standard SV description fields |
| `AC`, `AF`, `AN` | Joint-call frequencies |
| `N_HET`, `N_HOMALT` | Joint-call genotype counts |
| `MALE_AF`, `FEMALE_AF` | **Array-typed.** `AF_MALE`/`AF_FEMALE` accepted as an alternative spelling |

`MALE_AF`/`FEMALE_AF` deserve emphasis: `rearrange_annotations()` branches on `AF_MALE` and otherwise reads
`mt.info.MALE_AF` directly, so a VCF carrying neither raises. They must be arrays because
`filter_matrix_by_ac()` subscripts them as `mt.info.male_af[0]`, and that is the only filter a variant with a
null gnomAD AF is subject to. It is *not* the last line of defence for BNDs: the empty-`PREDICTED_LOF` filter
runs before `rearrange_annotations()`, so BNDs are already gone — see
[Breakends are dropped, not flooded](#breakends-are-dropped-not-flooded).

Defaulted if absent: `ALGORITHMS` (to `['gCNV']`), and `STATUS`/`CHR2`/`END2` (to a null string).

#### The chain supplies none of the joint-call fields

Of the above, this workflow produces only `PREDICTED_LOF`, `gnomad_v4.1_sv_AF` and `gnomad_v4.1_sv_SVID`.
`SVTYPE`, `SVLEN` and `END` are passed through from the input. `AC`, `AF`, `AN`, `N_HET`, `N_HOMALT`, `MALE_AF`
and `FEMALE_AF` are **never written by any step here** — see [what's left, item 3](#3-whether-real-input-vcfs-satisfy-the-output-contract).

### Matching rules

`-f/--minf` defaults to **0.001** — effectively "any single-basepair overlap matches" — so it must be set
explicitly. `params.sv_overlap_fraction` exists for this, defaulting to `0.5`.

What SVAFotate can express is narrower than ideal:

- matching is reciprocal overlap plus an **identical** `SVTYPE` string. No DUP↔CNV or DEL↔CNV cross-matching.
  For gnomAD rows this costs little — 533 `CNV` records in a 3M-row sample, against ~1.04M `DEL` and ~226k `DUP`
- `-f` takes one value **per source**, not per SVTYPE or size band. Tightening above 1 Mb, or falling back to
  breakpoint proximity below ~1 kb, is not expressible
- insertions match on near-identical coordinates regardless of `SVLEN` unless `--ins` is passed, which expands
  `END` by `SVLEN`. Inserted sequence content is never compared
- `-l/--lim` caps reference SV size, but only with `--cov` or `--uniq`

The invocation — note `-b` takes the [pre-filtered BED](#the-bed-is-pre-filtered-to-gnomad-rows), not the
shipped one:

```bash
svafotate annotate -v "${vcf}" -o "${cohort}_sv_popaf.vcf" -b new_af_file.bed.gz \
    -s gnomAD -f "${params.sv_overlap_fraction}" -a best --ins --cpu "${task.cpus}"
```

It writes an uncompressed VCF — the SVAFotate image has no htslib CLI, so `RenameSvAfFields` handles bgzip and
tabix in the Talos container.

### The BED is pre-filtered to gnomAD rows

`AnnotateSvWithSvafotate` does not hand SVAFotate the shipped BED. It first makes a single streaming pass —
`zcat | awk -F'\t' 'NR == 1 || $6 == "gnomAD"' | gzip` — and annotates against the result, to hold down the
footprint described in [what's left, item 2](#2-resource-sizing-and-svafotates-memory-floor). Rows are
selected; columns and contig names are untouched.

Four details of that one line have already cost debugging time:

- **do not reach for `zcat ... | head -1` to grab the header.** `head` exits after one line, `zcat` takes
  SIGPIPE on the rest of a multi-GB stream, and `set -o pipefail` turns that into exit 141
- **`NR == 1 ||` is load-bearing.** The header's own column names contain the string `gnomAD`, so a bare
  `grep gnomAD` emits it twice — once as a header, once as a data row
- **`$6` must be escaped as `\$6`.** Nextflow script blocks are Groovy GStrings, so a bare `$6` interpolates
  away and awk dies on `NR == 1 ||  == "gnomAD"`
- **`-F'\t'` is set explicitly** because awk's default separator collapses runs of whitespace, so a future BED
  carrying an empty field before `SOURCE` would shift the column numbering silently

The pass re-runs per cohort, per run, and it is a pure function of a static reference file — so it is a
candidate for the `CreateSequenceDictionary` treatment: computed once into `processed_annotations` and gated on
existence. Not done, because a ~30% trim that leaves the process at 20 GB may not be the right shape at all;
worth revisiting once the real memory floor is known.

## Reference data

All plain downloads, already in `large_files/gather_files.sh`. The only persistent derived file is the sequence
dictionary; the [gnomAD-only BED](#the-bed-is-pre-filtered-to-gnomad-rows) is rebuilt inside the task and never
leaves the work directory.

| Resource | Filename in `large_files/` | Size | Consumed by |
|----------|---------------------------|------|-------------|
| SVAFotate popAF BED | `SVAFotate_SV_popAFs.GRCh38.v4.1.bed.gz` | 469 MB (325 MB once filtered) | `AnnotateSvWithSvafotate` |
| MANE Ensembl GTF | `MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` | 8.6 MB | `AnnotateSvWithGatk` |
| Non-coding elements BED | `noncoding.sort.hg38.bed` | 64 MB | `AnnotateSvWithGatk` |
| Reference FASTA | `ref.fa` | 3.0 GB | `CreateSequenceDictionary` |
| Sequence dictionary | `ref.dict` (in `processed_annotations`) | 48 KB | `AnnotateSvWithGatk` |

Sources:

- SVAFotate BED — Zenodo record [11642574](https://zenodo.org/records/11642574), file
  `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`, saved locally without the `core_` infix
- GTF — `https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`
- Non-coding BED — `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed`

### The SVAFotate BED

Native GRCh38 gnomAD v4.1, no liftover, merging four callsets distinguished by a `SOURCE` column. A superset of
the GATK-SV frequency TSV, carrying `Het`, `HomAlt` and `PopMax_AF` columns that file lacks. Mandatory column
order: `CHROM`, `START`, `END`, `SVLEN`, `SVTYPE`, `SOURCE`, `SV_ID`, `AF`, then 169 optional columns.

Those eight are not actually sufficient in practice — `svafotate annotate` always writes
`Max_Het`/`Max_HomAlt`/`Max_PopMax_AF` and raises `KeyError` if `Het`, `HomAlt` and `PopMax_AF` are absent.

`START` is 0-based and VCF `POS` is 1-based, so a query matching a record exactly is `POS = START + 1`,
`END = END`.

### The MANE GTF

`SVAnnotate` wants primary or canonical transcripts, documented as "1-2 transcripts per gene only". NCBI's
`MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` satisfies this as downloaded — and despite the filename its contigs
are UCSC-style, matching `ref.fa`.

It is the same MANE release as `MANE.GRCh38.v1.5.summary.txt.gz`, which Talos already consumes. That matters:
`read_and_filter_mane_json()` maps `PREDICTED_LOF` symbols to ENSG IDs using the summary file, so sharing a
release guarantees the two gene vocabularies agree. GATK-SV's `gencode.v39.canonical.gtf.gz` would introduce a
second vocabulary and cause silent mapping losses.

Two properties look like they need filtering and do not: it carries `MANE_Plus_Clinical` alongside
`MANE_Select` so 72 genes have two transcripts and `GNAS` has three (a deletion spanning `GNAS` still annotates
correctly), and 35 of its 79 contigs are `_alt`/`_fix` scaffolds absent from `ref.fa` (`SVAnnotate` does not
object, even with `--sequence-dictionary`).

### The non-coding BED

Adds `PREDICTED_NONCODING_BREAKPOINT` and `PREDICTED_NONCODING_SPAN`, neither of which Talos reads — carried
for provenance and future use.

Its format contradicts SVAnnotate's own help, which describes six columns and a header; the shipped file has
**four columns and no header** and is read without complaint. Element classes are heavily unbalanced — `DNase`
is 2107358 of 2131313 rows (98.9%), with `Enhancer`, `Tommerup_TADanno` and `HAR` making up the rest. `DNase`
is therefore near-universal and near-worthless as a discriminator.

### The sequence dictionary

`large_files/` ships `ref.fa` and `ref.fa.fai` but no `ref.dict`, which `SVAnnotate` needs for contig ordering.
`CreateSequenceDictionary` builds it in **8 seconds** against the real 3.0 GB reference — an earlier revision
assumed this would be slow under `linux/amd64` emulation and used a hand-built dictionary instead. It is not.

!!! note "This is the pipeline's only `publishDir`"
    The dictionary belongs to the reference genome rather than to a run's results, so `params.ref_dict` points
    into `processed_annotations` — and Nextflow's workflow output block cannot write outside `outputDir`. This
    one process therefore uses `publishDir params.processed_annotations, mode: 'copy'`. Once published, the
    existence check skips it on every subsequent run.

    `publishDir` is documented as superseded by the output block, but the two coexist without complaint —
    verified on Nextflow 26.04.4.

## Containers

The main Talos image needs **no changes**. GATK is already published upstream (and bundling Java 17 would be
required otherwise); SVAFotate's dependencies are irreconcilable with Talos's.

| Container                    | Provenance                                                     |
|------------------------------|----------------------------------------------------------------|
| `params.container`           | existing Talos image                                           |
| `params.gatk_container`      | `broadinstitute/gatk:4.6.2.0`, pulled from the public registry |
| `params.svafotate_container` | built locally from `docker/SVAFotate_Dockerfile`               |

SVAFotate pins `cyvcf2==0.30.4` (Talos needs `>=0.30.18`) and `pandas==1.2.3`/`numpy==1.22.3` (Talos needs
newer, per `hail~=0.2.137`). Those pins are not reproduced in the image — current releases work, which avoids
being held to a Python 3.9 base. The sole constraint is [`pandas<3`](#pandas3-is-required-in-the-svafotate-image).

!!! note "The GATK image is amd64 only"
    No arm64 build, so on Apple Silicon it runs under emulation with a platform warning. It works, but
    developer-machine timings will not represent production.

### Process resources

In the `process { withName: ... }` block in `nextflow.config`:

| Process                     | Memory | CPUs |
|-----------------------------|--------|------|
| `CreateSequenceDictionary`  | 4 GB   | 1    |
| `AnnotateSvWithGatk`        | 4 GB   | 1    |
| `AnnotateSvWithSvafotate`   | 16 GB  | 4    |
| `RenameSvAfFields`          | 2 GB   | 1    |

Only the SVAFotate figure is informed by measurement, and only as a floor — see
[what's left, item 2](#2-resource-sizing-and-svafotates-memory-floor).

## Configuration

`src/talos/example_config.toml` carries:

```toml
[RunHailFilteringSv]
gnomad_population = 'gnomad_v4.1'
```

This key was already read by `run_hail_filtering_sv.py` but absent from the example config, so it defaulted
silently. `rename_sv_af_fields.py` reads the same key, which is the point.

Nextflow params — `svafotate_bed`, `mane_gtf`, `svannotate_noncoding_bed`, `ref_dict`, `sv_overlap_fraction`,
`gatk_container`, `svafotate_container` — are documented in
[NextflowConfiguration.md](NextflowConfiguration.md). `gnomad_sv_freq` and `sv_input_vcf` have both been
removed.

## Testing

| File                                     | Needs                 | Covers                                                                                                                                                                           |
|------------------------------------------|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `test/test_hail_filtering_sv.py`         | Hail only             | `rearrange_annotations()`, both blanket filters, the hemizygous fix, both `MALE_AF`/`AF_MALE` spellings                                                                          |
| `test/test_rename_sv_af_fields.py`       | nothing               | the `gnomAD_Count` gate, `Max_AF` over `Best_gnomAD_AF`, the all-zero warning                                                                                                    |
| `test/test_sv_annotation_integration.py` | `svafotate` on `PATH` | behaviours only the real tool exhibits, including the [prefixed-BED regression guard](#do-not-rewrite-the-svafotate-beds-contigs). Skipped when absent, so unit CI is unaffected |

To run the integration tests, put a shim on `PATH` that shells out to the container:

```bash
printf '#!/bin/sh\nexec docker run --rm -v "$PWD":"$PWD" -w "$(pwd)" svafotate:0.1.0 svafotate "$@"\n' > /tmp/bin/svafotate
chmod +x /tmp/bin/svafotate
PATH=/tmp/bin:$PATH pytest test/containerised_test_sv_annotation_integration.py --basetemp=.pytest-tmp
```

`--basetemp` inside the repo matters — pytest's default `tmp_path` is under `/var/folders` on macOS, which
Docker Desktop does not share.

### Workflow test data

`nextflow/inputs/generate_sv_test_data.py` writes an eight-variant `joint_sv.vcf.bgz`, and
`nextflow/inputs/test_sv.tsv` is the standard test input plus an `sv` column. `test.tsv` itself deliberately has
no `sv` column, so the default test run never pulls the SV reference data.

Every case is coordinate-accurate against the real reference files rather than invented. Two properties to
preserve if it is regenerated:

- **the joint-call frequencies describe a notional 2000-allele callset, not the three samples in the VCF.**
  Derived from a trio, every variant would sit above `callset_af_sv_recessive = 0.03` and `filter_matrix_by_ac`
  would empty the callset
- **exactly one variant is expected to survive.** If two survive, or none, that is the signal — not the absence
  of a crash

## What has been verified

The chain has been run process by process against the **real** reference data — the real MANE GTF, the real
non-coding BED, a `ref.dict` built from the real `ref.fa`, and a chr1 slice of the real SVAFotate BED — using
the eight-variant fixture, and on through `RunHailFilteringSv`:

| Variant                 | Outcome                                                                                                                                                     |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `lof_del_padi6`         | `PREDICTED_LOF=PADI6`, no qualifying gnomAD match, rare in callset — **the only survivor**, labelled `categorybooleansv1=1` with `lof_ensg=ENSG00000276747` |
| `common_del`            | `PREDICTED_LOF=OR4F16`, `gnomad_v4.1_sv_AF=0.2238` — dropped by `filter_matrix_by_af`                                                                       |
| `inv_common_in_callset` | `PREDICTED_LOF=TNFRSF18`, gnomAD AF 8e-06 but callset AF 0.2 — dropped by `filter_matrix_by_ac`                                                             |
| `small_del_in_large`    | `Max_AF=0`, `gnomAD_Count=0`, spurious `Best_gnomAD_ID` at OFP 0.003 withheld by the rename                                                                 |
| `dup_two_matches`       | four overlapping gnomAD DUPs, `Max_AF` took the highest (0.704181)                                                                                          |
| `intragenic_dup_gja9`   | `PREDICTED_INTRAGENIC_EXON_DUP=GJA9`, no `PREDICTED_LOF`, dropped                                                                                           |
| `ins_end_eq_pos`        | matched `gnomAD-SV_v3_INS_chr1_3f94b1dc` at OFP 1 via `--ins`, `PREDICTED_INTRONIC`, dropped                                                                |
| `bnd_1`                 | `PREDICTED_BREAKEND_EXONIC` and `PREDICTED_INTRONIC` but no `PREDICTED_LOF`, dropped                                                                        |

Earlier, under Nextflow with a harness mirroring `main.nf`: the annotating run chains all three processes and
publishes to `{cohort}_outputs/`; the gated run executes zero processes (`completed=0, cached=0`) with a clean
work directory, which rules out Nextflow's own caching; a `NO_SV` cohort is dropped and
`CreateSequenceDictionary` does not run; and re-publishing an in-place file leaves it byte-identical with a
usable index.

`CreateSequenceDictionary` has since been run under Nextflow against the real reference: it publishes to
`processed_annotations`, the gate skips it on a re-run with a wiped work directory, and the output matches a
hand-built dictionary.

## Risks

### Breakends are dropped, not flooded

An earlier revision predicted BNDs would evade frequency filtering and inflate output volume. The opposite is
true, and the mechanism is upstream of the AF filter entirely: `SVAnnotate` does not assign `PREDICTED_LOF` to a
BND.

| BND placement                                     | Annotation produced                                       |
|---------------------------------------------------|-----------------------------------------------------------|
| breakends flanking `PADI6`                        | `PREDICTED_INTERGENIC;PREDICTED_NEAREST_TSS=PADI6`        |
| both breakends inside the gene body               | `PREDICTED_INTRONIC=PADI6,RCC2`                           |
| one breakend intronic, one exonic in another gene | `PREDICTED_INTRONIC=PADI6;PREDICTED_BREAKEND_EXONIC=GJA9` |

Since rows with an empty `PREDICTED_LOF` are dropped outright, BNDs never reach the frequency filter, and the
`MISSING_INT = 0` behaviour in `filter_matrix_by_af()` — where a null AF passes — is never exercised for them.

Where a caller emits both a symbolic `<INV>`/`<DEL>` record and its constituent BND records, resolve to the
symbolic record before annotating.

### Input VCF conformance

`SVAnnotate` is strict — it needs `SVTYPE`, correct `END`/`SVLEN`, and `CHR2`/`END2` for breakends. It now runs
correctly against real reference data, so the risk has moved from "will the tool work" to "will real input
satisfy it".

This is a *separate* question from
[item 3](#3-whether-real-input-vcfs-satisfy-the-output-contract). `SVAnnotate` needs the SV description fields
to annotate at all; `run_hail_filtering_sv` additionally needs the joint-call frequency and genotype-count
fields, which nothing here adds. A gCNV-native VCF could satisfy the first and fail the second, at the very end
of the chain.

## Alternatives considered

| Tool | Verdict |
|------|---------|
| GATK-SV `AnnotateExternalAF` | Emits exactly the `{prefix}_AF`/`{prefix}_SVID` names Talos expects, so needs no rename step. Internally just `svtk vcf2bed` plus `bedtools intersect -r -f`. Rejected in favour of SVAFotate, whose BED is a superset of the GATK frequency TSV — but a reasonable fallback if the SVAFotate image proves hard to maintain |
| AnnotSV | Comprehensive clinical annotation, including gnomAD-SV, DGV, DDD and OMIM. Rejected because it emits TSV rather than VCF, which does not fit the pipeline shape |
| Truvari `collapse` | Good size- and sequence-aware matching, but awkward to repurpose for external AF transfer |
| `vcfanno` | Rejected. Plain interval overlap with no reciprocal-overlap or SVTYPE constraint, so it would annotate a 200 bp DEL with a 3 Mb gnomAD DEL's frequency |
