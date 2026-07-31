# Plan: SV Annotation Module

A Nextflow annotation workflow which annotates a joint-called Structural Variant VCF with:

- **gene consequences**, using [GATK `SVAnnotate`](https://gatk.broadinstitute.org/hc/en-us/articles/13832752531611-SVAnnotate)
- **population allele frequencies**, using [SVAFotate](https://github.com/fakedrtom/SVAFotate) against gnomAD-SV

**Current state: the annotation chain is implemented and verified. Nothing consumes its output yet.**
`SV_ANNOTATION` runs end to end and publishes an annotated VCF per cohort, but `RunHailFilteringSv` does not
exist, so no SV variant reaches a Talos report. See [Implementation status](#implementation-status) and
[Unresolved points](#unresolved-points).

!!! note
    Flags and output field names in this document have been verified by running both tools: `svafotate 0.1.0`
    in the locally built container against a slice of the real reference BED, and `SVAnnotate` from
    `broadinstitute/gatk:4.6.2.0` against the real MANE GTF. Everything describing the Talos repository has
    been verified against the code.

    Two caveats on that verification. The reference BED used was a chr1:0-2Mb slice, not the full 469 MB file,
    and the input VCFs were hand-written fixtures rather than real joint-called output — so nothing here
    exercises the [joint-call INFO fields](#the-chain-supplies-none-of-the-joint-call-fields) a real callset
    would carry.

## Motivation

`talos.run_hail_filtering_sv` already exists and applies the `CategoryBooleanSV1` label. However, it is only
wired into the CPG-internal `cpg-flow` path (`src/talos/cpg_internal_scripts/talos_stages.py`), where it
consumes an SV VCF that was annotated by a separate, upstream GATK-SV run.

The open-source Nextflow pipeline has no equivalent:

- `nextflow/modules/annotation/` contains no SV process
- `nextflow/modules/talos/` contains no `RunHailFilteringSv` process

This module makes the SV path self-contained, so that a raw joint-called SV VCF can be taken through to a
Talos report without depending on an external annotation pipeline.

## Output contract

The module's output must satisfy `talos.run_hail_filtering_sv`. The required INFO fields are derived from
`rearrange_annotations()` and the body of `main()` in `src/talos/run_hail_filtering_sv.py`.

### Hard requirements

`rearrange_annotations()` reads each of these unconditionally, so a missing field raises rather than
degrading:

| INFO field | Notes |
|------------|-------|
| `PREDICTED_LOF` | Array of gene symbols. Rows with an empty array are filtered out outright |
| `gnomad_v4.1_sv_AF` | Prefix is configurable — see [Configuration](#configuration) |
| `gnomad_v4.1_sv_SVID` | Prefix is configurable — see [Configuration](#configuration) |
| `SVTYPE`, `SVLEN`, `END` | Standard SV description fields |
| `AC`, `AF`, `AN` | Joint-call frequencies |
| `N_HET`, `N_HOMALT` | Joint-call genotype counts |
| `MALE_AF`, `FEMALE_AF` | **Array-typed.** `AF_MALE`/`AF_FEMALE` are accepted as an alternative spelling |

The `gnomad_v4.1` prefix comes from `config_retrieve(['RunHailFilteringSv', 'gnomad_population'], 'gnomad_v4.1')`.

`MALE_AF`/`FEMALE_AF` deserve emphasis, because an earlier revision of this plan listed them as optional and
defaulted. They are not. `rearrange_annotations()` branches on whether `AF_MALE` is present and otherwise
reads `mt.info.MALE_AF` directly, so a VCF carrying neither spelling raises. They must also be arrays —
`filter_matrix_by_ac()` subscripts them as `mt.info.male_af[0]`, and this is the only filter BNDs and
frequency-null variants are subject to.

### Soft requirements

Defaulted by `rearrange_annotations()` if absent:

| INFO field | Default behaviour |
|------------|-------------------|
| `ALGORITHMS` | Defaults to `['gCNV']` |
| `STATUS`, `CHR2`, `END2` | Inserted as a null string |

### The chain supplies none of the joint-call fields

Worth stating plainly, because it is easy to miss: of the hard requirements above, this workflow produces only
`PREDICTED_LOF`, `gnomad_v4.1_sv_AF` and `gnomad_v4.1_sv_SVID`. `SVTYPE`, `SVLEN` and `END` are passed through
from the input, and `AC`, `AF`, `AN`, `N_HET`, `N_HOMALT`, `MALE_AF` and `FEMALE_AF` are **never written by any
step here**. They have to already be present in the joint-called input VCF.

A full GATK-SV callset carries all of them. A bare gCNV or caller-native VCF may not, and would fail at
`rearrange_annotations()` rather than at annotation time — after both expensive steps have run. This is
untested, and is the single largest open risk. See [Unresolved points](#unresolved-points).

### The field naming mismatch

SVAFotate does **not** emit the field names Talos expects, so a third process is required to rename them.

By default `svafotate annotate` writes `Max_AF`, `Max_Het`, `Max_HomAlt`, `Max_PopMax_AF` and a
`[SOURCE]_Count` field per requested source. Passing `-a best` additionally emits `Best_[source]_ID`,
`Best_[source]_AF`, `Best_[source]_Het` and `Best_[source]_OFP` — where OFP is the Overlap Fraction Product
used to rank candidate matches.

The two renames required are therefore:

| SVAFotate field | Talos field |
|-----------------|-------------|
| `Max_AF` | `gnomad_v4.1_sv_AF` |
| `Best_gnomAD_ID` | `gnomad_v4.1_sv_SVID` |

Use `Max_AF`, **not** `Best_gnomAD_AF`. `Max_AF` is the maximum frequency across all qualifying matches,
which is the conservative choice for rare disease filtering — under-filtering a common variant is a worse
outcome than losing a rare one.

!!! danger "`Best_gnomAD_ID` is populated even when nothing matched"
    `-a best` reports the highest-ranked *candidate* overlap regardless of whether it cleared the `-f`
    threshold. A verified example: a 200 bp DEL with no qualifying match returned
    `Max_AF=0;gnomAD_Count=0` but still carried `Best_gnomAD_ID=gnomAD-SV_v3_DEL_chr1_f7b56087` with
    `Best_gnomAD_OFP=0.178508`, well below the `-f 0.5` cutoff.

    A naive rename would therefore stamp a spurious gnomAD SVID onto nearly every variant in the callset.
    **The rename must be gated on `gnomAD_Count > 0`**, emitting a null SVID otherwise.

## Module structure

Following the existing one-directory-per-process convention, with `set -euo pipefail`:

```text
nextflow/modules/annotation/
  AnnotateSvWithGatk/main.nf          # SVAnnotate -> gene consequences
  AnnotateSvWithSvafotate/main.nf     # SVAFotate -> population AFs
  RenameSvAfFields/main.nf            # Max_AF -> gnomad_v4.1_sv_AF
nextflow/modules/prep/
  CreateSequenceDictionary/main.nf    # ref.fa -> ref.dict
nextflow/modules/talos/
  RunHailFilteringSv/main.nf          # DOES NOT EXIST YET
nextflow/sv_annotation.nf             # SV_ANNOTATION workflow
nextflow/assets/NO_SV                 # sentinel for cohorts with no SV data
src/talos/annotation_scripts/
  rename_sv_af_fields.py              # invoked by RenameSvAfFields
```

Everything above exists except `RunHailFilteringSv`.

`CreateSequenceDictionary` lives under `prep/` but is invoked from `sv_annotation.nf`, not `preparation.nf`.
Keeping it out of the shared prep workflow means the SNV-only path never has to pull the GATK image.

The rename is a Python script rather than `bcftools annotate --rename-annots` because it has to read
`RunHailFilteringSv.gnomad_population` from the Talos config — the same key `run_hail_filtering_sv` reads — so
the written and expected field names cannot drift apart. It also needs conditional logic that a flat rename
cannot express, per the danger note above.

Each annotation process declares its own container: `params.gatk_container` for `AnnotateSvWithGatk`,
`params.svafotate_container` for `AnnotateSvWithSvafotate`, and the default `params.container` for
`RenameSvAfFields` and `RunHailFilteringSv`.

The annotation chain is `AnnotateSvWithGatk` -> `AnnotateSvWithSvafotate` -> `RenameSvAfFields`, mirroring
the way `nextflow/annotation.nf` chains `NormaliseAndRegionFilterVcf` -> `AnnotateWithEchtvar` ->
`AnnotateCsqWithBcftools`.

!!! note "No BED or GTF prep modules are required"
    Earlier revisions of this plan specified `BuildSvafotateBed` and `MakeSvAnnotateGtf` prep modules. Both
    have been dropped — the SVAFotate reference BED and the MANE GTF are both downloaded ready to use. See
    [Reference data](#reference-data).

### Why a separate workflow

This is a new `sv_annotation.nf` workflow rather than a branch inside `ANNOTATION`:

- the existing workflow's normalise, region-filter, split and echtvar steps are all SNV-specific, and would
  all need bypassing
- SV VCFs are small enough that the `vcf_split_n` sharding logic is unnecessary
- the SV input is a single joint-called VCF, so the `shards` / `ss_vcf_dir` / `vcf` input branching in
  `annotation.nf` has no SV equivalent — the input channel is simply `tuple(cohort, sv_vcf, config)`

### Input and gating

The SV VCF is supplied per cohort as an **`sv` column in the input TSV**, alongside the existing `mito`
column, and follows the same sentinel convention: a cohort with no SV data gets
`nextflow/assets/NO_SV`, which `main.nf` substitutes when the column is absent for that row.

`main.nf` reads the TSV header eagerly and only calls `SV_ANNOTATION` at all when an `sv` column exists. That
matters because the workflow body checks for the MANE GTF and SVAFotate BED and exits if they are missing —
an SNV-only user must never be made to download 478 MB of SV reference data.

Within the workflow, cohorts are branched three ways:

| Branch | Condition | Behaviour |
|--------|-----------|-----------|
| `sentinel` | the SV path is `NO_SV` | dropped, no annotation attempted |
| `complete` | `{outputDir}/{cohort}_outputs/{cohort}_sv_annotated.vcf.bgz` exists | the existing VCF is emitted as-is |
| `pending` | otherwise | run the full three-step chain |

The annotated VCF is published to `{cohort}_outputs/` by `main.nf`, matching the `path { id, ... ->
"${id}_outputs" }` convention used for `mts`, `labelled` and the rest. The gate therefore reads back exactly
where the previous run wrote, and both SVAnnotate and SVAFotate are skipped entirely on a re-run.

!!! note "The sequence dictionary is gated too"
    `CreateSequenceDictionary` is driven by a channel derived from the `pending` branch, not by
    `ch_ref_genome` directly. Without that it would rebuild the dictionary from the 3 GB reference even on a
    run where every cohort was already annotated — which it did, until the gate was tested. `first()` over an
    empty channel emits nothing, so with no pending work the process receives no input and never runs.

Publishing a file that is already in place is an established pattern here — `preparation.nf` does the same
thing with its ClinVar and PanelApp outputs when the current month's files already exist.

## Reference data

All three files are plain downloads, already added to `large_files/gather_files.sh` via the existing
`start_download` helper. No derived annotation files are needed beyond the sequence dictionary.

| Resource | Filename in `large_files/` | Size | Consumed by | Status |
|----------|---------------------------|------|-------------|--------|
| SVAFotate popAF BED | `SVAFotate_SV_popAFs.GRCh38.v4.1.bed.gz` | 469 MB | `AnnotateSvWithSvafotate` | downloaded |
| MANE Ensembl GTF | `MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` | 8.6 MB | `AnnotateSvWithGatk` | downloaded |
| Reference FASTA | `ref.fa` | 3.0 GB | `AnnotateSvWithGatk` | already present |
| Sequence dictionary | `ref.dict` (in `processed_annotations`) | small | `AnnotateSvWithGatk` | generated on first SV run |
| Non-coding elements BED | `noncoding.sort.hg38.bed` | 64 MB | `AnnotateSvWithGatk` | downloaded |

Sources:

- SVAFotate BED — Zenodo record [11642574](https://zenodo.org/records/11642574), file
  `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`, saved locally without the `core_` infix
- GTF — `https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`
- Non-coding BED — `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed`

!!! note "`gnomad_v4_SV.Freq.tsv.gz` was removed"
    `gather_files.sh` previously downloaded the GATK-SV team's `gnomad_v4_SV.Freq.tsv.gz` (83 MB), with a
    matching `gnomad_sv_freq` param. Both dated from an earlier revision in which a DIY `bedtools`
    implementation was preferred over SVAFotate, and the SVAFotate BED is a strict superset of that file. The
    download and the param have both been dropped; the 83 MB file may still be present in `large_files/` from
    an earlier run and can be deleted.

### The shipped BED is native gnomAD v4.1

An earlier revision of this plan claimed SVAFotate's published resources were gnomAD-SV v2.1 lifted over from
GRCh37. That was wrong. The file is `SVAFotate_SV_popAFs.GRCh38.v4.1.bed.gz` — native GRCh38 gnomAD v4.1, no
liftover — and it merges four callsets, distinguished by its `SOURCE` column: `gnomAD`, `CCDG`, `TOPMed` and
`ThousG`.

It is also a superset of the GATK-SV frequency TSV, carrying `Het`, `HomAlt` and `PopMax_AF` columns the
GATK file lacks. Mandatory column layout, in order: `CHROM`, `START`, `END`, `SVLEN`, `SVTYPE`, `SOURCE`,
`SV_ID`, `AF`, followed by 169 optional columns.

### `-s gnomAD` is mandatory

With no `-s`, SVAFotate uses every source in the BED and `Max_AF` becomes the maximum across all four
callsets. Writing that into a field named `gnomad_v4.1_sv_AF` would be mislabelling the data, and the
downstream Talos filter would attribute CCDG or TOPMed frequencies to gnomAD.

Restricting to `-s gnomAD` has a useful side effect. The BED's SVTYPE vocabulary is
`BND CNV CPX CTX DEL DUP INS INV MEI`, and SVAFotate requires an **identical** SVTYPE string to match — so a
separate `MEI` type would silently never match a query `INS`. Across a 3M-row sample, gnomAD-sourced rows use
only `DEL BND INS DUP CPX INV CNV CTX`, with no `MEI` at all: gnomAD calls mobile elements as `INS`. Scoping
to gnomAD removes that mismatch class entirely.

### Do not rewrite the BED's contigs

The BED uses Ensembl-style contig names (`1`, `2`, …) with no `chr` prefix, while Talos is UCSC-style
(`chr1`) throughout — the reference is UCSC hg38 and `run_hail_filtering_sv` runs Hail on GRCh38.

This looks like a silent-failure trap, and it is one — but **not** in the direction that intuition suggests.
SVAFotate normalises the *query VCF* by stripping the prefix (`svafotate_main.py:631`,
`v.CHROM[3:] if v.CHROM.startswith("chr")`). The source AF BED supplied via `-b` is never normalised. The
shipped Ensembl-style BED is therefore already correct, and "fixing" it would break annotation.

Verified by running one `chr1` test VCF against two slices of the real BED:

| Source BED contigs | `DEL_match` | `DUP_match` | `INS_match` |
|--------------------|-------------|-------------|-------------|
| `1` — as shipped | `Max_AF=0.11884` | `Max_AF=0.196779` | `Max_AF=0.053239` |
| `chr1` — rewritten | `Max_AF=0` | `Max_AF=0` | `Max_AF=0` |

So: use the BED exactly as downloaded, add no prefix-rewriting prep step, and pass the query VCF through
unmodified. The failure mode is entirely silent — a prefixed BED yields `Max_AF=0` for every variant with no
error raised, making the whole callset look rare — so this needs an explicit regression test.

The one exception is the optional `--target` BED, which *is* prefix-stripped on read
(`svafotate_main.py:590`), and so tolerates either style.

### The MANE GTF is used directly, but must be decompressed

`SVAnnotate` requires a GTF containing primary or canonical transcripts, documented as "1-2 transcripts per
gene only". NCBI's `MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` satisfies this without a prep step — despite the
filename, its contigs are UCSC-style (`chr1`), matching `ref.fa`.

The one handling requirement is decompression. `SVAnnotate` has no codec for a gzipped GTF and fails with
`Cannot read ... because no suitable codecs found`. `AnnotateSvWithGatk` therefore runs `gzip -dc` into a
scratch file before invoking the tool. Uncompressed, the file is read with `EnsemblGtfCodec`.

It is also the same MANE release as `MANE.GRCh38.v1.5.summary.txt.gz`, which Talos already consumes. That
matters because `run_hail_filtering_sv.read_and_filter_mane_json()` maps the gene symbols in `PREDICTED_LOF`
to ENSG IDs using the summary file — sharing a release guarantees the two gene vocabularies agree. Using
GATK-SV's `gencode.v39.canonical.gtf.gz` instead would introduce a second vocabulary and cause silent
mapping losses at that step.

Two properties of the file look like they should require filtering, and were tested. Neither does:

- **It is not strictly one transcript per gene.** The GTF carries `MANE_Plus_Clinical` transcripts alongside
  `MANE_Select`, so of 19437 transcripts, 72 genes have two and one — `GNAS` — has three, exceeding the
  documented limit. Filtering to `tag "MANE_Select"` would give an exact 1:1 across 19363 genes. It is not
  necessary: a deletion spanning `GNAS` annotates correctly as `PREDICTED_LOF=GNAS` with the unfiltered file.
- **35 of its 79 contigs are absent from `ref.fa`**, all `_alt` and `_fix` patch scaffolds. `SVAnnotate` does
  not object, even with `--sequence-dictionary` supplied.

Verified end to end: a deletion spanning `PADI6` yields `PREDICTED_LOF=PADI6`, and a duplication inside
`GJA9` yields `PREDICTED_INTRAGENIC_EXON_DUP=GJA9`. Symbols are emitted bare, as
`read_and_filter_mane_json()` expects.

### The non-coding BED is wired in, but nothing reads it

`--non-coding-bed` takes the GATK-SV public resource `noncoding.sort.hg38.bed` (64 MB, 2131313 rows, UCSC
contigs matching `ref.fa`). It adds two INFO fields: `PREDICTED_NONCODING_BREAKPOINT` and
`PREDICTED_NONCODING_SPAN`.

Its format contradicts SVAnnotate's own documentation, which is worth recording so nobody "fixes" it. The help
describes "BED file (with header). Columns: chrom, start, end, name, score (.), strand" — six columns and a
header. The shipped file has **four columns and no header**. `SVAnnotate` reads it without complaint.

The element classes are heavily unbalanced:

| Class | Rows | Share |
|-------|------|-------|
| `DNase` | 2107358 | 98.9% |
| `Enhancer` | 18209 | 0.9% |
| `Tommerup_TADanno` | 3012 | 0.1% |
| `HAR` | 2734 | 0.1% |

Because `DNase` is almost the entire file, `PREDICTED_NONCODING_BREAKPOINT=DNase` is close to universal — both
variants in a two-variant fixture picked it up. As a discriminator it is near-worthless; `Enhancer`, `HAR` and
`Tommerup_TADanno` are the informative classes.

!!! note "This changes no Talos output today"
    `PREDICTED_LOF` is byte-identical with and without the flag, and `run_hail_filtering_sv` reads no
    `PREDICTED_NONCODING_*` field. The annotation is carried for provenance and future use only. Verified
    surviving the full chain — both fields are present in the final published VCF, so SVAFotate and the rename
    step preserve them.

### The sequence dictionary is generated

`large_files/` contains `ref.fa` and `ref.fa.fai`, but no `ref.dict`, which `SVAnnotate` needs for contig
ordering. `CreateSequenceDictionary` builds it with `gatk CreateSequenceDictionary` into
`${params.processed_annotations}`, guarded on both the file not already existing and there being pending SV
work — see [Input and gating](#input-and-gating).

!!! warning "Not yet run against the real reference"
    Verification used a `.dict` hand-built from `ref.fa.fai` with `awk`, to avoid a slow pass over the 3.0 GB
    FASTA under `linux/amd64` emulation. The `CreateSequenceDictionary` process itself has therefore not been
    observed producing a dictionary from the real reference.

## Containers

The main Talos image needs **no changes**. Neither tool is added to it: GATK because the official image is
already published and Java 17 would otherwise have to be bundled, and SVAFotate because its dependencies are
irreconcilable with Talos's.

| Container | Provenance |
|-----------|------------|
| `params.container` | existing Talos image |
| `params.gatk_container` | `broadinstitute/gatk:4.6.2.0`, pulled from the public registry |
| `params.svafotate_container` | built locally from `docker/SVAFotate_Dockerfile` |

### Why SVAFotate needs its own image

SVAFotate's `requirements.txt` pins versions that directly conflict with Talos:

| SVAFotate pins | Talos requires |
|---|---|
| `cyvcf2==0.30.4` | `cyvcf2>=0.30.18` |
| `pandas==1.2.3`, `numpy==1.22.3` | newer, per `hail~=0.2.137` |

Those upstream pins are not reproduced in the image. Current releases of `cyvcf2`, `numpy`, `pyranges` and
`ncls` work, which avoids being held to a Python 3.9 base — `pandas 1.2.3` predates CPython 3.10 support and
has no wheels for it.

### `pandas<3` is a required pin

The single exception is `pandas`, which **must** be constrained to `<3`.

pandas 3.x returns the new `str` dtype for string columns. `pyranges 0.1.4` predates it: its `null_types()`
helper accepts only `object` or a dtype whose name contains `"string"`, and raises otherwise. Every
`svafotate annotate` run against pandas 3.0.5 therefore died with:

```text
Exception: Unknown dtype str in a column SVTYPE
```

Pinning `pandas<3` resolves to 2.3.3, which annotates cleanly alongside numpy 2.4.6. Note the pin must be
quoted in the `Dockerfile` (`"pandas<3"`) — unquoted, `<3` is a shell redirect from a file named `3`.

!!! warning "`RUN svafotate --version` is not a sufficient smoke test"
    The `--version` check passed on the broken image, because the failure is inside `pyranges.join` at
    annotation time rather than at import. A meaningful build-time test needs a few-line BED and VCF baked in
    and an actual `svafotate annotate` invocation.

### Process resources

Entries are in the `process { withName: ... }` block in `nextflow.config`:

| Process | Memory | CPUs |
|---------|--------|------|
| `CreateSequenceDictionary` | 4 GB | 1 |
| `AnnotateSvWithGatk` | 4 GB | 1 |
| `AnnotateSvWithSvafotate` | 16 GB | 4 |
| `RenameSvAfFields` | 2 GB | 1 |

These are starting guesses, not measurements. SVAFotate holds the reference BED in memory and the shipped file
is 469 MB compressed, hence the 16 GB; its `--cpu` argument is wired to `task.cpus`. All four need tuning
against a real callset — nothing here has been run on anything larger than a four-variant fixture.

!!! note "The GATK image is amd64 only"
    `broadinstitute/gatk:4.6.2.0` has no arm64 build, so on Apple Silicon it runs under emulation with a
    platform warning. It works, but is slow enough that timings taken on a developer machine will not
    represent production.

## Configuration

### Talos config

`src/talos/example_config.toml` now carries:

```toml
[RunHailFilteringSv]
gnomad_population = 'gnomad_v4.1'
```

This key was already read by `run_hail_filtering_sv.py` but absent from the example config, so it defaulted
silently. `rename_sv_af_fields.py` reads the same key, which is the point — the field names written and the
field names expected are derived from one source.

### Nextflow params

| Parameter | Description |
|-----------|-------------|
| `svafotate_bed` | Path in `large_files` to the downloaded SVAFotate popAF BED |
| `mane_gtf` | Path in `large_files` to the MANE Ensembl GTF |
| `svannotate_noncoding_bed` | Path in `large_files` to the GATK-SV non-coding elements BED |
| `ref_dict` | Path in `processed_annotations` to the generated sequence dictionary |
| `sv_overlap_fraction` | Reciprocal overlap threshold for AF matching. Default `0.5` |
| `gatk_container` | GATK image tag, `broadinstitute/gatk:4.6.2.0` |
| `svafotate_container` | Locally built SVAFotate image tag, `svafotate:0.1.0` |

`gnomad_sv_freq` has been removed, along with its download, per [Reference data](#reference-data).
`sv_input_vcf` has also been removed — the SV VCF is a per-cohort TSV column, not a global param.

New parameters should also be documented in `docs/NextflowConfiguration.md`.

## Matching rules

`-f/--minf` defaults to **0.001**, which is effectively "any single-basepair overlap matches" and would
annotate a 200 bp DEL with a 3 Mb gnomAD DEL's frequency. It **must** be set explicitly — `params.sv_overlap_fraction`
exists for this, defaulting to `0.5`.

What SVAFotate can actually express is narrower than is ideal:

- matching is reciprocal-overlap plus an **identical** `SVTYPE` string. There is no DUP↔CNV or DEL↔CNV
  cross-matching. For gnomAD-sourced rows this costs little — only 533 `CNV` records appeared in a 3M-row
  sample, against ~1.04M `DEL` and ~226k `DUP`
- `-f` takes a space-separated list, but one value **per source**, not per SVTYPE or per size band. Tightening
  the threshold above 1 Mb, or falling back to breakpoint proximity below ~1 kb, is not expressible
- insertions match on near-identical coordinates regardless of `SVLEN` unless `--ins` is passed, which expands
  `END` by `SVLEN` to allow reciprocal-overlap matching. Inserted sequence content is never compared
- `-l/--lim` caps the size of reference SVs considered, but only when `--cov` or `--uniq` is in use

The invocation `AnnotateSvWithSvafotate` uses is therefore:

```bash
svafotate annotate \
    -v "${vcf}" \
    -o "${cohort}_sv_popaf.vcf" \
    -b "${svafotate_bed}" \
    -s gnomAD \
    -f "${params.sv_overlap_fraction}" \
    -a best \
    --ins \
    --cpu "${task.cpus}"
```

It writes an uncompressed VCF: the SVAFotate image carries no htslib CLI, so compression and indexing happen
in `RenameSvAfFields`, which runs in the Talos container.

## Implementation status

### Done

| Component | Notes |
|-----------|-------|
| Containers | SVAFotate built locally from `docker/SVAFotate_Dockerfile`; GATK pulled. Talos image unchanged |
| `CreateSequenceDictionary` | Gated on both file absence and pending work |
| `AnnotateSvWithGatk` | Decompresses the GTF, then `SVAnnotate` with `--non-coding-bed` |
| `AnnotateSvWithSvafotate` | `-s gnomAD -f 0.5 -a best --ins` |
| `RenameSvAfFields` | Wraps `rename_sv_af_fields.py`, then bgzip + tabix |
| `sv_annotation.nf` | `SV_ANNOTATION`, with the three-way per-cohort branch |
| `main.nf` wiring | `sv` TSV column, `NO_SV` sentinel, publish to `{cohort}_outputs/` |
| `nextflow.config` | Params and `withName` resources; `sv_input_vcf` and `gnomad_sv_freq` removed |
| `example_config.toml` | `[RunHailFilteringSv]` block |
| `gather_files.sh` | SV downloads consolidated, redundant TSV removed |
| `NextflowConfiguration.md` | SV params documented, with the input, re-run and BED-contig notes |
| `mkdocs.yml` | This page added under Reference. `mkdocs build --strict` passes |

### Verified

The full chain has been run under Nextflow, using a harness mirroring `main.nf`'s publish and output blocks, a
four-variant SV VCF, and a chr1:0-2Mb slice of the reference BED:

- **annotating run** — all three processes execute and chain. The final VCF carries `PREDICTED_LOF=PADI6` from
  SVAnnotate and `gnomad_v4.1_sv_AF` from the rename on the same record, with both new INFO header lines
  declared. Published files land in `{cohort}_outputs/`
- **gated run** — with a clean work directory and no `-resume`, zero processes execute (`completed=0,
  cached=0`) and the previously published VCF is emitted. The empty work directory rules out Nextflow's own
  caching, so this is the gate
- **no-SV-data run** — a cohort on the `NO_SV` sentinel is dropped, and `CreateSequenceDictionary` does not run
- **re-publishing is non-destructive** — after the gated run republishes the in-place file to its own path, the
  bgzip file is byte-identical, still valid, and its tabix index still usable
- **rename edge cases**, by running the script directly against real SVAFotate output: the `gnomAD_Count` gate
  withholding a spurious SVID, the all-zero-AF warning firing, and the config-driven field prefix

### Not done

| Item | Notes |
|------|-------|
| `RunHailFilteringSv` | Absent. Needs PanelApp JSON, MANE JSON and pedigree, which `talos.nf` already assembles — so it belongs there, not in `sv_annotation.nf` |
| SV test data | No `generate_sv_test_data.py`. `nextflow/inputs/test.tsv` deliberately has no `sv` column, so the repo test workflow is unchanged |
| `docker/SVAFotate_Dockerfile` header | Comments still claim a Python 3.9 base and the upstream `pandas==1.2.3`/`numpy==1.22.3` pins. Neither is true; the image is `python:3.11-slim-bullseye` with only `pandas<3` constrained |
| Build-time smoke test | `RUN svafotate --version` does not catch the pandas failure mode |
| Stale download | `large_files/gnomad_v4_SV.Freq.tsv.gz` (83 MB) may remain on disk from before its removal |

## Unresolved points

Ordered by how much they could cost.

### 1. Nothing consumes the annotated VCF

`SV_ANNOTATION` publishes a VCF and stops. Until `RunHailFilteringSv` exists and `talos.nf` calls it, no SV
reaches a report, and the output contract above is unverified against the actual consumer. Everything below is
downstream of fixing this.

### 2. Whether real input VCFs satisfy the output contract

Per [the note above](#the-chain-supplies-none-of-the-joint-call-fields), `AC`, `AF`, `AN`, `N_HET`, `N_HOMALT`,
`MALE_AF` and `FEMALE_AF` must be in the input and are written by nothing here. `MALE_AF`/`FEMALE_AF` must also
be array-typed.

**Needs deciding:** what actually produces the cohort SV VCFs. If it is full GATK-SV, this is likely fine. If
it is gCNV or a caller-native VCF, a `NormaliseSvVcf` step is required, and the failure currently surfaces late
— at `rearrange_annotations()`, after both expensive annotation steps have run. A cheap upfront conformance
check on the input VCF header would convert that into a fast, clear failure.

### 3. Breakend coverage

BNDs never receive `PREDICTED_LOF` and are therefore dropped before frequency filtering — see
[Breakends are dropped, not flooded](#breakends-are-dropped-not-flooded). `--max-breakend-as-cnv-length` would
let them acquire consequences, but changes which variants reach the report.

**Needs deciding:** whether BND coverage is wanted at all, and if so what a "LoF breakend" means. Measure BND
volume in a real callset first — this may be a non-issue in practice.

### 4. Overlap threshold and matching fidelity

`sv_overlap_fraction = 0.5` is a convention, not a tuned value, and SVAFotate cannot express per-SVTYPE or
per-size-band thresholds — see [Matching rules](#matching-rules). Small SVs are where reciprocal overlap is
least stable.

**Needs deciding:** whether 0.5 is acceptable across the size range, or whether small variants need separate
handling. Requires a real callset to answer.

### 5. Resource sizing

Every `withName` value is a guess, taken against a four-variant fixture. The 16 GB for SVAFotate is inferred
from the BED being 469 MB compressed, not from observed usage.

### 6. Deferred and low-priority

- **Non-coding annotations are inert.** `PREDICTED_NONCODING_BREAKPOINT` and `PREDICTED_NONCODING_SPAN` are now
  produced and published, but `run_hail_filtering_sv` reads neither, so they affect nothing. Making use of them
  would mean extending the filtering script — and `DNase`, at 98.9% of the BED, is too common to filter on
- **`MANE_Plus_Clinical` transcripts.** The unfiltered GTF works, but 73 genes have more than one transcript
  and GNAS has three, over SVAnnotate's documented limit. If a gene-level annotation ever looks wrong, filtering
  to `tag "MANE_Select"` is the first thing to try
- **SVAFotate pinning.** No releases or tags upstream, so version bumps are a manual SHA change with no
  changelog

## Testing guidance

`nextflow/inputs/` already contains `generate_test_data.py` and `generate_mito_test_data.py` as patterns for a
new `generate_sv_test_data.py`. Make the SV test data adversarial rather than nominal:

- **a prefixed-BED regression test.** Assert that the shipped Ensembl-style BED annotates a `chr`-prefixed
  query VCF with non-zero AFs. This is the single most important test — it is the only guard against someone
  "normalising" the BED and silently zeroing every frequency in the callset
- a variant whose only candidate overlap falls below `-f`, asserting that `gnomad_v4.1_sv_SVID` is null and
  not the spurious `Best_gnomAD_ID`
- a 200 bp DEL abutting a 3 Mb DEL, catching a mis-set `-f`
- an INS where `END == POS`
- a DUP matching two overlapping gnomAD records with different AFs, asserting `Max_AF` takes the higher
- at least one each of DEL, DUP, INV, INS and BND

## Risks

### Breakends are dropped, not flooded

An earlier revision of this plan predicted that BNDs would evade frequency filtering and inflate output
volume. Testing shows the opposite, and the mechanism is upstream of the AF filter entirely.

`SVAnnotate` does not assign `PREDICTED_LOF` to a BND. Two cases were tested against the real GTF:

| BND placement | Annotation produced |
|---------------|---------------------|
| breakends flanking `PADI6` | `PREDICTED_INTERGENIC;PREDICTED_NEAREST_TSS=PADI6` |
| both breakends inside the gene body | `PREDICTED_INTRONIC=PADI6,RCC2` |

Neither yields `PREDICTED_LOF`. Since `run_hail_filtering_sv` drops rows with an empty `PREDICTED_LOF` array
outright, BNDs are discarded before the frequency filter is ever consulted. The `MISSING_INT = 0` behaviour in
`filter_matrix_by_af()` — where a null AF passes — is therefore not reached for breakends.

The consequence is the reverse of the original concern: BNDs contribute nothing to Talos SV output, whether or
not SVAFotate can match them. That is a coverage gap rather than a volume risk, and closing it would mean
deciding what a "LoF breakend" should mean, not tuning a frequency threshold.

`--max-breakend-as-cnv-length` is the obvious lever — it reinterprets sub-threshold BNDs as deletions or
duplications, which would let them acquire consequence annotations. It is not enabled here, because turning it
on changes which variants reach the report and that should be a deliberate, measured choice.

!!! note "A GATK documentation inconsistency"
    `gatk SVAnnotate --help` reports `--max-breakend-as-cnv-length` as `Default value: -1`, but passing `-1`
    explicitly is rejected with `minimum allowed value 0`. Omit the flag rather than passing its documented
    default.

Where a caller emits both a symbolic `<INV>`/`<DEL>` record and its constituent BND records, resolve to the
symbolic record before annotating.

### Input VCF conformance

`SVAnnotate` is strict about input VCF conformance — it needs `SVTYPE`, correct `END`/`SVLEN`, and
`CHR2`/`END2` for breakends. It was previously flagged as the main schedule risk; it now runs correctly against
hand-written fixtures, so the risk has moved from "will the tool work" to "will real input satisfy it".

Note this is a *separate* conformance question from the one in
[Unresolved points](#2-whether-real-input-vcfs-satisfy-the-output-contract). `SVAnnotate` needs the SV
description fields to annotate at all; `run_hail_filtering_sv` additionally needs the joint-call frequency and
genotype-count fields, which nothing in this workflow adds. A gCNV-native VCF could plausibly satisfy the first
and fail the second, and would do so only at the very end of the chain.

If the input callsets turn out to be gCNV rather than full GATK-SV, budget for a `NormaliseSvVcf` module
covering both sets of requirements.

### SVAFotate has no releases

The project is not published to PyPI and carries no releases or tags, so `docker/SVAFotate_Dockerfile` pins an
explicit commit SHA (`30b5004a0f4d26959c6b9a82f165651585293626`, dated 2026-07-16). It is actively maintained,
but any version bump is a manual SHA change with no changelog to consult.

## Alternatives considered

| Tool | Verdict |
|------|---------|
| GATK-SV `AnnotateExternalAF` | Emits exactly the `{prefix}_AF` / `{prefix}_SVID` names Talos expects, so needs no rename step. Internally just `svtk vcf2bed` plus `bedtools intersect -r -f`. Rejected in favour of SVAFotate, whose BED is a superset of the GATK frequency TSV, but a reasonable fallback if the SVAFotate image proves hard to maintain |
| AnnotSV | Comprehensive clinical annotation, including gnomAD-SV, DGV, DDD and OMIM. Rejected because it emits TSV rather than VCF, which does not fit the pipeline shape |
| Truvari `collapse` | Good size- and sequence-aware matching, but awkward to repurpose for external AF transfer |
| `vcfanno` | Rejected. Plain interval overlap with no reciprocal-overlap or SVTYPE constraint, so it would annotate a 200 bp DEL with a 3 Mb gnomAD DEL's frequency |
