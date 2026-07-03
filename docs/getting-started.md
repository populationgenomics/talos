# Getting Started

This guide walks you through preparing the environment, downloading the required reference data, and running your first Talos analysis.

Talos is implemented using **Nextflow**, with all dependencies containerised via Docker. The example workflows can be run locally or on a cluster.

There are two primary entry points:

- `preparation.nf` — downloads and formats data in preparation for Talos runs.
- `main.nf` — imports and executes the two sub-workflows (`annotation` + `Talos`) that comprise a full run.

And one secondary entrypoint:

- `talos_only.nf` - runs the analysis portion of the workflow only, requires a previous completion of the `annotation` wf.

---

## 1. Install requirements

You will need:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- Docker (or a compatible container runtime)

Build the Talos Docker image locally:

```bash
docker build -t talos:11.0.3 .
```

---

## 2. Download annotation resources

Talos requires several large external resources (reference genome, gnomAD, AlphaMissense, phenotype data). These are expected in a `large_files/` directory.

See [`large_files/README.md`](https://github.com/populationgenomics/talos/blob/main/large_files/README.md) in the repository for the full list of sources, and the helper script [`large_files/gather_files.sh`](https://github.com/populationgenomics/talos/blob/main/large_files/gather_files.sh), which performs the initial download.

---

## 3. Run the preparation workflow

In addition to the raw downloads, Talos requires the following annotation sources to be kept up to date:

- **ClinVar** data, formatted into Hail Tables.
- **PanelApp** data, dumped to JSON.
- **AlphaMissense** data, reformatted from the source TSV into a Hail Table.

The `preparation.nf` sub-workflow handles all three:

```bash
nextflow \
    -c nextflow.config \
    run preparation.nf \
    [--processed_annotations <path>] \
    [--large_files <path>]
```

The `processed_annotations` parameter should point to a static directory where Talos-generated files will be stored. Files written here are reusable across all future Talos runs; they are not tied to any individual cohort or callset.

!!! info "Recommended cadence"
    Re-run `preparation.nf` on a regular schedule (e.g. monthly) so that ClinVar and PanelApp evidence stays current between Talos analyses. Talos is set up to complain and exit if ClinVar and PanelApp data wasn't prepared this month, as determined by the date built into the file names.

---

## 4. Prepare your input TSV

From version `10.0.0` onwards, all per-cohort inputs are provided in a single TSV file via `--input_tsv`. Each row in the TSV represents one cohort, and the workflow runs them in parallel into separate output directories.

The optional columns (history, ext_ids, seqr_map, mito) can be omitted completely. If they are not provided, NextFlow defaults to a real but empty dummy file. See the provided example input file [here](https://github.com/populationgenomics/talos/blob/main/nextflow/inputs/test.tsv)

| Column     | Required | Description                                                                    |
|:-----------|:---------|:-------------------------------------------------------------------------------|
| `cohort`   | yes      | Collective name used in output paths and file naming.                          |
| `path`     | yes      | Path to the cohort's variant input.                                            |
| `type`     | yes      | One of `vcf`, `shards`, or `ss_vcf_dir` (see below).                           |
| `pedigree` | yes      | Path to the pedigree file. See [Pedigree & Phenotypes](Pedigree.md).           |
| `config`   | yes      | Path to the cohort's Talos TOML config. See [Configuration](Configuration.md). |
| `history`  | optional | Path to previous results, enabling reanalysis mode.                            |
| `ext_ids`  | optional | ID mapping for alternate IDs in the HTML report.                               |
| `seqr_map` | optional | ID mapping for Seqr hyperlinks in the HTML report.                             |
| `mito`     | optional | Path to a joint-called mitochondrial VCF.                                      |

### Input types

<table>
    <colgroup>
        <col style="width: 8rem">
        <col>
    </colgroup>
    <thead><tr><th>Type</th><th>Description</th></tr></thead>
    <tbody>
        <tr><td><code>vcf</code></td><td>A single multi-sample VCF; split into shards and processed in parallel.</td></tr>
        <tr><td><code>shards</code></td><td>A directory of pre-sharded multi-sample VCF fragments, each shard containing all samples.</td></tr>
        <tr><td><code>ss_vcf_dir</code></td><td>A directory of single-sample VCFs, merged then sharded inside the workflow. The glob extension is controlled by `params.input_vcf_extension` (defaults to `vcf.bgz`).</td></tr>
    </tbody>
</table>

---

## 5. Run the annotation + Talos workflow

```bash
nextflow \
    -c nextflow.config \
    run main.nf \
    --input_tsv nextflow/inputs/test.tsv \
    -output-dir <path_to_output_dir>
```

Results are written to `{workflow.outputDir}/{cohort}_outputs`. The annotation sub-workflow only needs to be run once per dataset — the resulting MatrixTables can be reused for every subsequent reanalysis cycle.

For subsequent cycles after the data has already been annotated, use the Talos-only entry point:

```bash
nextflow \
    -c nextflow.config \
    run talos_only.nf \
    --input_tsv nextflow/inputs/test.tsv \
    -output-dir <path_to_output_dir>
```

!!! tip "Reanalysis cadence"
    For best results, re-run the Talos workflow on a regular cadence — monthly or quarterly is typical. See the [Reanalysis](Reanalysis.md) page for how prior results are folded into new runs.

---

## 6. Check the outputs

A successful run produces:

- A `*.json` file listing all candidate variants for each proband, with variant- and gene-level evidence, inheritance checks, and phenotype tags.
- An optional **HTML report** for analyst/clinician review.
- An optional **simplified TSV** for Seqr ingestion via `MinimiseOutputForSeqr`.

Each variant carries reanalysis metadata (`first_seen`, `evidence_last_updated`) so downstream consumers can filter to newly actionable findings.

---

## What happens during startup

The first step of every Talos run is a `StartupChecks` module that validates inputs before any analysis begins:

1. Confirms a config file is present and that all required entries exist with the correct types.
2. Opens the MatrixTable and checks the schema and data types.
3. Parses the pedigree file, validating its format and that affected participants are present.
4. Checks the ClinVar data, ensuring it is recent and has sufficient entries.

If startup fails, Talos prints a collected list of all encountered errors. Fix these before restarting the workflow.

---

## Next steps

- Tune the rule-based logic via [Talos Configuration](Configuration.md).
- Tune Nextflow runtime and resource settings via [Nextflow Configuration](NextflowConfiguration.md).
- Read [Features](features.md) for an overview of the logic modules and reanalysis behaviour.
