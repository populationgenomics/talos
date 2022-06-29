# setting up a new run

CPG-specific guidance on how to generate input files for a new run

## Initial checks

1. Check that the [token map](https://github.com/populationgenomics/analysis-runner/blob/main/tokens/repository-map.json) allows AIP to run against the relevant project's bucket & contents
2. Check that you have permission to query DeepPan for the relevant project (e.g. you can run queries for the project on the [Swagger page](https://sample-metadata.populationgenomics.org.au/swagger))
3. Check that you have permission to upload relevant files to the project bucket (needs working out, current process involves uploading files to the permissive `-test` bucket due to ease)

## Generating files

(You'll want a cohort and run-specific folder to copy these files to. This can also be the output folder for the run)

1. Query DeepPan for the pedigree and internal -> external lookup files
   - [cohort_prep.sh](helpers/cohort_prep.sh) is a shortcut for this process
   - `bash helpers/cohort_prep.sh COHORT [% of cohort to retain]` will create files in `inputs/COHORT`
   - The file ending `.fam` is the pedigree file. Copy to GCP, and the path becomes an input argument
   - The file ending `external_lookup.json` is used to translate CPG to external IDs. Copy to GCP and enter the file path into the run config `output.external_lookup`

2. Get the seqr Family lookup
   - Optional, if provided the final HTML will hyperlink to Seqr variants and families
   - I realise that this would be far simpler to understand with pictures, but almost everything in screenshots is variant/family data...
   - Clunky process, aiming to improve soon
   - Go to a relevant project page in Seqr (clicking the project name)
   - Check that `output.seqr_instance` in the config file matches the seqr homepage
   - Identify the project name (e.g. in `https://testseqr.org.au/project/R00001_test_project/project_page`, get `R00001_test_project`) - add this to the config as `output.seqr_project`
   - (in Chrome) Burger -> More Tools -> Developer Tools (Opt+Cmd+I) & Open the Network tab
   - Refresh, then find the named request `get_overview`
   - On the right-side of the console, find the `Response` tab, which should show some JSON
   - Select all this JSON and save to a file
   - Run the [process_seqr_metadata script](helpers/process_seqr_metadata.py), with this file as input (`-i`) and your chosen path as output (`-o`)
   - Copy this file to GCP, and add the path to the config as `output.seqr_lookup`
3. Copy any additional files (vqsr header, csq header line) to the relevant bucket, and update paths in config
4. Check any other configuration values (e.g. thresholds)
5. Copy the config file to GCP

## Starting the run

Create a run script based on the example [run_reanalysis.sh](reanalysis/run_reanalysis.sh)

All paths referenced here must be GCP bucket paths, as no files will be copied to GCP at runtime. For input genomic data
this will involve finding the appropriate paths to use. For input files this will involve transcribing the paths you
previously created into this script as inputs.
