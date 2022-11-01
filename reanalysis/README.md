# Running on Azure Hail Batch

0) Get the right repos/branches cloned and installed all into a single conda env

populationgenomics/automated-interpretation-pipeline:syan
gregsmi/cpg-utils:main
jeremiahwander/analysis-runner-azcpg004:main

1) Set the hail backend to Azure

```bash
cat > ~/.hail/deploy-config.json <<EOF
{
  "location": "external",
  "default_namespace": "default",
  "domain": "azhail.popgen.rocks"
}
EOF
```

2) Get the Storage Account Connection String

Go to your Storage Account on Azure then:
<Storage Account> > Access keys > Connection string

Copy then set the following environment variable:
```bash
set -x AZURE_STORAGE_CONNECTION_STRING 'DefaultEndpointsProtocol=https;AccountName=blablah;AccountKey=longkyestringthing;EndpointSuffix=core.windows.net'
```

3) Create the cpg_config toml file for the run and set the env
   variable appropriately

```bash
cat > master.toml <<EOF
[buckets]
web_suffix = "web"
tmp_suffix = "tmp"
analysis_suffix = "analysis"

[workflow]
dataset_path = "cpgfewgenomes00fd84e4/test"
output_prefix = "output"
path_scheme = "hail-az"
image_registry_prefix = "cpghailimages.azurecr.io"

[hail]
billing_project = "fewgenomes"
bucket = "hail-az://cpgfewgenomes00fd84e4/test"

[images]
hail = "hailgenetics/hail:0.2.93"
EOF
set -x CPG_CONFIG_PATH (realpath master.toml)
```

4) Run the following command

> Note: The storage account does not need to be in the string since
> the AZURE_STORAGE_CONNECTION_STRING encodes that information

```bash
./aip/reanalysis/interpretation_runner.py \
    -i az://test/986d792a448c66a8a5cfba65434e7d1ce9b1ff_1051-validation.mt \
    --pedigree az://test/pedigree.fam \
    --skip_annotation
````
