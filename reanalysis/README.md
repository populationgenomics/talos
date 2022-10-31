# Running on Azure Hail Batch

1) Get the Storage Account Connection String

Go to your Storage Account on Azure then:
<Storage Account> > Access keys > Connection string

Copy then set the following environment variable:
```bash
set -x AZURE_STORAGE_CONNECTION_STRING 'DefaultEndpointsProtocol=https;AccountName=blablah;AccountKey=longkyestringthing;EndpointSuffix=core.windows.net'
```

2) Run the following command

> Note: The storage account does not need to be in the string since
> the AZURE_STORAGE_CONNECTION_STRING encodes that information

```bash
./aip/reanalysis/interpretation_runner.py \
    -i az://test/annotated-file.mt \
    --pedigree az://test/pedigree.fam \
    --skip_annotation
````
