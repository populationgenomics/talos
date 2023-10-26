# Annotation

AIP leverages numerous functions within the Hail Batch and Hail Query libraries to organise workflows and efficiently
query large datasets. Unfortunately at time of writing, Hail's ability to run annotation using VEP is limited to GCS
DataProc instances. To allow for a more flexible implementation, the CPG has created an annotation workaround,
within the [cpg_workflows package](https://github.com/populationgenomics/production-pipelines/tree/main/cpg_workflows).
This acts directly on a VCF, fragmenting the raw data and annotating in parallelised jobs, forming the annotated data
back into a [Hail MatrixTable](https://hail.is/docs/0.2/hail.MatrixTable.html), which is the AIP starting point.
