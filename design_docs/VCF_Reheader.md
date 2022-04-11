# VCF Re-Header

Following the labelling stage (currently in Hail Query), A VCF is generated which contains a number of new fields. For
most of those fields, the header content contains a descriptive enough statement, e.g.

```bash
##INFO=<ID=exac_ac_het,Number=1,Type=Integer,Description="">
##INFO=<ID=gnomad_af,Number=1,Type=Float,Description="">
##INFO=<ID=Category1,Number=1,Type=Integer,Description="">
```

These fields are annotated with the name and type, which is typically descriptive enough. One annotation which requires
separate explanation is the `CSQ` field. This is created as an aggregation of separate values (present or missing) from
a number of component fields. For additional complexity the CSQ annotation can contain multiple sections,
comma-delimited, which represent a separate group of annotations all relating to a separate transcript.

A CSQ string description may look like this:

```bash
Allele|Consequence|SYMBOL|Gene|Feature|MANE_select|BIOTYPE|EXON|HGVSc|HGVSp|cDNA_position|CDS_position|
Protein_position|Amino_acids|Codons|ALLELE_NUM|VARIANT_CLASS|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|
SIFT|PolyPhen|LoF|LoF_filter|LoF_flags
```

This description is obtained from the configuration file, and used within the Hail-Query code in the following way:

1. For each variant, create a key:value dictionary of each struct within the VEP.transcript_consequences section
2. Parse the configuration string, identifying the values & ordering to keep
3. For each transcript consequence, create a "|"-delimited string of all the values to keep
4. Add the list of new Strings into the MatrixTable as an annotation

In order to use the VCF when separated from the configuration file, we have to update the VCF header to include the
values & value-order encoded into the CSQ string. To do this we use the BCFtools reheader command:

1. Use BCFtools to read all header lines
2. Pipe through SED to edit the description line
3. Use BCFtools reheader to apply the new header to the original file

Using the example above, this changes

```bash
##INFO=<ID=CSQ,Number=.,Type=String,Description="">
```

to

```bash
##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|SYMBOL|Gene|Feature|MANE_select|BIOTYPE|
EXON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|VARIANT_CLASS|TSL|APPRIS|
CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|LoF|LoF_filter|LoF_flags">
```
