"""
script used to create an artificial VCF for testing
"""

import hail as hl

from talos.data_model import BaseFields, Entry, SneakyTable, TXFields, VepVariant

hl.init(default_reference='GRCh38')

comp_het = {'CPG255208': Entry('1/1'), 'CPG255224': Entry('0/1'), 'CPG255216': Entry('0/1')}
de_novo = {'CPG255208': Entry('0/1'), 'CPG255224': Entry('0/0'), 'CPG255216': Entry('0/0')}
mat_inherited = {'CPG255208': Entry('0/1'), 'CPG255224': Entry('0/1'), 'CPG255216': Entry('0/0')}
pat_inherited = {'CPG255208': Entry('0/1'), 'CPG255224': Entry('0/0'), 'CPG255216': Entry('0/1')}

# dull example, single tx consequence, all default values, i.e. only gene symbol & ID
t = TXFields('a', 'ensga')

# and for each of these Sample entries take the default schema
sample_schema = {k: v.get_schema_entry() for k, v in comp_het.items()}

# Cat 1, 3, dominant, WT1
v1 = VepVariant(BaseFields('chr11:32392032', ['G', 'A']), [t], sample_data=de_novo)

# cat 6, dominant, DAAM2
v2 = VepVariant(BaseFields('chr6:39887558', ['C', 'T']), [t], sample_data=de_novo)

# cats 1 5 3, POC1B, AR
v3 = VepVariant(BaseFields('chr12:89470359', ['A', 'C']), [t], sample_data=comp_het)

# cat pm5, dominant, PKHD1
v4 = VepVariant(BaseFields('chr6:51659900', ['T', 'C']), [t], sample_data=de_novo)

# cat 1, recessive, HFE
v5 = VepVariant(BaseFields('chr6:26090951', ['C', 'G']), [t], sample_data=comp_het)

loci = [
    hl.Locus('chr11', 32392032),
    hl.Locus('chr6', 39887558),
    hl.Locus('chr12', 89470359),
    hl.Locus('chr6', 51659900),
    hl.Locus('chr6', 26090951),
]

# create a SneakyTable object, which will take a list of VepVariant objects
sn = SneakyTable([v1, v2, v3, v4, v5], '.', sample_schema)

# once the table is parsed using the JSON schema, convert to a Hail MatrixTable
# using the row_major functionality - this takes the list of sample IDs, converts those
# to columns, and converts the calls from Strings to proper GT values
mt = sn.to_hail()
mt.describe()

# this is a valid VCF, so we can export it (drops all 'VEP' annotations)
hl.export_vcf(mt, 'test_1.vcf.bgz', tabix=True)
# endregion