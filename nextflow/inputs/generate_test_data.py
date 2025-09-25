"""
script used to create an artificial VCF for testing
"""

import hail as hl

from talos.data_model import BaseFields, Entry, SneakyTable, TXFields, VepVariant

hl.default_reference('GRCh38')

comp_het = {'proband': Entry('1/1'), 'mother': Entry('0/1'), 'father': Entry('0/1')}
de_novo = {'proband': Entry('0/1'), 'mother': Entry('0/0'), 'father': Entry('0/0')}
mat_inherited = {'proband': Entry('0/1'), 'mother': Entry('0/1'), 'father': Entry('0/0')}
pat_inherited = {'proband': Entry('0/1'), 'mother': Entry('0/0'), 'father': Entry('0/1')}

# this is for a parsing check that we can accurately process hemizygous alleles
pat_inherited_hemi = {'proband': Entry('1'), 'mother': Entry('0/0'), 'father': Entry('1')}
mat_inherited_hemi = {'proband': Entry('1'), 'mother': Entry('0/1'), 'father': Entry('0')}

# dull example, single tx consequence, all default values, i.e. only gene symbol & ID
t = TXFields('a', 'ensga')

# and for each of these Sample entries take the default schema
sample_schema = {k: v.get_schema_entry() for k, v in comp_het.items()}

# cat 3, no canonical consequences, here to check normalisation
vnorm = VepVariant(BaseFields('chr1:21706892', ['GA', 'GAA']), [t], sample_data=comp_het)

# cat 1, recessive, HFE - not retained by the filtering logic when HFE is a phenotype-match only gene
v1 = VepVariant(BaseFields('chr6:26090951', ['C', 'G']), [t], sample_data=comp_het)

# PKHD1 (AR) | v4a pm5 | v4b cat1 | both should be reported
v2a = VepVariant(BaseFields('chr6:52043699', ['T', 'A']), [t], sample_data=de_novo)
v2b = VepVariant(BaseFields('chr6:52043102', ['C', 'G']), [t], sample_data=mat_inherited)

# DARS (AR) | v3a cat3 | v3b cat6 | v3a should be reported, v3b supporting only
v3a = VepVariant(BaseFields('chr2:135920591', ['G', 'C']), [t], sample_data=pat_inherited)
v3b = VepVariant(BaseFields('chr2:135912503', ['G', 'A']), [t], sample_data=mat_inherited)

# cat 6, dominant, DAAM2 - will be filtered out
v4 = VepVariant(BaseFields('chr6:39887558', ['C', 'T']), [t], sample_data=de_novo)

# Cat 1, 3, dominant, WT1
v5 = VepVariant(BaseFields('chr11:32392032', ['G', 'A']), [t], sample_data=de_novo)

# cats 1 5 3, POC1B, AR
v6 = VepVariant(BaseFields('chr12:89470359', ['A', 'C']), [t], sample_data=comp_het)

# IL2RG (hemi/bi in females) | cat 1 | should be reported
v7 = VepVariant(BaseFields('chrX:71109321', ['G', 'A']), [t], sample_data=mat_inherited_hemi)

# rnu4-2 (AD/AR) | cat 1, de novo | should be reported
v8 = VepVariant(BaseFields('chr12:120291834', ['A', 'G']), [t], sample_data=de_novo)

# create a SneakyTable object, which will take a list of VepVariant objects
sn = SneakyTable([vnorm, v1, v2a, v2b, v3a, v3b, v4, v5, v6, v7, v8], '', sample_schema)

# once the table is parsed using the JSON schema, convert to a Hail MatrixTable
# using the row_major functionality - this takes the list of sample IDs, converts those
# to columns, and converts the calls from Strings to proper GT values
mt = sn.to_hail()
mt.describe()

# this is a valid VCF, so we can export it (drops all 'VEP' annotations)
hl.export_vcf(mt, 'test_1.vcf.bgz', tabix=True)

for sample_id in ['mother', 'father', 'proband']:
    sample_mt = mt.filter_cols(mt.s == sample_id)
    # this will create a VCF with the sample IDs as the column names
    hl.export_vcf(sample_mt, f'{sample_id}.vcf.bgz', tabix=True)
