"""
script used to create an artificial VCF for testing
"""

import hail as hl

from talos.data_model import BaseFields, Entry, SneakyTable, TXFields, VepVariant

hl.default_reference('GRCh38')

comp_het = {'proband': Entry('1/1'), 'mother': Entry('0/1'), 'father': Entry('0/1')}
de_novo = {'proband': Entry('0/1'), 'mother': Entry('0/0'), 'father': Entry('0/0')}

mat_inherited = {'proband': Entry('0/1'), 'mother': Entry('0/1'), 'father': Entry('0/0')}

# dull example, single tx consequence, all default values, i.e. only gene symbol & ID
t = TXFields('a', 'ensga')

# and for each of these Sample entries take the default schema
sample_schema = {k: v.get_schema_entry() for k, v in comp_het.items()}

# cat 3, no canonical consequences, here to check normalisation
v1 = VepVariant(BaseFields('chrM:3243', ['A', 'G']), [t], sample_data=mat_inherited)

v2 = VepVariant(BaseFields('chrM:8528', ['T', 'C']), [t], sample_data=mat_inherited)

v3 = VepVariant(BaseFields('chrM:12278', ['T', 'C']), [t], sample_data=mat_inherited)

# create a SneakyTable object, which will take a list of VepVariant objects
sn = SneakyTable([v1, v2, v3], '', sample_schema)

# once the table is parsed using the JSON schema, convert to a Hail MatrixTable
# using the row_major functionality - this takes the list of sample IDs, converts those
# to columns, and converts the calls from Strings to proper GT values
mt = sn.to_hail()
mt.describe()

# this is a valid VCF, so we can export it (drops all 'VEP' annotations)
hl.export_vcf(mt, 'joint_mito.vcf.bgz', tabix=True)
