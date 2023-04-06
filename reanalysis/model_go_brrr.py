"""
example script for using the data model

Almost all elements of the data model are optional,
Other than locus/alleles, and gene ID/symbol.

The data model is designed to be flexible, and can be used to create
a variety of different data structures, from a single variant to a
multi-sample VCF.

This example script shows how to create a single variant, with a single
transcript consequence, and a trio of samples.
"""

import hail as hl

from reanalysis.data_model import BaseFields, TXFields, VepVariant, SneakyTable, Entry


# for this minimal example, we'll just use a single transcript consequence with all
# default values, i.e. no consequence, just a gene symbol & ID
t = TXFields('a', 'ensga')

# let's make a weirdly phased trio
sample_data = {'sam1': Entry('0|1'), 'sam2': Entry('0|1'), 'sam3': Entry('1|0')}

# and for each of these Sample entries take the default schema
sample_schema = {k: v.get_schema_entry() for k, v in sample_data.items()}

# create a single VepVariant object, with the sample data and TX CSQ
v = VepVariant(BaseFields('chr1:123456', ['A', 'C']), [t], sample_data=sample_data)

# create a SneakyTable object, which will take a list of VepVariant objects
sn = SneakyTable([v], sample_schema, '.')

# once the table is parsed using the JSON schema, convert to a Hail MatrixTable
# using the row_major functionality - this takes the list of sample IDs, converts those
# to columns, and converts the calls from Strings to proper GT values
mt = sn.to_hail()

# blast off! (this was a copilot suggestion, couldn't resist)
mt.describe()
mt.entries().show()

# this is a valid VCF, so we can export it (drops all 'VEP' annotations)
hl.export_vcf(mt, 'mt.vcf.bgz', tabix=True)
