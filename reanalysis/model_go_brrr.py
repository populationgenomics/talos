"""
example script for using the data model
"""

from reanalysis.data_model import BaseFields, TXFields, VepVariant, SneakyTable


b = BaseFields('chr1:123456', ['A', 'C'])
t = TXFields('a', 'ensga')
v = VepVariant([t], b)
sn = SneakyTable([v], '.')
ht = sn.to_hail()

ht.show()
