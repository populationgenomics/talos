"""
import the HPO obo tree
"""

from obonet import read_obo

# import networkx

graph = read_obo('hpo_terms.obo')

# this is just a stub, the panelapp content is not ready
# relevant usage guide:
# https://github.com/dhimmel/obonet/blob/main/examples/go-obonet.ipynb
