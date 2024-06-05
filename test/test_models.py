"""
any methods for testing model functionality
"""

from reanalysis.utils import make_flexible_pedigree
from test.test_utils import FIVE_EXPECTED, THREE_EXPECTED, TWO_EXPECTED


def test_flexi_pedigree(test_input_path):
    input_ped = str(test_input_path / 'peds.ped')
    flexi_ped = make_flexible_pedigree(input_ped)
    assert len(flexi_ped.by_family) == THREE_EXPECTED
    assert len(flexi_ped.by_id) == FIVE_EXPECTED
    assert len(flexi_ped.members) == FIVE_EXPECTED
    trio_proband = flexi_ped.by_id['CPGABC1']
    assert trio_proband.father == 'CPGABC2'
    assert trio_proband.mother == 'CPGABC3'
    assert trio_proband.ext_id == 'EXT1'
    assert len(trio_proband.hpo_terms) == TWO_EXPECTED
    assert sorted([hpo.id for hpo in trio_proband.hpo_terms]) == ['HPO:1', 'HPO:2']

    # this has no HPO terms
    assert len(flexi_ped.by_id['CPGABC5'].hpo_terms) == 0
