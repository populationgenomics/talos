"""
any methods for testing model functionality
"""

from reanalysis.utils import make_flexible_pedigree


def test_flexi_pedigree(test_input_path):
    input_ped = str(test_input_path / 'peds.ped')
    flexi_ped = make_flexible_pedigree(input_ped)
    print(flexi_ped)
    assert len(flexi_ped.by_family) == 3
    assert len(flexi_ped.by_id) == 5
    assert len(flexi_ped.members) == 5
    trio_proband = flexi_ped.by_id['CPGABC1']
    assert trio_proband.father == 'CPGABC2'
    assert trio_proband.mother == 'CPGABC3'
    assert trio_proband.ext_id == 'EXT1'
    assert len(trio_proband.hpo_terms) == 2
    assert sorted([hpo.id for hpo in trio_proband.hpo_terms]) == ['HPO:1', 'HPO:2']

    # this has no HPO terms
    assert len(flexi_ped.by_id['CPGABC5'].hpo_terms) == 0
