"""
any methods for testing model functionality
"""

from os.path import join
from test.test_utils import FIVE_EXPECTED, ONE_EXPECTED, THREE_EXPECTED

from talos.models import PanelApp
from talos.utils import make_flexible_pedigree


def test_flexi_pedigree(test_input_path):
    # build the phenotype matched panels object
    sample_panels = PanelApp(
        participants={
            'CPGABC1': {'panels': {1, 3}, 'external_id': 'EXT1', 'hpo_terms': [{'id': 'HPO:1', 'label': 'Boneitis!'}]},
            'female': {
                'panels': {1, 2},
                'external_id': 'female',
                'hpo_terms': [{'id': 'HPF', 'label': 'HPFemale'}],
            },
        },
    )

    input_ped = join(test_input_path, 'peds.ped')

    # requires ext IDs in a PhenoMatchedPanels object
    flexi_ped = make_flexible_pedigree(input_ped, sample_panels)
    assert len(flexi_ped.by_family) == THREE_EXPECTED
    assert len(flexi_ped.by_id) == FIVE_EXPECTED
    assert len(flexi_ped.members) == FIVE_EXPECTED
    trio_proband = flexi_ped.by_id['CPGABC1']
    assert trio_proband.father == 'CPGABC2'
    assert trio_proband.mother == 'CPGABC3'
    assert trio_proband.ext_id == 'EXT1'
    assert len(trio_proband.hpo_terms) == ONE_EXPECTED
    assert [hpo.id for hpo in trio_proband.hpo_terms] == ['HPO:1']

    # this has no HPO terms
    assert len(flexi_ped.by_id['CPGABC5'].hpo_terms) == 0


def test_flexi_pedigree_no_packets(test_input_path):
    input_ped = join(test_input_path, 'peds.ped')

    # requires ext IDs in a PhenoMatchedPanels object
    flexi_ped = make_flexible_pedigree(input_ped)
    assert len(flexi_ped.by_family) == THREE_EXPECTED
    assert len(flexi_ped.by_id) == FIVE_EXPECTED
    assert len(flexi_ped.members) == FIVE_EXPECTED
    trio_proband = flexi_ped.by_id['CPGABC1']
    assert trio_proband.father == 'CPGABC2'
    assert trio_proband.mother == 'CPGABC3'
    assert trio_proband.ext_id == 'Missing'
    assert len(trio_proband.hpo_terms) == 0

    # this has no HPO terms
    assert len(flexi_ped.by_id['CPGABC5'].hpo_terms) == 0
    assert flexi_ped.by_id['CPGABC5'].family == 'FAM3'
