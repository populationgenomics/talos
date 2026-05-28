from talos.download_panelapp import (
    PANELS_ENDPOINT,
    get_panels_and_hpo_terms,
    parse_panel,
    parse_panel_activity,
)
from talos.liftover.lift_2_2_0_to_2_3_0 import dl_panelapp as dl_pa_220_to_230
from talos.models import (
    CURRENT_VERSION,
    DownloadedPanelApp,
    DownloadedPanelAppGenePanelDetail,
    HpoTerm,
    lift_up_model_version,
)


def test_panel_hpo_query(httpx_mock, panels_and_hpos):
    """check that the default parsing delivers correct data"""

    httpx_mock.add_response(url=PANELS_ENDPOINT, json=panels_and_hpos)

    parsed_response = get_panels_and_hpo_terms()

    assert parsed_response == {
        3149: [HpoTerm(id='HP:0011516', label='')],
        4059: [
            HpoTerm(id='HP:0001638', label=''),
            HpoTerm(id='HP:0001637', label=''),
            HpoTerm(id='HP:0011675', label=''),
        ],
        3302: [],
    }


def test_activity_parser(panel_activities):
    """
    check that we correctly parse the activities JSON
    """

    activity_dict = parse_panel_activity(panel_activities)

    # this should be absent
    assert 'NOT_GENE' not in activity_dict

    assert 'GENE1' in activity_dict
    assert activity_dict['GENE1'] == '2022-02-01'

    assert 'GENE2' in activity_dict
    assert activity_dict['GENE2'] == '2024-04-25'

    assert 'GENE3' in activity_dict
    assert activity_dict['GENE3'] == '2023-08-15'


def test_parse_panel(latest_mendeliome, panel_activities):
    result = parse_panel(panel_data=latest_mendeliome, panel_activities=panel_activities)
    assert result == {
        'ENSG00ABCD': {
            'symbol': 'ABCD',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'biallelic',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
        'ENSG00EFGH': {
            'symbol': 'EFGH',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'monoallelic',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
        'ENSG00IJKL': {
            'symbol': 'IJKL',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'both',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
    }


def test_parse_panel_with_mane(latest_mendeliome, panel_activities):
    """check that the default parsing delivers correct data"""
    fake_mane_symbols = {'ABCD': 'EasyAs123D'}
    result = parse_panel(panel_data=latest_mendeliome, panel_activities=panel_activities, symbol_dict=fake_mane_symbols)
    assert result == {
        'EasyAs123D': {
            'symbol': 'ABCD',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'biallelic',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
        'ENSG00ABCD': {
            'symbol': 'ABCD',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'biallelic',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
        'ENSG00EFGH': {
            'symbol': 'EFGH',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'monoallelic',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
        'ENSG00IJKL': {
            'symbol': 'IJKL',
            'chrom': '1',
            'mane_symbol': '',
            'moi': 'both',
            'green_date': '1970-01-01',
            'confidence_level': 3,
        },
    }


_GENE_TEMPLATE = {
    'gene_data': {
        'ensembl_genes': {
            'GRch38': {
                '90': {
                    'ensembl_id': 'ENSG{symbol}',
                    'location': '1:',
                },
            },
        },
    },
    'entity_type': 'gene',
    'mode_of_inheritance': 'Biallelic',
}


def _make_gene(symbol: str, confidence: int) -> dict:
    return {
        **_GENE_TEMPLATE,
        'entity_name': symbol,
        'confidence_level': str(confidence),
        'gene_data': {
            'ensembl_genes': {
                'GRch38': {
                    '90': {
                        'ensembl_id': f'ENSG{symbol}',
                        'location': '1:',
                    },
                },
            },
        },
    }


def test_parse_panel_includes_all_when_threshold_one(panel_activities):
    """setting GENE_CONFIDENCE to 1 admits red, amber, and green genes"""
    panel_data = [
        _make_gene('GREEN', 3),
        _make_gene('AMBER', 2),
        _make_gene('RED', 1),
    ]
    result = parse_panel(panel_data=panel_data, panel_activities=panel_activities)
    assert 'ENSGGREEN' in result
    assert 'ENSGAMBER' in result
    assert 'ENSGRED' in result
    assert result['ENSGRED']['confidence_level'] == 1


def test_downloaded_panel_app_gene_panel_detail_confidence():
    """DownloadedPanelAppGenePanelDetail must persist the confidence field"""
    detail = DownloadedPanelAppGenePanelDetail(moi='biallelic', date='2024-01-01', confidence=2)
    assert detail.confidence == 2

    detail_default = DownloadedPanelAppGenePanelDetail(moi='monoallelic')
    assert detail_default.confidence == 0


def test_liftover_dl_panelapp_220_to_230():
    """liftover function must add confidence=3 to every panel entry and bump version"""
    data = {
        'version': '2.2.0',
        'genes': {
            'ENSG001': {
                'panels': {
                    137: {'moi': 'biallelic', 'date': '2023-01-01'},
                    42: {'moi': 'monoallelic', 'date': '2022-06-15'},
                },
            },
            'ENSG002': {
                'panels': {
                    137: {'moi': 'both', 'date': '2021-03-10'},
                },
            },
        },
    }
    result = dl_pa_220_to_230(data)
    assert result['version'] == '2.3.0'
    for gene_data in result['genes'].values():
        for panel_data in gene_data['panels'].values():
            assert panel_data['confidence'] == 3


def test_liftover_downloaded_panelapp_via_model():
    """lift_up_model_version must walk the chain from 2.2.0 to current and produce a valid model"""
    data = {
        'version': '2.2.0',
        'genes': {
            'ENSG001': {
                'symbol': 'GENE1',
                'chrom': '1',
                'mane_symbol': '',
                'ensg': 'ENSG001',
                'panels': {
                    137: {'moi': 'biallelic', 'date': '2023-01-01'},
                },
            },
        },
        'versions': [],
        'hpos': {},
    }
    lifted = lift_up_model_version(data, model=DownloadedPanelApp)
    assert lifted['version'] == CURRENT_VERSION
    parsed = DownloadedPanelApp.model_validate(lifted)
    assert parsed.genes['ENSG001'].panels[137].confidence == 3
