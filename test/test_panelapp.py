"""
tests for the PanelApp parser
"""

from copy import deepcopy

import pytest

from talos.models import PanelApp, PanelDetail
from talos.QueryPanelapp import get_best_moi, get_panel, parse_panel_activity

empty_gene_dict = PanelApp(genes={})


def test_activity_parser(panel_activities):
    """
    check that we correctly parse the activities JSON
    """

    activity_dict = parse_panel_activity(panel_activities)

    # this should be absent
    assert 'NOT_GENE' not in activity_dict

    assert 'GENE1' in activity_dict
    assert activity_dict['GENE1'].strftime('%Y-%m-%d') == '2024-02-01'

    assert 'GENE2' in activity_dict
    assert activity_dict['GENE2'].strftime('%Y-%m-%d') == '2024-04-25'

    assert 'GENE3' in activity_dict
    assert activity_dict['GENE3'].strftime('%Y-%m-%d') == '2024-09-15'


@pytest.fixture(name='fake_panelapp')
def fixture_fake_panelapp(requests_mock, latest_mendeliome, latest_incidentalome):
    """
    prepares the web requests mock to serve as stand-in panelapp
    Args:
        requests_mock ():
        latest_mendeliome ():
        latest_incidentalome ():
    """

    requests_mock.register_uri('GET', 'https://panelapp.agha.umccr.org/api/v1/panels/137', json=latest_mendeliome)
    requests_mock.register_uri('GET', 'https://panelapp.agha.umccr.org/api/v1/panels/137/activities', json=[])
    requests_mock.register_uri('GET', 'https://panelapp.agha.umccr.org/api/v1/panels/126', json=latest_incidentalome)
    requests_mock.register_uri('GET', 'https://panelapp.agha.umccr.org/api/v1/panels/126/activities', json=[])


def test_panel_query(fake_panelapp):  # noqa: ARG001
    """
    check that the default parsing delivers correct data
    :param fake_panelapp: fake web hook mock
    """

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, forbidden_genes=set())
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert gd.genes['ENSG00EFGH'].all_moi == {'monoallelic'}


def test_panel_query_removal(fake_panelapp):  # noqa: ARG001
    """
    check that the default parsing delivers correct data
    :param fake_panelapp: fake web hook mock
    """

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, blacklist=['ENSG00EFGH'])
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_forbidden(fake_panelapp):  # noqa: ARG001
    """
    check that the default parsing delivers correct data
    :param fake_panelapp: fake web hook mock
    """
    gd = deepcopy(empty_gene_dict)
    get_panel(gd, forbidden_genes={'ENSG00EFGH'})
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_removal_2(fake_panelapp):  # noqa: ARG001
    """
    check skipping by symbol works as well
    :param fake_panelapp: fake web hook mock
    """

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, blacklist=['EFGH'])
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_forbidden_2(fake_panelapp):  # noqa: ARG001
    """
    check skipping by symbol works as well
    :param fake_panelapp: fake web hook mock
    """

    gd = deepcopy(empty_gene_dict)
    get_panel(gd, forbidden_genes={'EFGH'})
    assert gd.genes['ENSG00ABCD'].all_moi == {'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137}
    assert 'ENSG00EFGH' not in gd.genes


def test_panel_query_addition(fake_panelapp: pytest.fixture):  # noqa: ARG001
    """
    check that the default parsing delivers correct data
    oof, this was a tricky one

    Args:
        fake_panelapp (): fake web hook mock
    """

    # assumed data we already gathered
    gd = PanelApp(
        metadata=[{'version': '0.11088', 'name': 'Mendeliome', 'id': 137}],
        genes={
            'ENSG00ABCD': {
                'symbol': 'ABCD',
                'all_moi': {'monoallelic'},
                'new': [],
                'panels': [137],
            },
            'ENSG00IJKL': {
                'symbol': 'IJKL',
                'all_moi': {'both'},
                'new': [137],
                'panels': [123, 137],
            },
        },
    )

    # should query for and integrate the incidentalome content
    get_panel(gd, panel_id=126)
    assert gd.genes['ENSG00ABCD'].all_moi == {'monoallelic', 'biallelic'}
    assert gd.genes['ENSG00ABCD'].panels == {137, 126}
    assert gd.genes['ENSG00IJKL'].all_moi == {'both'}
    assert gd.genes['ENSG00IJKL'].panels == {123, 137}
    assert 'ENSG00EFGH' not in gd.genes


def test_get_best_moi_empty():
    """
    check that the MOI summary works
    """

    d = {'ensg1': PanelDetail(all_moi=set(), chrom='1', symbol='ensg1')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Biallelic'

    d = {'ensg1': PanelDetail(all_moi=set(), chrom='X', symbol='ensgX')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Hemi_Bi_In_Female'


def test_get_best_moi_mono():
    """
    check that the MOI summary works
    """

    d = {'ensg1': PanelDetail(all_moi={'monoallelic'}, chrom='1', symbol='ensg1')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Monoallelic'


def test_get_best_moi_mono_and_biallelic():
    """
    check that the MOI summary works
    """

    d = {'ensg1': PanelDetail(all_moi={'monoallelic', 'biallelic'}, chrom='1', symbol='ensg1')}
    get_best_moi(d)
    assert d['ensg1'].moi == 'Mono_And_Biallelic'


def test_get_best_moi_1():
    """
    check that the MOI summary works
    """

    d = {
        'ensg1': PanelDetail(
            all_moi={'Monoallelic', 'Biallelic', 'both'},
            chrom='1',
            symbol='ensg1',
        ),
    }
    get_best_moi(d)
    assert d['ensg1'].moi == 'Mono_And_Biallelic'


def test_get_best_moi_x():
    """
    check that the MOI summary works
    """

    d = {
        'ensg1': PanelDetail(
            all_moi={'x-linked biallelic', 'x-linked'},
            chrom='X',
            symbol='ensgX',
        ),
    }
    get_best_moi(d)
    assert d['ensg1'].moi == 'Hemi_Mono_In_Female'
