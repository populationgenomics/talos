"""
unit testing collection for the hail MT methods
"""

import pytest

import hail as hl

from talos.models import PanelApp, PanelDetail
from talos.run_hail_filtering import (
    filter_matrix_by_ac,
    filter_on_quality_flags,
    green_gene_intervals,
    parse_gene_location,
    populate_callset_frequencies,
)


def test_no_freq_replacement(make_a_mt, caplog):
    """check that when AC/AF/AN are all present, no replacement is attempted"""
    matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            AC=[1],
            AN=10,
            AF=[0.1],
        ),
    )
    _matrix = populate_callset_frequencies(matrix)
    assert 'AC, AN, AF already present, skipping annotation' in caplog.text, caplog.text


def test_af_replacement(make_a_mt, caplog):
    """check that when AC and AN are present, AF is replaced"""
    matrix = make_a_mt.annotate_rows(
        info=hl.struct(
            AC=[1],
            AN=10,
        ),
    )
    matrix_out = populate_callset_frequencies(matrix)
    assert 'AC, AN present, deriving AF from existing annotations' in caplog.text, caplog.text
    assert matrix_out.info.AF.collect()[0] == [0.1]


def test_full_freq_replacement(make_a_mt, caplog):
    """check that when AC and AN annotations are present, AF is replaced"""
    matrix = make_a_mt.annotate_rows(
        info=hl.struct(),
    )
    matrix_out = populate_callset_frequencies(matrix)
    assert 'Adding AC/AN/AF annotations to MT based on this callset alone' in caplog.text, caplog.text
    assert 'This is unlikely to provide meaningful variant filtering unless this is a huge callset' in caplog.text, (
        caplog.text
    )

    assert matrix_out.info.AF.collect()[0] == [0.5]
    assert matrix_out.info.AC.collect()[0] == [1]
    assert matrix_out.info.AN.collect()[0] == 2


@pytest.mark.parametrize(  # needs clinvar
    'ac,af,clinvar,threshold,rows',
    [
        (1, 0.0, 0, 0.01, 1),
        (6, 0.1, 0, 0.01, 0),
        (6, 0.1, 1, 0.01, 1),
        (50, 0.001, 0, 0.01, 1),
        (50, 0.2, 0, 0.01, 0),
        (50, 0.2, 1, 0.01, 1),
    ],
)
def test_ac_filter_no_filt(
    ac: int,
    af: float,
    clinvar: int,
    threshold: float,
    rows: int,
    make_a_mt: hl.MatrixTable,
):
    """
    run tests on the ac filtering method
    """
    matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            clinvar_talos=clinvar,
            AC=[ac],
            AF=[af],
        ),
    )

    assert filter_matrix_by_ac(mt=matrix, af_threshold=threshold).count_rows() == rows


@pytest.mark.parametrize(
    'filters,clinvar,length',
    [
        (hl.empty_set(hl.tstr), 0, 1),
        (hl.literal({'fail'}), 0, 0),
        (hl.literal({'fail'}), 1, 1),
        (hl.literal({'VQSR'}), 0, 0),
        (hl.literal({'VQSR'}), 1, 1),
    ],
)
def test_filter_on_quality_flags(
    filters: hl.set,
    clinvar: int,
    length: int,
    make_a_mt: hl.MatrixTable,
):
    """
    annotate filters and run tests
    """
    # to add new alleles, we need to scrub alleles from the key fields
    anno_matrix = make_a_mt.key_rows_by('locus')
    anno_matrix = anno_matrix.annotate_rows(
        filters=filters,
        info=anno_matrix.info.annotate(clinvar_talos=clinvar),
    )
    assert filter_on_quality_flags(anno_matrix).count_rows() == length


@pytest.fixture(name='panelapp_from_locations')
def fixture_panelapp_from_locations():
    """build a minimal PanelApp object from a dict of ENSG: location"""

    def _panelapp(locations: dict[str, str]) -> PanelApp:
        return PanelApp(
            genes={
                ensg: PanelDetail(symbol=ensg, chrom=location.split(':')[0], location=location)
                for ensg, location in locations.items()
            },
        )

    return _panelapp


@pytest.mark.parametrize(
    ('location', 'expected'),
    [
        ('1:100-200', ('chr1', 100, 200)),
        ('X:100-200', ('chrX', 100, 200)),
        # Ensembl calls the mitochondrion MT, GRCh38 in Hail calls it chrM
        ('MT:100-200', ('chrM', 100, 200)),
        # already prefixed locations are left alone
        ('chr1:100-200', ('chr1', 100, 200)),
        # genes added through config have no coordinates
        ('Unknown', None),
        ('', None),
        ('1:100', None),
    ],
)
def test_parse_gene_location(location: str, expected: tuple[str, int, int] | None):
    """PanelApp locations become GRCh38 contigs and coordinates, or None when unusable"""
    assert parse_gene_location(location) == expected


def test_green_gene_intervals(panelapp_from_locations):
    """overlapping and abutting regions merge, coordinates are clamped inside the contig"""
    panelapp = panelapp_from_locations(
        {
            # flanking runs off the start of the contig, and this gene overlaps GENE2
            'ENSG01': '1:1-18000',
            'ENSG02': '1:17001-28000',
            # once flanked, this abuts the end of GENE2
            'ENSG03': '1:32001-38000',
            # separated from the ENSG01-03 block
            'ENSG04': '1:500001-600000',
            'ENSG05': '2:3001-4000',
            # flanking runs off the end of the contig
            'ENSG06': 'MT:1-16000',
        },
    )
    intervals = green_gene_intervals(panelapp)
    assert [str(interval) for interval in intervals] == [
        '[chr1:1-40000]',
        '[chr1:498001-602000]',
        '[chr2:1001-6000]',
        '[chrM:1-16569]',
    ]


def test_green_gene_intervals_skips_unusable_regions(panelapp_from_locations, caplog):
    """genes with no location, or on unknown contigs, are dropped with a warning"""
    panelapp = panelapp_from_locations(
        {
            'ENSG01': '2:3001-4000',
            'ENSG02': 'BANANA:1-100',
            'ENSG03': 'Unknown',
        },
    )
    intervals = green_gene_intervals(panelapp)
    assert [str(interval) for interval in intervals] == ['[chr2:1001-6000]']
    assert 'contig chrBANANA is not in GRCh38' in caplog.text
    assert '2 genes of interest had no usable location' in caplog.text


def test_green_gene_intervals_no_locations(panelapp_from_locations):
    """if none of the genes of interest have a location, that's fatal - we'd read an empty MT"""
    panelapp = panelapp_from_locations({'ENSG_NOPE': 'Unknown'})
    with pytest.raises(ValueError, match='None of the 1 genes of interest had a usable location'):
        green_gene_intervals(panelapp)
