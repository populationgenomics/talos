"""
unit testing collection for the hail MT methods
"""

import pandas as pd
import pytest

import hail as hl

from talos.models import PanelApp
from talos.RunHailFiltering import (
    annotate_category_alphamissense,
    annotate_category_high_impact,
    annotate_category_spliceai,
    annotate_clinvarbitration,
    filter_to_categorised,
    filter_to_population_rare,
    green_from_panelapp,
    split_rows_by_gene_and_filter_to_green,
)
from test.test_utils import ONE_EXPECTED, TWO_EXPECTED

hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'clinvar_talos,consequence_terms,classified',
    [
        (0, 'frameshift', ONE_EXPECTED),
        (1, 'frameshift', ONE_EXPECTED),
    ],
)
def test_highimpact(clinvar_talos, consequence_terms, classified, make_a_mt):
    """"""

    anno_matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(clinvar_talos=clinvar_talos),
        transcript_consequences=hl.array(
            [
                hl.Struct(consequence=consequence_terms),
            ],
        ),
    )

    anno_matrix = annotate_category_high_impact(anno_matrix)
    assert anno_matrix.info.categorybooleanhighimpact.collect() == [classified]


@pytest.mark.parametrize(
    'splice_score,classified',
    [
        (0.0, 0),
        (0.499, 0),
        (1.0, ONE_EXPECTED),
    ],
)
def test_spliceai(splice_score, classified, make_a_mt):
    """"""

    anno_matrix = make_a_mt.annotate_rows(
        splice_ai=hl.Struct(delta_score=splice_score),
    )

    anno_matrix = annotate_category_spliceai(anno_matrix)
    assert anno_matrix.info.categorybooleanspliceai.collect() == [classified]


@pytest.mark.skip(reason='category 5 currently inactive')
@pytest.mark.parametrize(
    'spliceai_score,flag',
    [(0.1, 0), (0.11, 0), (0.3, 0), (0.49, 0), (0.5, 1), (0.69, 1), (0.9, 1)],
)
def test_category_5_assignment(spliceai_score: float, flag: int, make_a_mt):
    """
    we don't currently use category 5, but test the logic here anyway
    """

    matrix = make_a_mt.annotate_rows(info=make_a_mt.info.annotate(splice_ai_delta=spliceai_score))
    assert matrix.info.categorybooleanspliceai.collect() == [flag]


@pytest.mark.parametrize(
    'am_class,classified',
    [
        ('likely_pathogenic', 1),
        ('not_pathogenic', 0),
        ('', 0),
        (hl.missing('tstr'), 0),
    ],
)
def test_alphamissense_assignment(am_class, classified, make_a_mt):
    """"""

    anno_matrix = make_a_mt.annotate_rows(
        transcript_consequences=hl.array(
            [
                hl.Struct(am_class=am_class),
            ],
        ),
    )

    anno_matrix = annotate_category_alphamissense(anno_matrix)
    anno_matrix.rows().show()
    assert anno_matrix.info.categorybooleanalphamissense.collect() == [classified]


def annotate_c6_missing(make_a_mt, caplog):
    """
    test what happens if the am_class attribute is missing

    Args:
        make_a_mt ():
        caplog (fixture): pytest fixture, captures all logged text
    """
    anno_matrix = make_a_mt.annotate_rows(
        transcript_consequences=hl.array(
            [
                hl.Struct(not_am='a value'),
            ],
        ),
    )

    anno_matrix = annotate_category_alphamissense(anno_matrix)
    assert anno_matrix.info.categorybooleanalphamissense.collect() == [0]
    assert 'AlphaMissense class not found, skipping annotation' in caplog.text


def test_green_from_panelapp():
    """
    check that the set expressions from panelapp data are correct
    this is collection of ENSG names from panelapp
    """
    mendeliome_data = {
        'ENSG00ABCD': {'new': [1], 'symbol': 'ABCD'},
        'ENSG00EFGH': {'new': [], 'symbol': 'EFHG'},
        'ENSG00IJKL': {'new': [2], 'symbol': 'IJKL'},
    }
    mendeliome = PanelApp.model_validate({'genes': mendeliome_data})
    green_expression = green_from_panelapp(mendeliome)

    # check types
    assert isinstance(green_expression, hl.SetExpression)

    # check content by collecting
    assert sorted(green_expression.collect()[0]) == ['ENSG00ABCD', 'ENSG00EFGH', 'ENSG00IJKL']


@pytest.mark.parametrize(
    'genomes,clinvar,length',
    [
        [0, 0, 1],
        [0.0001, 0, 1],
        [0.0001, 1, 1],
    ],
)
def test_filter_rows_for_rare(
    genomes: float,
    clinvar: int,
    length: int,
    make_a_mt,
):
    """
    annotate categories and test for retention
    """
    anno_matrix = make_a_mt.annotate_rows(
        gnomad=hl.Struct(gnomad_AF=genomes),
        info=make_a_mt.info.annotate(clinvar_talos=clinvar),
    )
    matrix = filter_to_population_rare(anno_matrix)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'gene_ids,length',
    [
        ({'not_green'}, 0),
        ({'green'}, 1),
        ({'gene'}, 1),
        ({'gene', 'not_green'}, 1),
        ({'green', 'gene'}, 2),
        ({hl.missing(t=hl.tstr)}, 0),
    ],
)
def test_filter_to_green_genes_and_split(gene_ids, length, make_a_mt):
    """

    Args:
        gene_ids ():
        length ():
        make_a_mt ():
    """
    green_genes = hl.literal({'green', 'gene'})
    anno_matrix = make_a_mt.annotate_rows(
        gene_ids=hl.literal(gene_ids),
        transcript_consequences=hl.array([hl.Struct(gene_id='gene', biotype='protein_coding', mane_id='')]),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == length


def test_filter_to_green_genes_and_split__consequence(make_a_mt):
    """

    Args:
        make_a_mt ():
    """

    green_genes = hl.literal({'green'})
    anno_matrix = make_a_mt.annotate_rows(
        gene_ids=green_genes,
        transcript_consequences=hl.array(
            [
                hl.Struct(gene_id='green', biotype='protein_coding', mane_id=''),
                hl.Struct(gene_id='green', biotype='batman', mane_id='NM_Bane'),
                hl.Struct(gene_id='green', biotype='non_coding', mane_id=''),
                hl.Struct(gene_id='NOT_GREEN', biotype='protein_coding', mane_id=''),
            ],
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == 1
    matrix = matrix.filter_rows(hl.len(matrix.transcript_consequences) == TWO_EXPECTED)
    assert matrix.count_rows() == 1


@pytest.mark.parametrize(
    'one,three,four,five,six,pm5,svdb,exomiser,zerostar,newgene,length',
    [
        (0, 0, 'missing', 0, 0, 'missing', 0, 'missing', 0, 0, 0),
        (0, 0, 'missing', 0, 0, 'missing', 0, 'present', 0, 0, 1),
        (1, 0, 'missing', 0, 0, 'missing', 0, 'missing', 0, 0, 1),
        (0, 1, 'missing', 0, 0, 'missing', 0, 'missing', 0, 0, 1),
        (0, 0, 'present', 0, 0, 'missing', 0, 'missing', 0, 0, 1),
        (0, 0, 'missing', 0, 1, 'missing', 0, 'missing', 0, 0, 1),
        (0, 0, 'missing', 0, 0, 'present', 0, 'missing', 0, 0, 1),
        (0, 0, 'missing', 0, 0, 'missing', 1, 'missing', 0, 0, 1),
        (0, 0, 'missing', 0, 0, 'missing', 0, 'missing', 1, 0, 1),
        (0, 0, 'missing', 0, 0, 'missing', 0, 'missing', 0, 1, 1),
    ],
)
def test_filter_to_classified(
    one: int,
    three: int,
    four: str,
    five: int,
    six: int,
    pm5: str,
    svdb: int,
    exomiser: str,
    zerostar: int,
    newgene: int,
    length: int,
    make_a_mt: hl.MatrixTable,  # via a pytest fixture
):
    """
    one argument per category
    """
    anno_matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            categorybooleanclinvarplp=one,
            categorybooleanhighimpact=three,
            categorysampledenovo=four,
            categorybooleanspliceai=five,
            categorybooleanalphamissense=six,
            categorydetailspm5=pm5,
            categorybooleansvdb=svdb,
            categorydetailsexomiser=exomiser,
            categorybooleanclinvar0star=zerostar,
            categorybooleanclinvar0starnewgene=newgene,
        ),
    )
    matrix = filter_to_categorised(anno_matrix)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'rating,stars,rows,regular,strong',
    [
        ('benign', 0, 1, 0, 0),
        ('benign', 1, 0, 0, 0),
        ('other', 7, 1, 0, 0),
        ('Pathogenic/Likely Pathogenic', 0, 1, 1, 0),
        ('Pathogenic/Likely Pathogenic', 1, 1, 1, 1),
    ],
)
def test_annotate_talos_clinvar(rating, stars, rows, regular, strong, tmp_path, make_a_mt):
    """
    Test intention
    - take a VCF of two variants w/default clinvar annotations
    - create a single variant annotation table with each run
    - apply the parametrized annotations to the table
    """

    table = hl.Table.from_pandas(
        pd.DataFrame(
            [
                {
                    'locus': hl.Locus(contig='chr1', position=12345),
                    'alleles': ['A', 'G'],
                    'clinical_significance': rating,
                    'gold_stars': stars,
                    'allele_id': 1,
                },
            ],
        ),
        key=['locus', 'alleles'],
    )

    table_path = str(tmp_path / 'anno.ht')
    table.write(table_path)

    returned_table = annotate_clinvarbitration(make_a_mt, clinvar=table_path, new_genes=hl.literal({'ensga'}))
    assert returned_table.count_rows() == rows
    assert len([x for x in returned_table.info.clinvar_talos.collect() if x == 1]) == regular
    assert len([x for x in returned_table.info.categorybooleanclinvarplp.collect() if x == 1]) == strong


@pytest.mark.parametrize(
    'rating,stars,expected_flag',
    [
        ('Pathogenic/Likely Pathogenic', 0, 1),
        ('Pathogenic/Likely Pathogenic', 1, 0),
        ('other', 3, 0),
        ('benign', 0, 0),  # note: benign with stars>0 would be filtered out entirely
    ],
)
def test_annotate_clinvar_0star_category(rating, stars, expected_flag, tmp_path, make_a_mt):
    """Verify categorybooleanclinvar0star is set only for P/LP with 0 stars."""

    table = hl.Table.from_pandas(
        pd.DataFrame(
            [
                {
                    'locus': hl.Locus(contig='chr1', position=12345),
                    'alleles': ['A', 'G'],
                    'clinical_significance': rating,
                    'gold_stars': stars,
                    'allele_id': 1,
                },
            ],
        ),
        key=['locus', 'alleles'],
    )
    table_path = str(tmp_path / 'anno_0star.ht')
    table.write(table_path)

    returned_table = annotate_clinvarbitration(make_a_mt, clinvar=table_path, new_genes=hl.literal({'ensga'}))
    assert returned_table.info.categorybooleanclinvar0star.collect() == [expected_flag]


@pytest.mark.parametrize(
    'rating,stars,new_set,expected_flag',
    [
        ('Pathogenic/Likely Pathogenic', 0, {'ensga'}, 1),  # matches fixture gene and is 0-star
        ('Pathogenic/Likely Pathogenic', 0, {'NOT_MATCH'}, 0),
        ('Pathogenic/Likely Pathogenic', 1, {'ensga'}, 0),  # not 0-star
        ('other', 0, {'ensga'}, 0),
    ],
)
def test_annotate_clinvar_0star_newgene_category(rating, stars, new_set, expected_flag, tmp_path, make_a_mt):
    """Verify categorybooleanclinvar0starnewgene requires P/LP, 0 stars, and gene in the provided new_genes set."""

    table = hl.Table.from_pandas(
        pd.DataFrame(
            [
                {
                    'locus': hl.Locus(contig='chr1', position=12345),
                    'alleles': ['A', 'G'],
                    'clinical_significance': rating,
                    'gold_stars': stars,
                    'allele_id': 1,
                },
            ],
        ),
        key=['locus', 'alleles'],
    )
    table_path = str(tmp_path / 'anno_0star_newgene.ht')
    table.write(table_path)

    returned_table = annotate_clinvarbitration(make_a_mt, clinvar=table_path, new_genes=hl.literal(set(new_set)))
    assert returned_table.info.categorybooleanclinvar0starnewgene.collect() == [expected_flag]
