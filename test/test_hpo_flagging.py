import pytest

from talos import config, hpo_flagging
from talos.hpo_flagging import annotate_phenotype_matches, find_genes_in_these_results
from talos.models import (
    Coordinates,
    HpoTerm,
    ParticipantMeta,
    ParticipantResults,
    ReportVariant,
    ResultData,
    SmallVariant,
)
from talos.static_values import get_granular_date

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
SMALL_1 = SmallVariant(
    info={'gene_id': 'ENSG1'},
    coordinates=TEST_COORDS,
    transcript_consequences=[],
)
SV_1 = SmallVariant(
    info={'gene_id': 'ENSG2', 'lof': 'ENSG3,ENSG4'},
    coordinates=TEST_COORDS,
    transcript_consequences=[],
)


def test_find_genes_in_these_results():
    pm = ParticipantMeta(ext_id='male', family_id='family_1')
    test_results = ResultData(
        results={
            'male': ParticipantResults(metadata=pm, variants=[ReportVariant(sample='male', var_data=SMALL_1)]),
            'female': ParticipantResults(metadata=pm, variants=[ReportVariant(sample='male', var_data=SV_1)]),
        },
    )
    identified = find_genes_in_these_results(test_results)
    assert identified == {'ENSG1', 'ENSG2', 'ENSG3', 'ENSG4'}


@pytest.fixture
def hpo_config(monkeypatch):
    """swap in a config holding only an HPOFlagging section - note: no min_similarity"""

    def _set(semantic_match: bool) -> None:
        monkeypatch.setattr(config, '_config', {'HPOFlagging': {'semantic_match': semantic_match}})

    return _set


def one_variant_result(date_of_phenotype_match: str | None = None) -> tuple[ResultData, ReportVariant]:
    """a single participant with HP:1 and HP:2, carrying one variant in ENSG1"""
    variant = ReportVariant(sample='male', var_data=SMALL_1, date_of_phenotype_match=date_of_phenotype_match)
    meta = ParticipantMeta(
        ext_id='male',
        family_id='family_1',
        phenotypes=[HpoTerm(id='HP:1', label='one'), HpoTerm(id='HP:2', label='two')],
    )
    result = ResultData(results={'male': ParticipantResults(metadata=meta, variants=[variant])})
    return result, variant


def test_exact_match_labelled_and_dated(hpo_config):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result()
    annotate_phenotype_matches(result, {'ENSG1': {'HP:1', 'HP:9'}})
    assert variant.phenotype_labels == {'HP:1: one'}
    assert variant.date_of_phenotype_match == get_granular_date()


def test_exact_match_keeps_earlier_date(hpo_config):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result(date_of_phenotype_match='2020-01-01')
    annotate_phenotype_matches(result, {'ENSG1': {'HP:1'}})
    assert variant.phenotype_labels == {'HP:1: one'}
    assert variant.date_of_phenotype_match == '2020-01-01'


def test_no_match_no_labels_no_date(hpo_config):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result()
    annotate_phenotype_matches(result, {'ENSG1': {'HP:9'}})
    assert variant.phenotype_labels == set()
    assert variant.date_of_phenotype_match is None


def test_semantic_matches_add_to_exact(hpo_config, monkeypatch):
    """semantic matches are unioned with exact ones; min_similarity absent from config falls back to the default"""
    hpo_config(semantic_match=True)
    seen_similarity: list[float] = []

    def fake_semantic(_participant: set[str], _gene: set[str], min_similarity: float) -> set[str]:
        seen_similarity.append(min_similarity)
        return {'HP:5: five'}

    monkeypatch.setattr(hpo_flagging, 'semantic_matches', fake_semantic)
    result, variant = one_variant_result()
    annotate_phenotype_matches(result, {'ENSG1': {'HP:1', 'HP:5'}})
    assert variant.phenotype_labels == {'HP:1: one', 'HP:5: five'}
    assert variant.date_of_phenotype_match == get_granular_date()
    assert seen_similarity == [hpo_flagging.DEFAULT_MIN_SIMILARITY]
