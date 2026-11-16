from pathlib import Path

import pytest

from talos import config
from talos.hpo_flagging import HpoOntology, annotate_phenotype_matches, find_genes_in_these_results
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

# All
# ├── Phenotypic abnormality
# │   ├── Abnormality of the nervous system > Abnormal nervous system physiology > Seizure (alt_id HP:0001275)
# │   │                                                                           ├── Focal-onset seizure
# │   │                                                                           └── Bilateral tonic-clonic seizure
# │   └── Abnormality of the eye > Visual impairment
# └── Mode of inheritance > Autosomal dominant inheritance
# plus obsolete HP:0009999, replaced by Seizure
OBO = str(Path(__file__).parent / 'input' / 'hpo_phenotype_match.obo')

SEIZURE = 'HP:0001250'
FOCAL = 'HP:0007359'
TONIC_CLONIC = 'HP:0002069'
VISUAL = 'HP:0000505'
AD_INHERITANCE = 'HP:0000006'

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


@pytest.fixture(name='ontology', scope='module')
def fixture_ontology() -> HpoOntology:
    return HpoOntology(OBO)


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


def test_ontology_counts_phenotypic_abnormality_subtree(ontology):
    assert ontology.descendant_count[SEIZURE] == 3
    assert ontology.descendant_count['HP:0000118'] == 8
    assert AD_INHERITANCE not in ontology.ancestors


def test_ontology_filters_and_redirects(ontology):
    """inheritance and the root are dropped, alt_ids and obsolete IDs with a replacement are redirected"""
    terms = {AD_INHERITANCE, 'HP:0000118', 'HP:0001275', 'HP:0009999', VISUAL, 'HP:1234567'}
    assert ontology.phenotypic_abnormalities(terms) == {SEIZURE, VISUAL}


@pytest.fixture
def hpo_config(monkeypatch):
    """swap in a config holding only an HPOFlagging section - note: no max_ancestor_descendants unless given"""

    def _set(semantic_match: bool, **extra: int) -> None:
        monkeypatch.setattr(config, '_config', {'HPOFlagging': {'semantic_match': semantic_match, **extra}})

    return _set


def one_variant_result(*hpos: str, date_of_phenotype_match: str | None = None) -> tuple[ResultData, ReportVariant]:
    """a single participant with these HPO terms, carrying one variant in ENSG1"""
    variant = ReportVariant(sample='male', var_data=SMALL_1, date_of_phenotype_match=date_of_phenotype_match)
    meta = ParticipantMeta(
        ext_id='male',
        family_id='family_1',
        phenotypes=[HpoTerm(id=hpo, label='participant label') for hpo in hpos],
    )
    result = ResultData(results={'male': ParticipantResults(metadata=meta, variants=[variant])})
    return result, variant


def test_exact_match_labelled_and_dated(hpo_config, ontology):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result(SEIZURE, VISUAL)
    annotate_phenotype_matches(result, {'ENSG1': {SEIZURE, FOCAL}}, ontology)
    assert variant.phenotype_labels == {'HP:0001250: Seizure'}
    assert variant.date_of_phenotype_match == get_granular_date()


def test_exact_match_keeps_earlier_date(hpo_config, ontology):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result(SEIZURE, date_of_phenotype_match='2020-01-01')
    annotate_phenotype_matches(result, {'ENSG1': {SEIZURE}}, ontology)
    assert variant.phenotype_labels == {'HP:0001250: Seizure'}
    assert variant.date_of_phenotype_match == '2020-01-01'


def test_no_match_no_labels_no_date(hpo_config, ontology):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result(VISUAL)
    annotate_phenotype_matches(result, {'ENSG1': {SEIZURE}}, ontology)
    assert variant.phenotype_labels == set()
    assert variant.date_of_phenotype_match is None


def test_exact_match_via_alt_id(hpo_config, ontology):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result('HP:0001275')
    annotate_phenotype_matches(result, {'ENSG1': {SEIZURE}}, ontology)
    assert variant.phenotype_labels == {'HP:0001250: Seizure'}


@pytest.mark.parametrize('semantic_match', [False, True])
def test_inheritance_never_matches(hpo_config, ontology, semantic_match):
    hpo_config(semantic_match=semantic_match)
    result, variant = one_variant_result(AD_INHERITANCE)
    annotate_phenotype_matches(result, {'ENSG1': {AD_INHERITANCE}}, ontology)
    assert variant.phenotype_labels == set()
    assert variant.date_of_phenotype_match is None


def test_semantic_match_labels_gene_term(hpo_config, ontology):
    """siblings under Seizure (3 descendants) match at the default threshold"""
    hpo_config(semantic_match=True)
    result, variant = one_variant_result(TONIC_CLONIC)
    annotate_phenotype_matches(result, {'ENSG1': {FOCAL, VISUAL}}, ontology)
    assert variant.phenotype_labels == {'HP:0007359: Focal-onset seizure'}
    assert variant.date_of_phenotype_match == get_granular_date()


def test_semantic_match_threshold(hpo_config, ontology):
    """Seizure has 3 descendants, so a threshold of 2 rejects the sibling match"""
    hpo_config(semantic_match=True, max_ancestor_descendants=2)
    result, variant = one_variant_result(TONIC_CLONIC)
    annotate_phenotype_matches(result, {'ENSG1': {FOCAL}}, ontology)
    assert variant.phenotype_labels == set()


def test_semantic_match_off(hpo_config, ontology):
    hpo_config(semantic_match=False)
    result, variant = one_variant_result(TONIC_CLONIC)
    annotate_phenotype_matches(result, {'ENSG1': {FOCAL}}, ontology)
    assert variant.phenotype_labels == set()


def test_semantic_match_only_for_unmatched_terms(hpo_config, ontology):
    """an exact match isn't joined by a semantic one for the same participant term"""
    hpo_config(semantic_match=True)
    result, variant = one_variant_result(SEIZURE)
    annotate_phenotype_matches(result, {'ENSG1': {SEIZURE, FOCAL}}, ontology)
    assert variant.phenotype_labels == {'HP:0001250: Seizure'}
