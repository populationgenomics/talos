"""
test class for the utils collection
"""

from copy import deepcopy

from cyvcf2 import VCFReader
from mendelbrot.pedigree_parser import PedigreeParser

from talos.models import (
    Coordinates,
    ReportVariant,
    ResultData,
    SmallVariant,
)
from talos.static_values import get_granular_date
from talos.utils import (
    annotate_variant_dates_using_prior_results,
    find_comp_hets,
    gather_gene_dict_from_contig,
    get_non_ref_samples,
)

ZERO_EXPECTED = 0
ONE_EXPECTED = 1
TWO_EXPECTED = 2
THREE_EXPECTED = 3
FOUR_EXPECTED = 4
FIVE_EXPECTED = 5


class FakePanelApp:
    def __init__(self):
        self.str_genes: set[str] = set()


def test_coord_sorting():
    """
    check that coord sorting methods work
    """
    coord_1 = Coordinates(chrom='4', pos=20, ref='A', alt='C')
    coord_1b = Coordinates(chrom='4', pos=21, ref='A', alt='C')
    coord_1c = Coordinates(chrom='4', pos=21, ref='A', alt='C')
    coord_2 = Coordinates(chrom='5', pos=20, ref='A', alt='C')
    assert coord_1 < coord_2
    assert coord_1 < coord_1b
    assert not coord_1b < coord_1c


def test_abs_var_sorting(two_trio_abs_variants: list[SmallVariant]):
    """
    test sorting and equivalence at the AbsVar level
    """

    var1, var2 = two_trio_abs_variants
    assert var1 < var2
    assert sorted([var2, var1]) == [var1, var2]
    # not sure if I should be able to just override the chrom...
    var1.coordinates.chrom = 'HLA1234'
    assert var1 > var2


def test_reported_variant_ordering(trio_abs_variant: SmallVariant):
    """
    test that equivalence between Report objects works as exp.
    """
    report_1 = ReportVariant(
        sample='1',
        family='1',
        gene='2',
        var_data=deepcopy(trio_abs_variant),
        reasons='test',
        genotypes={},
    )
    report_2 = ReportVariant(
        sample='1',
        family='1',
        gene='2',
        var_data=deepcopy(trio_abs_variant),
        reasons='test',
        genotypes={},
    )
    assert report_1 == report_2
    # alter sample ID, expected mismatch
    report_1.sample = '2'
    assert report_1 != report_2
    report_2.sample = '2'
    report_1.var_data.coordinates.chrom = '1'
    report_2.var_data.coordinates.chrom = '11'
    assert report_1 < report_2


def test_get_non_ref_samples(cyvcf_example_variant):
    """
    this simple test can be done without the use of a cyvcf2 object
    :return:
    """

    samples = ['male', 'father', 'mother']
    het, hom = get_non_ref_samples(variant=cyvcf_example_variant, samples=samples)
    assert het == {'male'}
    assert not hom


def test_av_categories(trio_abs_variant: SmallVariant):
    """
    Cat. 3, and Cat. 4 for PROBAND only:
    """
    assert trio_abs_variant.info.get('categoryboolean3')
    assert not trio_abs_variant.info.get('categoryboolean1')
    assert not trio_abs_variant.info.get('categoryboolean2')
    assert trio_abs_variant.sample_category_check('male')

    for sample_cat in trio_abs_variant.sample_categories:
        assert 'father_1' not in trio_abs_variant.info[sample_cat]


def test_av_categories_support(trio_abs_variant: SmallVariant):
    """
    Cat. 3, and Cat. 4 for PROBAND only
    """
    assert trio_abs_variant.info.get('categoryboolean3')
    assert trio_abs_variant.sample_category_check('male')

    # now make the categories support-only
    trio_abs_variant.support_categories.update({'3', '4'})
    assert trio_abs_variant.sample_category_check('male')
    assert not trio_abs_variant.sample_category_check('male', allow_support=False)


def test_av_phase(trio_abs_variant: SmallVariant):
    """
    nothing here yet
    """
    assert trio_abs_variant.phased == {}


def test_gene_dict(two_trio_variants_vcf):
    """
    gene = ENSG00000075043
    """
    reader = VCFReader(two_trio_variants_vcf)
    var_dict = gather_gene_dict_from_contig(contig='chr20', variant_sources={'small': reader}, panelapp=FakePanelApp())
    assert len(var_dict) == 1
    assert 'ENSG00000075043' in var_dict
    assert len(var_dict['ENSG00000075043']) == TWO_EXPECTED


def test_comp_hets(two_trio_abs_variants: list[SmallVariant], pedigree_path):
    """
    {
        'male': {
            '20-63406931-C-CGG': [Variant()],
            '20-63406991-C-CGG': [Variant()]
        }
    }
    """
    ch_dict = find_comp_hets(two_trio_abs_variants, pedigree=PedigreeParser(pedigree_path))
    assert 'male' in ch_dict
    results = ch_dict.get('male')
    assert isinstance(results, dict)
    assert len(results) == TWO_EXPECTED
    key_1, key_2 = list(results.keys())
    assert results[key_1][0].coordinates.string_format == key_2
    assert results[key_2][0].coordinates.string_format == key_1


def test_phased_dict(phased_vcf_path):
    """
    gene = ENSG00000075043
    """
    reader = VCFReader(phased_vcf_path)
    var_dict = gather_gene_dict_from_contig(contig='chr20', variant_sources={'small': reader}, panelapp=FakePanelApp())
    assert len(var_dict) == ONE_EXPECTED
    assert 'ENSG00000075043' in var_dict
    assert len(var_dict['ENSG00000075043']) == TWO_EXPECTED
    var_pair = var_dict['ENSG00000075043']
    for variant in var_pair:
        assert 'mother_1' in variant.phased
        assert variant.phased['mother_1'] == {420: '0|1'}


def test_phased_comp_hets(phased_variants: list[SmallVariant], pedigree_path: str):
    """
    phased variants shouldn't form a comp-het
    'mother_1' is het for both variants, but phase-set is same for both
    """
    ch_dict = find_comp_hets(phased_variants, pedigree=PedigreeParser(pedigree_path))
    assert len(ch_dict) == ZERO_EXPECTED


OLD_DATE = '2020-01-01'
DATE_VAR = SmallVariant(
    coordinates=Coordinates(chrom='1', pos=100, ref='A', alt='C'),
    info={},
    transcript_consequences=[],
)


def _single_variant_results(
    max_confidence: int,
    date: str,
    newly_green: str = '',
    **kwargs: str | bool,
) -> ResultData:
    """Build a ResultData with one sample carrying one variant, tagged with a single category on the given date."""
    variant = ReportVariant(
        sample='sam1',
        var_data=DATE_VAR,
        gene='ENSG1',
        categories={'2': date},
        max_confidence=max_confidence,
        first_tagged=date,
        evidence_last_updated=date,
        newly_green_date=newly_green,
        **kwargs,
    )
    return ResultData(
        results={'sam1': {'metadata': {'ext_id': 'sam1', 'family_id': 'fam1'}, 'variants': [variant]}},
    )


def test_annotate_dates_confidence_increase():
    """A genuine jump in panel confidence flags the variant and re-dates the evidence to today."""
    today = get_granular_date()
    old = _single_variant_results(max_confidence=2, date=OLD_DATE)
    new = _single_variant_results(max_confidence=3, date=today)

    annotate_variant_dates_using_prior_results(new, old)

    variant = new.results['sam1'].variants[0]
    assert variant.confidence_increase
    assert variant.evidence_last_updated == today
    # the category itself was seen in the old run, but a transition to green is marked
    assert variant.first_tagged == today
    assert variant.newly_green_date == today


def test_newly_green_date_persists():
    """A genuine jump in panel confidence flags the variant and re-dates the evidence to today."""
    today = get_granular_date()
    between_super_old_and_now = '2022-01-01'
    old = _single_variant_results(max_confidence=3, date=OLD_DATE, newly_green=between_super_old_and_now)

    new = _single_variant_results(max_confidence=3, date=today)

    annotate_variant_dates_using_prior_results(new, old)

    variant = new.results['sam1'].variants[0]
    assert variant.evidence_last_updated == between_super_old_and_now
    assert variant.first_tagged == between_super_old_and_now


def test_annotate_dates_confidence_placeholder_ignored():
    """
    -1 is the liftover placeholder for results which pre-date confidence tracking
    it must not count as an increase, and must not short-circuit the rest of the date carry-forward
    """
    today = get_granular_date()
    old = _single_variant_results(max_confidence=-1, date=OLD_DATE, date_of_phenotype_match=OLD_DATE)
    new = _single_variant_results(max_confidence=3, date=today)

    annotate_variant_dates_using_prior_results(new, old)

    variant = new.results['sam1'].variants[0]
    assert not variant.confidence_increase
    assert variant.evidence_last_updated == OLD_DATE
    assert variant.first_tagged == OLD_DATE
    # this carry-forward sits after the confidence check, and was previously skipped for placeholder values
    assert variant.date_of_phenotype_match == OLD_DATE


def test_annotate_dates_confidence_unchanged_or_lower():
    """Same or lower confidence is not an increase, and leaves the evidence date alone."""
    today = get_granular_date()
    for new_confidence in (3, 2):
        old = _single_variant_results(max_confidence=3, date=OLD_DATE)
        new = _single_variant_results(max_confidence=new_confidence, date=today)

        annotate_variant_dates_using_prior_results(new, old)

        variant = new.results['sam1'].variants[0]
        assert not variant.confidence_increase
        assert variant.evidence_last_updated == OLD_DATE


def test_annotate_dates_confidence_reset_when_not_found():
    """A variant absent from the current run is carried over, with any prior increase flag cleared."""
    old = _single_variant_results(max_confidence=3, date=OLD_DATE, confidence_increase=True)
    new = ResultData(results={'sam1': {'metadata': {'ext_id': 'sam1', 'family_id': 'fam1'}, 'variants': []}})

    annotate_variant_dates_using_prior_results(new, old)

    assert len(new.results['sam1'].variants) == ONE_EXPECTED
    variant = new.results['sam1'].variants[0]
    assert not variant.found_in_current_run
    assert not variant.confidence_increase
