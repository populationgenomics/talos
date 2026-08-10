"""
Tests for the streaming small-variant filtering process (run_streaming_filtering.py).

These port the behavioural spec from the Hail tests (test_hail_categories.py,
test_hail_filters.py, test_de_novo.py) onto the cyvcf2 implementation.
"""

import json

import pytest
from cyvcf2 import VCF

from talos.models import PanelApp, PanelDetail
from talos.run_streaming_filtering import (
    AUTOSOME_OR_PAR,
    HET,
    HOM_ALT,
    HOM_REF,
    MITO,
    X_NONPAR,
    Y_NONPAR,
    DeNovoConfig,
    any_category_assigned,
    candidate_configuration,
    category_alphamissense,
    category_avi,
    category_high_impact,
    category_spliceai,
    clinvar_category_flags,
    consequences_for_gene,
    filter_by_consequence,
    is_confident_benign,
    passes_ac_filter,
    passes_population_rare,
    pm5_matches,
    read_clinvar_fields,
    resolve_trio_entry,
    variant_region,
)
from talos.run_streaming_filtering import (
    main as streaming_main,
)
from talos.vcf_streaming import MISSING_STRING, PATHOGENIC

CRITICAL_CSQS = {'frameshift', 'splice_acceptor', 'splice_donor', 'start_lost', 'stop_gained', 'stop_lost'}
DN_RELEVANT_CSQS = CRITICAL_CSQS | {'missense', 'inframe_deletion', 'inframe_insertion'}

# the BCSQ-era csq_string, matching example_config.toml
CSQ_STRING_CONFIG = [
    'consequence',
    'gene_id',
    'gene',
    'transcript',
    'mane_id',
    'mane',
    'biotype',
    'dna_change',
    'amino_acid_change',
    'codon',
    'ensp',
    'am_class',
    'am_pathogenicity',
]

BCSQ_HEADER_LINE = (
    '##INFO=<ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq, '
    'Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">'
)

VCF_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr1,length=248956422>\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FILTER=<ID=FAIL,Description="Failed">\n'
    f'{BCSQ_HEADER_LINE}\n'
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n'
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
    '##INFO=<ID=gnomad_AC_joint,Number=1,Type=Integer,Description="gnomAD AC">\n'
    '##INFO=<ID=gnomad_AF_joint,Number=1,Type=Float,Description="gnomAD AF">\n'
    '##INFO=<ID=gnomad_AC_joint_XY,Number=1,Type=Integer,Description="gnomAD AC XY">\n'
    '##INFO=<ID=gnomad_HomAlt_joint,Number=1,Type=Integer,Description="gnomAD HomAlt">\n'
    '##INFO=<ID=clinvar_significance,Number=1,Type=String,Description="ClinvArbitration significance">\n'
    '##INFO=<ID=clinvar_stars,Number=1,Type=Integer,Description="ClinvArbitration stars">\n'
    '##INFO=<ID=clinvar_allele,Number=1,Type=Integer,Description="ClinvArbitration allele ID">\n'
    '##INFO=<ID=am_transcript,Number=1,Type=String,Description="AlphaMissense transcript">\n'
    '##INFO=<ID=am_class,Number=1,Type=String,Description="AlphaMissense class">\n'
    '##INFO=<ID=am_score,Number=1,Type=Float,Description="AlphaMissense score">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale\tfather_1\tmother_1\n'
)

BASE_INFO = 'AC=1;AN=6;AF=0.001;gnomad_AC_joint=1;gnomad_AF_joint=0.0001;gnomad_AC_joint_XY=0;gnomad_HomAlt_joint=0'
PLP_INFO = 'clinvar_significance=Pathogenic/Likely_Pathogenic;clinvar_stars=1;clinvar_allele=25'

TRIO_GTS = ('0/1:60:30', '0/0:60:30', '0/0:60:30')

TWO_EXPECTED = 2


def bcsq(
    *,
    consequence='missense',
    gene='GENE1',
    transcript='ENST01',
    biotype='protein_coding',
    aa='123P',
    dna='456A>G',
):
    return f'{consequence}|{gene}|{transcript}|{biotype}|+|{aa}|{dna}'


def make_row(info: str, *, gts=TRIO_GTS, filters='PASS', pos=12345, ref='A', alt='G', gt_format='GT:GQ:DP'):
    samples = '\t'.join(gts)
    return f'chr1\t{pos}\t.\t{ref}\t{alt}\t60\t{filters}\t{info}\t{gt_format}\t{samples}'


def fake_config_retrieve(key, default=None):
    """Stand-in config: BCSQ-era csq_string, defaults for everything else."""
    if key == ['RunHailFiltering', 'csq_string']:
        return CSQ_STRING_CONFIG
    return default


def make_pedigree(tmp_path):
    ped_path = tmp_path / 'pedigree.ped'
    ped_path.write_text(
        'family_1\tmale\tfather_1\tmother_1\t1\t2\nfamily_1\tfather_1\t0\t0\t1\t1\nfamily_1\tmother_1\t0\t0\t2\t1\n',
    )
    return str(ped_path)


def run_streaming(tmp_path, monkeypatch, rows: list[str], new_gene: bool = False, pm5: dict | None = None):
    """Run the full streaming process over constructed rows, returning the labelled output variants."""
    monkeypatch.setattr('talos.run_streaming_filtering.config_retrieve', fake_config_retrieve)
    monkeypatch.setattr('talos.vcf_streaming.config_retrieve', fake_config_retrieve)

    vcf_path = tmp_path / 'input.vcf'
    vcf_path.write_text(VCF_HEADER + ''.join(f'{row}\n' for row in rows))

    panel = PanelApp(
        genes={
            'ENSG_GREEN': PanelDetail(symbol='GENE1', chrom='1', new=[1] if new_gene else []),
            'ENSG_OTHER': PanelDetail(symbol='GENE2', chrom='1'),
        },
    )
    panel_path = tmp_path / 'panelapp.json'
    panel_path.write_text(panel.model_dump_json())

    mane_path = tmp_path / 'mane.json'
    mane_path.write_text(json.dumps({'ENST01': {'mane_status': 'MANE Select', 'ensp': 'ENSP01', 'mane_id': 'NM_001'}}))

    pm5_path = None
    if pm5 is not None:
        pm5_path = str(tmp_path / 'pm5.json')
        with open(pm5_path, 'w') as handle:
            json.dump(pm5, handle)

    output_path = str(tmp_path / 'labelled.vcf.gz')
    streaming_main(
        vcf_path=str(vcf_path),
        panelapp_path=str(panel_path),
        pedigree=make_pedigree(tmp_path),
        vcf_out=output_path,
        mane_path=str(mane_path),
        pm5=pm5_path,
    )

    return [dict(variant.INFO) for variant in VCF(output_path)]


# --- pure-function ports of the Hail unit tests ---


@pytest.mark.parametrize(
    'consequence,classified',
    [
        ('frameshift', 1),
        ('synonymous', 0),
        ('stop_gained&frameshift', 1),
    ],
)
def test_highimpact(consequence, classified):
    assert category_high_impact([{'consequence': consequence}], CRITICAL_CSQS) == classified


@pytest.mark.parametrize(
    'splice_score,classified',
    [
        (0.0, 0),
        (0.499, 0),
        (1.0, 1),
        (None, 0),
    ],
)
def test_spliceai(splice_score, classified):
    assert category_spliceai(splice_score, threshold=0.5) == classified


@pytest.mark.parametrize(
    'am_pathogenicity,classified',
    [
        (1, 1),
        (0, 0),
        (0.88, 1),
        (None, 0),
    ],
)
def test_alphamissense_assignment(am_pathogenicity, classified):
    assert category_alphamissense([{'am_pathogenicity': am_pathogenicity}], am_threshold=0.564) == classified


@pytest.mark.parametrize(
    'avi,expected',
    [
        (0.0, 0),
        (0.74, 0),
        (0.75, 1),
        (0.99, 1),
        (99.9, 1),
    ],
)
def test_annotate_avi_category(avi, expected):
    assert category_avi(avi, threshold=0.75) == expected


def test_avi_category_without_threshold():
    assert category_avi(99.9, threshold=None) == 0


@pytest.mark.parametrize(
    'genomes,clinvar,kept',
    [
        (0, 0, True),
        (0.0001, 0, True),
        (0.0001, 1, True),
        (0.1, 0, False),
        (0.1, 1, True),
        (None, 0, True),
    ],
)
def test_filter_rows_for_rare(genomes, clinvar, kept):
    assert passes_population_rare(genomes, clinvar, af_threshold=0.01) == kept


@pytest.mark.parametrize(
    'ac,af,clinvar,kept',
    [
        (1, 0.0, 0, True),
        (6, 0.1, 0, False),
        (6, 0.1, 1, True),
        (50, 0.001, 0, True),
        (50, 0.2, 0, False),
        (50, 0.2, 1, True),
    ],
)
def test_ac_filter(ac, af, clinvar, kept):
    assert passes_ac_filter({'AC': ac, 'AF': af}, clinvar, af_threshold=0.01) == kept


def test_consequences_for_gene():
    """port of test_filter_to_green_genes_and_split__consequence"""
    consequences = [
        {'gene_id': 'green', 'biotype': 'protein_coding', 'mane_id': ''},
        {'gene_id': 'green', 'biotype': 'batman', 'mane_id': 'NM_Bane'},
        {'gene_id': 'green', 'biotype': 'non_coding', 'mane_id': ''},
        {'gene_id': 'NOT_GREEN', 'biotype': 'protein_coding', 'mane_id': ''},
    ]
    retained = consequences_for_gene(consequences, 'green')
    assert len(retained) == TWO_EXPECTED
    assert {csq['biotype'] for csq in retained} == {'protein_coding', 'batman'}


UNCATEGORISED = {
    'categorybooleanclinvarplp': 0,
    'categorybooleanhighimpact': 0,
    'categorysampledenovo': 'missing',
    'categorybooleanspliceai': 0,
    'categorybooleanalphamissense': 0,
    'categorydetailspm5': 'missing',
    'categorybooleanclinvar0star': 0,
    'categorybooleanclinvar0starnewgene': 0,
    'categorybooleanavi': 0,
}


@pytest.mark.parametrize(
    'overrides,expected',
    [
        ({}, False),
        ({'categorybooleanavi': 1}, True),
        ({'categorybooleanclinvarplp': 1}, True),
        ({'categorybooleanhighimpact': 1}, True),
        ({'categorysampledenovo': 'present'}, True),
        ({'categorybooleanalphamissense': 1}, True),
        ({'categorydetailspm5': 'present'}, True),
        ({'categorybooleanclinvar0star': 1}, True),
        ({'categorybooleanclinvar0starnewgene': 1}, True),
        ({'categorybooleanspliceai': 1}, True),
    ],
)
def test_filter_to_classified(overrides, expected):
    assert any_category_assigned(UNCATEGORISED | overrides) == expected


@pytest.mark.parametrize(
    'rating,stars,dropped,regular,strong',
    [
        ('benign', 0, False, 0, 0),
        ('benign', 1, True, 0, 0),
        ('other', 7, False, 0, 0),
        (PATHOGENIC, 0, False, 1, 0),
        (PATHOGENIC, 1, False, 1, 1),
    ],
)
def test_annotate_talos_clinvar(rating, stars, dropped, regular, strong):
    """port of test_annotate_talos_clinvar"""
    assert is_confident_benign(rating, stars) == dropped
    flags = clinvar_category_flags(rating, stars, gene_ids={'ensga'}, new_genes={'ensga'})
    assert flags['clinvar_talos'] == regular
    assert flags['categorybooleanclinvarplp'] == strong


@pytest.mark.parametrize(
    'rating,stars,expected_flag',
    [
        (PATHOGENIC, 0, 1),
        (PATHOGENIC, 1, 0),
        ('other', 3, 0),
        ('benign', 0, 0),
    ],
)
def test_annotate_clinvar_0star_category(rating, stars, expected_flag):
    flags = clinvar_category_flags(rating, stars, gene_ids={'ensga'}, new_genes={'ensga'})
    assert flags['categorybooleanclinvar0star'] == expected_flag


@pytest.mark.parametrize(
    'rating,stars,new_set,expected_flag',
    [
        (PATHOGENIC, 0, {'ensga'}, 1),
        (PATHOGENIC, 0, {'NOT_MATCH'}, 0),
        (PATHOGENIC, 1, {'ensga'}, 0),
        ('other', 0, {'ensga'}, 0),
    ],
)
def test_annotate_clinvar_0star_newgene_category(rating, stars, new_set, expected_flag):
    flags = clinvar_category_flags(rating, stars, gene_ids={'ensga'}, new_genes=new_set)
    assert flags['categorybooleanclinvar0starnewgene'] == expected_flag


def test_read_clinvar_fields_absent():
    assert read_clinvar_fields({}) == (MISSING_STRING, 0, 0)


def test_read_clinvar_fields_echtvar_sentinels():
    """echtvar fills unannotated integer fields with -1, and encodes spaces as underscores"""
    info = {'clinvar_significance': 'Pathogenic/Likely_Pathogenic', 'clinvar_stars': -1, 'clinvar_allele': -1}
    assert read_clinvar_fields(info) == (PATHOGENIC, 0, 0)


def test_pm5_matches():
    lookup = {'ENST01::123': '25::2+55::1'}
    consequences = [{'consequence': 'missense', 'transcript': 'ENST01', 'codon': 123}]
    assert pm5_matches(True, consequences, lookup) == '25::2+55::1'
    # non-SNV, non-missense, unknown codon, and no-lookup all miss
    assert pm5_matches(False, consequences, lookup) == MISSING_STRING
    assert pm5_matches(True, [{'consequence': 'synonymous', 'transcript': 'ENST01', 'codon': 123}], lookup) == (
        MISSING_STRING
    )
    assert pm5_matches(True, [{'consequence': 'missense', 'transcript': 'ENST01', 'codon': None}], lookup) == (
        MISSING_STRING
    )
    assert pm5_matches(True, consequences, None) == MISSING_STRING


def test_filter_by_consequence():
    consequences = [
        {'consequence': 'frameshift', 'biotype': 'protein_coding'},
        {'consequence': 'synonymous', 'biotype': 'protein_coding'},
        {'consequence': 'synonymous', 'biotype': 'snRNA'},
    ]
    retained = filter_by_consequence(consequences, DN_RELEVANT_CSQS)
    assert len(retained) == TWO_EXPECTED
    assert not filter_by_consequence([{'consequence': 'synonymous', 'biotype': 'protein_coding'}], DN_RELEVANT_CSQS)


# --- de novo machinery ---


@pytest.mark.parametrize(
    'chrom,pos,region',
    [
        ('chr1', 12345, AUTOSOME_OR_PAR),
        ('chr22', 1, AUTOSOME_OR_PAR),
        # X PAR1 is start-inclusive, end-exclusive
        ('chrX', 10001, AUTOSOME_OR_PAR),
        ('chrX', 2781478, AUTOSOME_OR_PAR),
        ('chrX', 2781479, X_NONPAR),
        ('chrX', 5000000, X_NONPAR),
        ('chrX', 155701383, AUTOSOME_OR_PAR),
        ('chrY', 5000000, Y_NONPAR),
        ('chrY', 10001, AUTOSOME_OR_PAR),
        ('chrM', 100, MITO),
        ('chrUn_KI270742v1', 100, None),
    ],
)
def test_variant_region(chrom, pos, region):
    assert variant_region(chrom, pos) == region


@pytest.mark.parametrize(
    'region,is_female,trio_gts,expected',
    [
        (AUTOSOME_OR_PAR, False, (HET, HOM_REF, HOM_REF), True),
        (AUTOSOME_OR_PAR, False, (HOM_ALT, HOM_REF, HOM_REF), False),
        (AUTOSOME_OR_PAR, False, (HET, HET, HOM_REF), False),
        (AUTOSOME_OR_PAR, False, (HET, None, HOM_REF), False),
        (X_NONPAR, False, (HOM_ALT, HOM_REF, HOM_REF), True),
        (X_NONPAR, False, (HET, HOM_REF, HOM_REF), True),
        (X_NONPAR, True, (HOM_ALT, HOM_REF, HOM_REF), False),
        (X_NONPAR, True, (HET, HOM_REF, HOM_REF), True),
        (X_NONPAR, None, (HET, HOM_REF, HOM_REF), False),
        (Y_NONPAR, False, (HOM_ALT, HOM_REF, None), True),
        (Y_NONPAR, False, (HOM_ALT, HET, None), False),
        (MITO, False, (HOM_ALT, None, HOM_REF), True),
        (MITO, False, (HET, None, HOM_REF), False),
    ],
)
def test_candidate_configuration(region, is_female, trio_gts, expected):
    kid, dad, mom = trio_gts
    assert candidate_configuration(region, is_female, kid, dad, mom) == expected


def test_resolve_trio_entry():
    conf = DeNovoConfig()

    # a clean entry passes through
    assert resolve_trio_entry(HET, 60, 30, affected=True, conf=conf) == (HET, 60)

    # below the all-sample GQ floor, the entry is removed
    assert resolve_trio_entry(HET, 10, 30, affected=False, conf=conf) == (None, None)

    # ... unless that filter is disabled
    lax = DeNovoConfig(apply_min_all_sample_gq=False)
    assert resolve_trio_entry(HET, 10, 30, affected=False, conf=lax) == (HET, 10)

    # a missing GT with a convincing GQ is recovered as HomRef
    assert resolve_trio_entry(2, 60, 30, affected=False, conf=conf) == (HOM_REF, 60)

    # a missing GT with missing GQ stays missing
    assert resolve_trio_entry(2, -1, 30, affected=False, conf=conf) == (None, None)

    # depth out of bounds removes the entry
    assert resolve_trio_entry(HET, 60, 2, affected=True, conf=conf) == (None, None)
    assert resolve_trio_entry(HET, 60, 2000, affected=True, conf=conf) == (None, None)

    # affected members need GQ above the proband floor; unaffected members do not
    assert resolve_trio_entry(HET, 25, 30, affected=True, conf=conf) == (None, None)
    assert resolve_trio_entry(HET, 25, 30, affected=False, conf=conf) == (HET, 25)


# --- end-to-end, through main ---


def test_clinvar_plp_is_labelled_and_de_novo(tmp_path, monkeypatch):
    """a ClinVar P/LP missense in a green gene is kept, reshaped, and the trio flags a de novo"""
    output = run_streaming(tmp_path, monkeypatch, [make_row(f'BCSQ={bcsq()};{BASE_INFO};{PLP_INFO}')])
    assert len(output) == 1
    info = output[0]

    assert info['categorybooleanclinvarplp'] == 1
    assert info['clinvar_talos'] == 1
    assert info['clinvar_significance'] == PATHOGENIC
    assert info['clinvar_stars'] == 1
    assert info['clinvar_allele'] == 25
    assert info['gene_id'] == 'ENSG_GREEN'

    # missense is a de novo relevant consequence, and this trio configuration is a candidate
    assert info['categorysampledenovo'] == 'male'

    # BCSQ and the raw annotations were reshaped
    assert 'BCSQ' not in info
    assert 'gnomad_AF_joint' not in info
    assert info['gnomad_AC'] == 1
    assert pytest.approx(info['gnomad_AF'], abs=1e-6) == 0.0001

    # the csq string is rebuilt against config ordering, with MANE and gene ID folded in
    csq = info['csq'].split('|')
    csq_dict = dict(zip(CSQ_STRING_CONFIG, csq, strict=True))
    assert csq_dict['consequence'] == 'missense'
    assert csq_dict['gene_id'] == 'ENSG_GREEN'
    assert csq_dict['transcript'] == 'ENST01'
    assert csq_dict['mane_id'] == 'NM_001'
    assert csq_dict['ensp'] == 'ENSP01'
    assert csq_dict['codon'] == '123'


def test_high_impact_without_clinvar(tmp_path, monkeypatch):
    output = run_streaming(tmp_path, monkeypatch, [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO}')])
    assert len(output) == 1
    info = output[0]
    assert info['categorybooleanhighimpact'] == 1
    assert info['categorybooleanclinvarplp'] == 0
    assert info['clinvar_significance'] == MISSING_STRING
    assert info['categorysampledenovo'] == 'male'


def test_uncategorised_row_is_dropped(tmp_path, monkeypatch):
    """a synonymous variant with no ClinVar rating earns no categories"""
    assert not run_streaming(tmp_path, monkeypatch, [make_row(f'BCSQ={bcsq(consequence="synonymous")};{BASE_INFO}')])


def test_alphamissense_category(tmp_path, monkeypatch):
    """AlphaMissense at the matching transcript labels the row; the de novo gate stays closed for synonymous"""
    info_string = (
        f'BCSQ={bcsq(consequence="synonymous")};{BASE_INFO};am_transcript=ENST01;am_class=likely_pathogenic;'
        'am_score=0.99'
    )
    output = run_streaming(tmp_path, monkeypatch, [make_row(info_string)])
    assert len(output) == 1
    info = output[0]
    assert info['categorybooleanalphamissense'] == 1
    assert info['categorysampledenovo'] == MISSING_STRING
    csq_dict = dict(zip(CSQ_STRING_CONFIG, info['csq'].split('|'), strict=True))
    assert csq_dict['am_class'] == 'likely_pathogenic'
    assert 'am_transcript' not in info


def test_alphamissense_other_transcript_not_labelled(tmp_path, monkeypatch):
    info_string = (
        f'BCSQ={bcsq(consequence="synonymous")};{BASE_INFO};am_transcript=ENST99;am_class=likely_pathogenic;'
        'am_score=0.99'
    )
    assert not run_streaming(tmp_path, monkeypatch, [make_row(info_string)])


def test_common_in_gnomad_dropped_unless_clinvar(tmp_path, monkeypatch):
    common = BASE_INFO.replace('gnomad_AF_joint=0.0001', 'gnomad_AF_joint=0.05')
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{common}')]
    assert not run_streaming(tmp_path, monkeypatch, rows)

    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{common};{PLP_INFO}')]
    output = run_streaming(tmp_path, monkeypatch, rows)
    assert len(output) == 1
    assert output[0]['categorybooleanclinvarplp'] == 1


def test_filter_failed_dropped_unless_clinvar(tmp_path, monkeypatch):
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO}', filters='FAIL')]
    assert not run_streaming(tmp_path, monkeypatch, rows)

    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO};{PLP_INFO}', filters='FAIL')]
    assert len(run_streaming(tmp_path, monkeypatch, rows)) == 1


def test_common_in_callset_dropped_unless_clinvar(tmp_path, monkeypatch):
    common = BASE_INFO.replace('AC=1;AN=6;AF=0.001', 'AC=50;AN=100;AF=0.5')
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{common}')]
    assert not run_streaming(tmp_path, monkeypatch, rows)

    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{common};{PLP_INFO}')]
    assert len(run_streaming(tmp_path, monkeypatch, rows)) == 1


def test_confidently_benign_dropped(tmp_path, monkeypatch):
    benign = 'clinvar_significance=benign;clinvar_stars=2;clinvar_allele=25'
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO};{benign}')]
    assert not run_streaming(tmp_path, monkeypatch, rows)


def test_non_green_gene_dropped(tmp_path, monkeypatch):
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift", gene="NOT_A_PANEL_GENE")};{BASE_INFO}')]
    assert not run_streaming(tmp_path, monkeypatch, rows)


def test_two_green_genes_exploded(tmp_path, monkeypatch):
    """a variant in two green genes is written once per gene, with per-gene csq content"""
    two_genes = f'{bcsq(consequence="frameshift")},{bcsq(consequence="frameshift", gene="GENE2", transcript="ENST02")}'
    output = run_streaming(tmp_path, monkeypatch, [make_row(f'BCSQ={two_genes};{BASE_INFO}')])
    assert [info['gene_id'] for info in output] == ['ENSG_GREEN', 'ENSG_OTHER']
    for info in output:
        csq_dict = dict(zip(CSQ_STRING_CONFIG, info['csq'].split('|'), strict=True))
        assert csq_dict['gene_id'] == info['gene_id']


def test_pm5_annotation(tmp_path, monkeypatch):
    output = run_streaming(
        tmp_path,
        monkeypatch,
        [make_row(f'BCSQ={bcsq()};{BASE_INFO};{PLP_INFO}')],
        pm5={'ENST01::123': '25::2+55::1'},
    )
    assert len(output) == 1
    assert output[0]['categorydetailspm5'] == '25::2+55::1'


def test_clinvar_0star_newgene(tmp_path, monkeypatch):
    zero_star = 'clinvar_significance=Pathogenic/Likely_Pathogenic;clinvar_stars=0;clinvar_allele=25'
    output = run_streaming(tmp_path, monkeypatch, [make_row(f'BCSQ={bcsq()};{BASE_INFO};{zero_star}')], new_gene=True)
    assert len(output) == 1
    info = output[0]
    assert info['categorybooleanclinvarplp'] == 0
    assert info['categorybooleanclinvar0star'] == 1
    assert info['categorybooleanclinvar0starnewgene'] == 1


def test_star_allele_dropped(tmp_path, monkeypatch):
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO};{PLP_INFO}', alt='*')]
    assert not run_streaming(tmp_path, monkeypatch, rows)


def test_de_novo_requires_wt_parents(tmp_path, monkeypatch):
    """an inherited het is not flagged, but the row survives on its other category"""
    gts = ('0/1:60:30', '0/1:60:30', '0/0:60:30')
    output = run_streaming(
        tmp_path, monkeypatch, [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO}', gts=gts)]
    )
    assert len(output) == 1
    assert output[0]['categorysampledenovo'] == MISSING_STRING


def test_de_novo_low_gq_parent_blocks_call(tmp_path, monkeypatch):
    """a parent below the all-sample GQ floor removes the WT evidence"""
    gts = ('0/1:60:30', '0/0:10:30', '0/0:60:30')
    output = run_streaming(
        tmp_path, monkeypatch, [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO}', gts=gts)]
    )
    assert len(output) == 1
    assert output[0]['categorysampledenovo'] == MISSING_STRING


def test_de_novo_with_ad_fallback_depth(tmp_path, monkeypatch):
    """port of test_dn_bch_one - depth comes from summed AD when DP is absent"""
    gts = ('0/1:60:6,6', '0/0:60:30,0', '0/0:60:30,0')
    rows = [make_row(f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO}', gts=gts, gt_format='GT:GQ:AD')]
    output = run_streaming(tmp_path, monkeypatch, rows)
    assert len(output) == 1
    assert output[0]['categorysampledenovo'] == 'male'


def test_af_derived_when_absent(tmp_path, monkeypatch):
    """with AC/AN but no AF in the header, AF is derived and used for filtering"""
    header_no_af = VCF_HEADER.replace('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n', '')
    vcf_path = tmp_path / 'no_af.vcf'
    info_string = f'BCSQ={bcsq(consequence="frameshift")};{BASE_INFO}'.replace('AF=0.001;', '')
    rows = [make_row(info_string, gts=TRIO_GTS)]
    vcf_path.write_text(header_no_af + ''.join(f'{row}\n' for row in rows))

    monkeypatch.setattr('talos.run_streaming_filtering.config_retrieve', fake_config_retrieve)
    monkeypatch.setattr('talos.vcf_streaming.config_retrieve', fake_config_retrieve)

    panel = PanelApp(genes={'ENSG_GREEN': PanelDetail(symbol='GENE1', chrom='1')})
    panel_path = tmp_path / 'panelapp.json'
    panel_path.write_text(panel.model_dump_json())
    mane_path = tmp_path / 'mane.json'
    mane_path.write_text(json.dumps({}))

    output_path = str(tmp_path / 'labelled.vcf.gz')
    streaming_main(
        vcf_path=str(vcf_path),
        panelapp_path=str(panel_path),
        pedigree=make_pedigree(tmp_path),
        vcf_out=output_path,
        mane_path=str(mane_path),
    )
    output = [dict(variant.INFO) for variant in VCF(output_path)]
    assert len(output) == 1
    assert pytest.approx(output[0]['AF'], abs=1e-6) == 1 / 6


def test_no_shared_samples_raises(tmp_path, monkeypatch):
    monkeypatch.setattr('talos.run_streaming_filtering.config_retrieve', fake_config_retrieve)
    monkeypatch.setattr('talos.vcf_streaming.config_retrieve', fake_config_retrieve)

    vcf_path = tmp_path / 'input.vcf'
    vcf_path.write_text(VCF_HEADER + make_row(f'BCSQ={bcsq()};{BASE_INFO}') + '\n')

    panel = PanelApp(genes={'ENSG_GREEN': PanelDetail(symbol='GENE1', chrom='1')})
    panel_path = tmp_path / 'panelapp.json'
    panel_path.write_text(panel.model_dump_json())
    mane_path = tmp_path / 'mane.json'
    mane_path.write_text(json.dumps({}))

    ped_path = tmp_path / 'unrelated.ped'
    ped_path.write_text('family_9\tstranger\t0\t0\t1\t2\n')

    with pytest.raises(ValueError, match='No samples shared between pedigree and VCF'):
        streaming_main(
            vcf_path=str(vcf_path),
            panelapp_path=str(panel_path),
            pedigree=str(ped_path),
            vcf_out=str(tmp_path / 'labelled.vcf.gz'),
            mane_path=str(mane_path),
        )
