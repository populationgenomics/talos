"""
Unit tests for the filtering helpers in talos.run_sv_filtering.

These are the ported spec of the Hail-based tests this streaming implementation replaces.
They cover the field rearrangement, the two blanket frequency filters, the hemizygous-call
fix, and the end-to-end labelling process.
"""

import pytest
from cyvcf2 import VCF
from loguru import logger

from talos.models import GeneDetail, PanelApp, RegionDetail
from talos.run_sv_filtering import (
    diploidise_genotypes,
    matching_regions,
    passes_af_filter,
    passes_callset_af_filter,
    rearrange_annotations,
)
from talos.run_sv_filtering import (
    main as sv_main,
)


def _annotated_sv_info(use_af_male_spelling: bool = False) -> dict:
    """
    build an INFO dict carrying every field rearrange_annotations reads unconditionally

    Args:
        use_af_male_spelling (bool): if True, provide AF_MALE/AF_FEMALE instead of MALE_AF/FEMALE_AF
    """
    info = {
        'PREDICTED_LOF': 'GENE1',
        'SVTYPE': 'DEL',
        'SVLEN': 5000,
        'END': 6000,
        'AC': 1,
        'AF': 0.01,
        'AN': 100,
        'N_HET': 1,
        'N_HOMALT': 0,
        'gnomad_v4.1_sv_AF': 0.001,
        'gnomad_v4.1_sv_SVID': 'gnomad-del-1',
    }
    if use_af_male_spelling:
        info['AF_MALE'] = (0.03,)
        info['AF_FEMALE'] = (0.04,)
    else:
        info['MALE_AF'] = (0.03,)
        info['FEMALE_AF'] = (0.04,)
    return info


def test_rearrange_surfaces_gnomad_and_defaults():
    """the population fields are surfaced under the names Talos reads, and absent SV fields are defaulted"""
    out = rearrange_annotations(_annotated_sv_info(), {'GENE1': 'ENSG1'}, uses_af_male_spelling=False)

    assert out['gnomad_sv_AF'] == 0.001
    assert out['gnomad_sv_ID'] == 'gnomad-del-1'
    assert out['svtype'] == 'DEL'

    # ALGORITHMS was absent, so it defaults to gCNV; STATUS/CHR2/END2 are inserted as null
    assert out['algorithms'] == 'gCNV'
    assert out['status'] is None
    assert out['chr2'] is None
    assert out['end2'] is None


def test_rearrange_maps_symbols_to_ensg():
    """PREDICTED_LOF symbols are mapped to gene IDs, and unknown symbols pass through unchanged"""
    out = rearrange_annotations(_annotated_sv_info(), {'GENE1': 'ENSG1'}, uses_af_male_spelling=False)
    assert out['lof'] == {'GENE1'}
    assert out['lof_ensg'] == {'ENSG1'}
    assert out['gene_id'] == {'ENSG1'}


def test_rearrange_symbol_without_mapping_is_kept():
    """a symbol missing from the mapping is retained as itself, not dropped"""
    out = rearrange_annotations(_annotated_sv_info(), {'OTHER': 'ENSG9'}, uses_af_male_spelling=False)
    assert out['lof_ensg'] == {'GENE1'}


def test_rearrange_multiple_lof_genes():
    """a comma-joined PREDICTED_LOF string becomes a set, each symbol independently mapped"""
    info = _annotated_sv_info() | {'PREDICTED_LOF': 'GENE1,GENE2'}
    out = rearrange_annotations(info, {'GENE1': 'ENSG1'}, uses_af_male_spelling=False)
    assert out['lof'] == {'GENE1', 'GENE2'}
    assert out['gene_id'] == {'ENSG1', 'GENE2'}


def test_rearrange_reads_male_af():
    """MALE_AF/FEMALE_AF are read through as male_af/female_af"""
    out = rearrange_annotations(_annotated_sv_info(), {'GENE1': 'ENSG1'}, uses_af_male_spelling=False)
    assert out['male_af'] == (0.03,)
    assert out['female_af'] == (0.04,)


def test_rearrange_accepts_af_male_alternative():
    """the AF_MALE/AF_FEMALE spelling is accepted as an alternative to MALE_AF/FEMALE_AF"""
    out = rearrange_annotations(
        _annotated_sv_info(use_af_male_spelling=True),
        {'GENE1': 'ENSG1'},
        uses_af_male_spelling=True,
    )
    assert out['male_af'] == (0.03,)
    assert out['female_af'] == (0.04,)


def test_af_filter_removes_common_but_keeps_missing():
    """AF above threshold is removed; a null AF is treated as 0 and survives"""
    assert not passes_af_filter({'gnomad_sv_AF': 0.5}, af_threshold=0.03)
    assert passes_af_filter({'gnomad_sv_AF': None}, af_threshold=0.03)
    assert passes_af_filter({'gnomad_sv_AF': 0.001}, af_threshold=0.03)


def test_callset_af_filter_needs_both_sexes_below_threshold():
    """a variant survives only when both male and female callset AFs are at or below the threshold"""
    assert passes_callset_af_filter({'male_af': (0.01,), 'female_af': (0.02,)}, ac_threshold=0.03)
    assert not passes_callset_af_filter({'male_af': (0.9,), 'female_af': (0.02,)}, ac_threshold=0.03)
    assert not passes_callset_af_filter({'male_af': (0.01,), 'female_af': (0.9,)}, ac_threshold=0.03)
    # a missing frequency fails, matching the Hail filter this replaces
    assert not passes_callset_af_filter({'male_af': None, 'female_af': (0.02,)}, ac_threshold=0.03)
    # overall AF fallback, and no callset frequency at all
    assert passes_callset_af_filter({'AF': 0.01}, ac_threshold=0.03, af_fields=('AF',))
    assert not passes_callset_af_filter({'AF': 0.9}, ac_threshold=0.03, af_fields=('AF',))
    assert passes_callset_af_filter({}, ac_threshold=0.03, af_fields=())


def test_diploidise_genotypes():
    """haploid calls are recast biallelic (alt -> 1/1, ref -> 0/0); diploid and missing calls are untouched"""
    genotypes, modified = diploidise_genotypes(
        [
            [1, False],  # haploid alt
            [0, False],  # haploid ref
            [0, 1, False],  # diploid het
            [-1, False],  # missing haploid
        ],
    )
    assert modified
    assert genotypes == [
        [1, 1, False],
        [0, 0, False],
        [0, 1, False],
        [-1, False],
    ]

    unchanged, modified = diploidise_genotypes([[0, 1, False], [1, 1, True]])
    assert not modified
    assert unchanged == [[0, 1, False], [1, 1, True]]


def test_matching_regions():
    """a CNV matches a region if the type matches, and it overlaps the required % of the region, or is inside it and
    covers at least regulatory_region_percentage of it
    """
    regions = [
        RegionDetail(name='LOSS', chrom='1', start=1000, end=2000, variant_type='cnv_loss', required_overlap=80),
        RegionDetail(name='BOTH', chrom='1', start=1000, end=2000, variant_type='cnv_both', required_overlap=50),
    ]
    # DEL 90% of region, 100% of CNV -> both
    assert matching_regions('chr1', 1000, 1900, 'DEL', regions) == {'LOSS', 'BOTH'}
    # DEL 100% of region, only 50% of CNV -> both, CNV size doesn't matter
    assert matching_regions('chr1', 1000, 3000, 'DEL', regions) == {'LOSS', 'BOTH'}
    # DEL 70% of region, hanging off the start -> BOTH only, too small for LOSS's 80%
    assert matching_regions('chr1', 300, 1700, 'DEL', regions) == {'BOTH'}
    # DEL 20% of region, fully inside it -> both, containment only needs regulatory_region_percentage (15)
    assert matching_regions('chr1', 1400, 1600, 'DEL', regions) == {'LOSS', 'BOTH'}
    # DEL 10% of region, fully inside it -> none, below regulatory_region_percentage
    assert not matching_regions('chr1', 1400, 1500, 'DEL', regions)
    # DUP 90% of region -> BOTH only, wrong type for LOSS
    assert matching_regions('1', 1100, 2000, 'DUP', regions) == {'BOTH'}
    # elsewhere, wrong contig, or no END -> none
    assert not matching_regions('chr1', 5000, 6000, 'DEL', regions)
    assert not matching_regions('chr2', 1000, 1900, 'DEL', regions)
    assert not matching_regions('chr1', 1000, None, 'DEL', regions)


SV_VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=FAIL,Description="Failed">
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=PREDICTED_LOF,Number=.,Type=String,Description="LoF genes">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=N_HET,Number=1,Type=Integer,Description="Het count">
##INFO=<ID=N_HOMALT,Number=1,Type=Integer,Description="HomAlt count">
##INFO=<ID=MALE_AF,Number=A,Type=Float,Description="Male AF">
##INFO=<ID=FEMALE_AF,Number=A,Type=Float,Description="Female AF">
##INFO=<ID=gnomad_v4.1_sv_AF,Number=1,Type=Float,Description="gnomAD SV AF">
##INFO=<ID=gnomad_v4.1_sv_SVID,Number=1,Type=String,Description="gnomAD SV ID">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale
"""

COMMON_INFO = 'SVTYPE=DEL;SVLEN=5000;END=6000;AC=1;AF=0.01;AN=100;N_HET=1;N_HOMALT=0;MALE_AF=0.01;FEMALE_AF=0.02'


class TestSvMain:
    """End-to-end tests for the streaming SV labelling process."""

    def run_sv(self, tmp_path, vcf_rows, header=SV_VCF_HEADER):
        vcf_path = tmp_path / 'input.vcf'
        vcf_path.write_text(header + ''.join(f'{row}\n' for row in vcf_rows))

        panel = PanelApp(
            genes={'ENSG1': GeneDetail(symbol='GENE1', chrom='1')},
            regions={
                'REGION': RegionDetail(
                    name='REGION', chrom='1', start=10000, end=15000, variant_type='cnv_loss', required_overlap=80
                ),
            },
        )
        pa_path = tmp_path / 'panelapp.json'
        pa_path.write_text(panel.model_dump_json())

        ped_path = tmp_path / 'pedigree.ped'
        ped_path.write_text('family_1\tmale\t0\t0\t1\t2\n')

        out_path = str(tmp_path / 'labelled.vcf.bgz')
        sv_main(
            vcf_path=str(vcf_path),
            panelapp_path=str(pa_path),
            pedigree=str(ped_path),
            vcf_out=out_path,
        )
        return list(VCF(out_path))

    def test_rare_lof_sv_is_labelled(self, tmp_path):
        variants = self.run_sv(
            tmp_path,
            [
                (
                    'chr1\t1000\tsv_1\tN\t<DEL>\t50\tPASS\t'
                    f'PREDICTED_LOF=GENE1;gnomad_v4.1_sv_AF=0.001;gnomad_v4.1_sv_SVID=gnomad-del-1;{COMMON_INFO}\tGT\t0/1'
                ),
            ],
        )
        assert len(variants) == 1
        variant = variants[0]
        assert variant.INFO['categorybooleansv1'] == 1
        assert variant.INFO['gene_id'] == 'ENSG1'
        assert variant.INFO['gnomad_sv_AF'] == pytest.approx(0.001)
        assert variant.INFO['gnomad_sv_ID'] == 'gnomad-del-1'
        assert variant.INFO['svtype'] == 'DEL'
        assert variant.INFO['lof'] == 'GENE1'
        # raw fields are not carried through
        assert variant.INFO.get('PREDICTED_LOF') is None
        assert variant.INFO.get('SVTYPE') is None

    def test_frequency_fields_and_svlen_optional(self, tmp_path):
        """a caller with no AC/AF/AN/N_HET/N_HOMALT/per-sex AF/SVLEN is still labelled, with a warning"""
        optional = tuple(
            f'##INFO=<ID={field},' for field in ('SVLEN', 'AC', 'AF', 'AN', 'N_HET', 'N_HOMALT', 'MALE_AF', 'FEMALE_AF')
        )
        header = ''.join(line for line in SV_VCF_HEADER.splitlines(keepends=True) if not line.startswith(optional))
        messages: list[str] = []
        handler = logger.add(messages.append, format='{message}')
        try:
            variants = self.run_sv(
                tmp_path,
                ['chr1\t1000\tsv_1\tN\t<DEL>\t50\tPASS\tPREDICTED_LOF=GENE1;SVTYPE=DEL;END=6000\tGT\t0/1'],
                header=header,
            )
        finally:
            logger.remove(handler)

        (variant,) = variants
        assert variant.INFO['categorybooleansv1'] == 1
        assert variant.INFO.get('svlen') is None
        assert variant.INFO.get('AC') is None
        assert any('SV VCF has no AC, AF, AN, N_HET, N_HOMALT' in message for message in messages)
        assert any('callset AF filter is skipped' in message for message in messages)

    @pytest.mark.parametrize(('af', 'kept'), [(0.01, 1), (0.9, 0)])
    def test_overall_af_fallback(self, tmp_path, af, kept):
        """without per-sex frequencies, the callset AF filter falls back to overall AF"""
        header = ''.join(
            line
            for line in SV_VCF_HEADER.splitlines(keepends=True)
            if not line.startswith(('##INFO=<ID=MALE_AF,', '##INFO=<ID=FEMALE_AF,'))
        )
        info = COMMON_INFO.replace(';MALE_AF=0.01;FEMALE_AF=0.02', '').replace('AF=0.01', f'AF={af}')
        variants = self.run_sv(
            tmp_path,
            [f'chr1\t1000\tsv_1\tN\t<DEL>\t50\tPASS\tPREDICTED_LOF=GENE1;gnomad_v4.1_sv_AF=0.001;{info}\tGT\t0/1'],
            header=header,
        )
        assert len(variants) == kept

    def test_region_match_gets_its_own_row(self, tmp_path):
        """a LoF DEL matching a region is written once per gene and once per region, each with only its category"""
        # the other tests' CNVs are at 1000-6000, clear of the region at 10000-15000
        region_info = COMMON_INFO.replace('END=6000', 'END=15000')
        variants = self.run_sv(
            tmp_path,
            [
                (
                    'chr1\t10000\tsv_1\tN\t<DEL>\t50\tPASS\t'
                    f'PREDICTED_LOF=GENE1;gnomad_v4.1_sv_AF=0.001;gnomad_v4.1_sv_SVID=gnomad-del-1;{region_info}\tGT\t0/1'
                ),
                # region match without any LoF consequence is retained, but a common one is not
                f'chr1\t10000\tsv_2\tN\t<DEL>\t50\tPASS\tgnomad_v4.1_sv_AF=0.001;{region_info}\tGT\t0/1',
                f'chr1\t10000\tsv_3\tN\t<DEL>\t50\tPASS\tgnomad_v4.1_sv_AF=0.5;{region_info}\tGT\t0/1',
            ],
        )
        assert [
            (
                v.ID,
                v.INFO['gene_id'],
                v.INFO['categorybooleansv1'],
                v.INFO['categorybooleanregulatorycnv'],
            )
            for v in variants
        ] == [('sv_1', 'ENSG1', 1, 0), ('sv_1', 'REGION', 0, 1), ('sv_2', 'REGION', 0, 1)]

    def test_common_in_gnomad_is_dropped(self, tmp_path):
        variants = self.run_sv(
            tmp_path,
            [
                (
                    'chr1\t1000\tsv_1\tN\t<DEL>\t50\tPASS\t'
                    f'PREDICTED_LOF=GENE1;gnomad_v4.1_sv_AF=0.5;gnomad_v4.1_sv_SVID=gnomad-del-1;{COMMON_INFO}\tGT\t0/1'
                ),
            ],
        )
        assert not variants

    def test_non_lof_and_filtered_svs_are_dropped(self, tmp_path):
        variants = self.run_sv(
            tmp_path,
            [
                # no PREDICTED_LOF
                (
                    'chr1\t1000\tsv_1\tN\t<DEL>\t50\tPASS\t'
                    f'gnomad_v4.1_sv_AF=0.001;gnomad_v4.1_sv_SVID=gnomad-del-1;{COMMON_INFO}\tGT\t0/1'
                ),
                # filter failed
                (
                    'chr1\t2000\tsv_2\tN\t<DEL>\t50\tFAIL\t'
                    f'PREDICTED_LOF=GENE1;gnomad_v4.1_sv_AF=0.001;gnomad_v4.1_sv_SVID=gnomad-del-2;{COMMON_INFO}\tGT\t0/1'
                ),
                # LoF in a non-green gene
                (
                    'chr1\t3000\tsv_3\tN\t<DEL>\t50\tPASS\t'
                    f'PREDICTED_LOF=UNLISTED;gnomad_v4.1_sv_AF=0.001;gnomad_v4.1_sv_SVID=gnomad-del-3;{COMMON_INFO}\tGT\t0/1'
                ),
            ],
        )
        assert not variants

    def test_no_shared_samples_raises(self, tmp_path):
        vcf_path = tmp_path / 'input.vcf'
        vcf_path.write_text(SV_VCF_HEADER)
        panel = PanelApp(genes={'ENSG1': GeneDetail(symbol='GENE1', chrom='1')})
        pa_path = tmp_path / 'panelapp.json'
        pa_path.write_text(panel.model_dump_json())
        ped_path = tmp_path / 'pedigree.ped'
        ped_path.write_text('family_1\tsomeone_else\t0\t0\t1\t2\n')

        with pytest.raises(ValueError, match='No samples shared'):
            sv_main(
                vcf_path=str(vcf_path),
                panelapp_path=str(pa_path),
                pedigree=str(ped_path),
                vcf_out=str(tmp_path / 'labelled.vcf.bgz'),
            )
