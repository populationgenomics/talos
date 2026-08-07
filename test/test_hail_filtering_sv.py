"""
Unit tests for the filtering helpers in talos.run_hail_filtering_sv.

These exercise the pure Hail transforms, which need no external tools (unlike the SVAnnotate/SVAFotate steps).
They cover the field rearrangement, the two blanket frequency filters, and the hemizygous-call fix.
"""

import json

import hail as hl

from talos.run_hail_filtering_sv import (
    filter_matrix_by_ac,
    filter_matrix_by_af,
    fix_hemi_calls,
    read_and_filter_mane_json,
    rearrange_annotations,
)


def _annotated_sv_mt(use_af_male_spelling: bool = False) -> hl.MatrixTable:
    """
    build a single-row SV MatrixTable carrying every INFO field rearrange_annotations reads unconditionally

    Args:
        use_af_male_spelling (bool): if True, provide AF_MALE/AF_FEMALE instead of MALE_AF/FEMALE_AF
    """
    mt = hl.utils.range_matrix_table(n_rows=1, n_cols=1)
    info = {
        'PREDICTED_LOF': hl.literal(['GENE1']),
        'SVTYPE': hl.literal('DEL'),
        'SVLEN': hl.int32(5000),
        'END': hl.int32(6000),
        'AC': hl.array([hl.int32(1)]),
        'AF': hl.array([hl.float64(0.01)]),
        'AN': hl.int32(100),
        'N_HET': hl.int32(1),
        'N_HOMALT': hl.int32(0),
        'gnomad_v4.1_sv_AF': hl.float64(0.001),
        'gnomad_v4.1_sv_SVID': hl.literal('gnomad-del-1'),
    }
    if use_af_male_spelling:
        info['AF_MALE'] = hl.array([hl.float64(0.03)])
        info['AF_FEMALE'] = hl.array([hl.float64(0.04)])
    else:
        info['MALE_AF'] = hl.array([hl.float64(0.03)])
        info['FEMALE_AF'] = hl.array([hl.float64(0.04)])
    return mt.annotate_rows(info=hl.struct(**info))


def test_read_and_filter_mane_json(tmp_path):
    """the MANE JSON is read into a symbol -> ENSG lookup"""
    mane = tmp_path / 'mane.json'
    mane.write_text(
        json.dumps(
            {
                'ENST1': {'symbol': 'GENE1', 'ensg': 'ENSG1'},
                'ENST2': {'symbol': 'GENE2', 'ensg': 'ENSG2'},
            },
        ),
    )
    mapping = read_and_filter_mane_json(str(mane))
    assert hl.eval(mapping.get('GENE1')) == 'ENSG1'
    assert hl.eval(mapping.get('GENE2')) == 'ENSG2'


def test_rearrange_surfaces_gnomad_and_defaults():
    """the population fields are surfaced under the names Talos reads, and absent SV fields are defaulted"""
    mt = _annotated_sv_mt()
    out = rearrange_annotations(mt, hl.literal({'GENE1': 'ENSG1'}))
    row = out.rows().collect()[0]

    assert row.info.gnomad_sv_AF == 0.001
    assert row.info.gnomad_sv_ID == 'gnomad-del-1'
    assert row.info.svtype == 'DEL'

    # ALGORITHMS was absent, so it defaults to ['gCNV']; STATUS/CHR2/END2 are inserted as null
    assert row.info.algorithms == ['gCNV']
    assert row.info.status is None
    assert row.info.chr2 is None
    assert row.info.end2 is None


def test_rearrange_maps_symbols_to_ensg():
    """PREDICTED_LOF symbols are mapped to gene IDs, and unknown symbols pass through unchanged"""
    mt = _annotated_sv_mt()
    out = rearrange_annotations(mt, hl.literal({'GENE1': 'ENSG1'}))
    row = out.rows().collect()[0]
    assert set(row.info.lof) == {'GENE1'}
    assert set(row.info.lof_ensg) == {'ENSG1'}
    assert set(row.info.gene_id) == {'ENSG1'}


def test_rearrange_symbol_without_mapping_is_kept():
    """a symbol missing from the MANE map is retained as itself, not dropped"""
    mt = _annotated_sv_mt()
    out = rearrange_annotations(mt, hl.literal({'OTHER': 'ENSG9'}))
    row = out.rows().collect()[0]
    assert set(row.info.lof_ensg) == {'GENE1'}


def test_rearrange_reads_male_af():
    """MALE_AF/FEMALE_AF are read through as the array-typed male_af/female_af"""
    out = rearrange_annotations(_annotated_sv_mt(), hl.literal({'GENE1': 'ENSG1'}))
    row = out.rows().collect()[0]
    assert row.info.male_af == [0.03]
    assert row.info.female_af == [0.04]


def test_rearrange_accepts_af_male_alternative():
    """the AF_MALE/AF_FEMALE spelling is accepted as an alternative to MALE_AF/FEMALE_AF"""
    out = rearrange_annotations(_annotated_sv_mt(use_af_male_spelling=True), hl.literal({'GENE1': 'ENSG1'}))
    row = out.rows().collect()[0]
    assert row.info.male_af == [0.03]
    assert row.info.female_af == [0.04]


def test_filter_matrix_by_af_removes_common_but_keeps_missing():
    """AF above threshold is removed; a null AF is treated as 0 and survives"""
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=1)
    mt = mt.annotate_rows(
        info=hl.struct(
            gnomad_sv_AF=hl.if_else(
                mt.row_idx == 0,
                hl.float64(0.5),
                hl.if_else(mt.row_idx == 1, hl.missing(hl.tfloat64), hl.float64(0.001)),
            ),
        ),
    )
    out = filter_matrix_by_af(mt, af_threshold=0.03)
    # row 0 (0.5) dropped; row 1 (missing -> 0) kept; row 2 (0.001) kept
    assert out.count_rows() == 2
    assert set(out.info.gnomad_sv_AF.collect()) == {None, 0.001}


def test_filter_matrix_by_ac_needs_both_sexes_below_threshold():
    """a variant survives only when both male and female callset AFs are at or below the threshold"""
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=1)
    mt = mt.annotate_rows(
        info=hl.struct(
            male_af=hl.if_else(mt.row_idx == 1, hl.array([hl.float64(0.9)]), hl.array([hl.float64(0.01)])),
            female_af=hl.if_else(mt.row_idx == 2, hl.array([hl.float64(0.9)]), hl.array([hl.float64(0.02)])),
        ),
    )
    out = filter_matrix_by_ac(mt, ac_threshold=0.03)
    # row 0 kept; row 1 dropped (male high); row 2 dropped (female high)
    assert out.count_rows() == 1


def test_fix_hemi_calls_diploidises_haploid_genotypes():
    """haploid calls are recast biallelic (alt -> 1/1, ref -> 0/0); diploid calls are untouched"""
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=1)
    mt = mt.annotate_entries(
        GT=hl.if_else(
            mt.row_idx == 0,
            hl.call(1),  # haploid alt
            hl.if_else(mt.row_idx == 1, hl.call(0), hl.call(0, 1)),  # haploid ref, then diploid het
        ),
    )
    out = fix_hemi_calls(mt)
    out = out.annotate_entries(dip=out.GT.is_diploid(), n_alt=out.GT.n_alt_alleles())
    by_idx = {entry.row_idx: entry for entry in out.entries().collect()}

    assert by_idx[0].dip and by_idx[0].n_alt == 2  # 1 -> 1/1
    assert by_idx[1].dip and by_idx[1].n_alt == 0  # 0 -> 0/0
    assert by_idx[2].dip and by_idx[2].n_alt == 1  # 0/1 unchanged
