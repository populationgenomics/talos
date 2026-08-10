"""
Tests for the cyvcf2 streaming helpers, and the streaming mito labelling process built on them.

These encode the same behavioural spec as the Hail implementations they replace.
"""

import gzip
from typing import ClassVar

import pytest
from cyvcf2 import VCF

from talos.models import PanelApp, PanelDetail
from talos.reformat_and_label_mito_vcf import main as mito_main
from talos.vcf_streaming import (
    PATHOGENIC,
    consequences_to_csq_string,
    normalise_chrom,
    parse_bcsq_entries,
    read_clinvar_decisions_tsv,
    split_csq_header,
)

BCSQ_FIELDS = ['consequence', 'gene', 'transcript', 'biotype', 'strand', 'amino_acid_change', 'dna_change']

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

MITO_VCF_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chrM,length=16569>\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FILTER=<ID=FAIL,Description="Failed">\n'
    f'{BCSQ_HEADER_LINE}\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale\n'
)

CLINVAR_TSV_HEADER = '\t'.join(
    ['contig', 'position', 'reference', 'alternate', 'clinical_significance', 'gold_stars', 'allele_id'],
)


def make_vcf(tmp_path, rows: list[str], name: str = 'input.vcf'):
    vcf_path = tmp_path / name
    vcf_path.write_text(MITO_VCF_HEADER + ''.join(f'{row}\n' for row in rows))
    return str(vcf_path)


def make_clinvar_tsv(tmp_path, rows: list[str], gzipped: bool = False):
    content = CLINVAR_TSV_HEADER + '\n' + ''.join(f'{row}\n' for row in rows)
    if gzipped:
        tsv_path = tmp_path / 'clinvar.tsv.gz'
        with gzip.open(tsv_path, 'wt') as handle:
            handle.write(content)
    else:
        tsv_path = tmp_path / 'clinvar.tsv'
        tsv_path.write_text(content)
    return str(tsv_path)


def make_panelapp(tmp_path, genes: dict[str, tuple[str, str]]):
    """genes: {ensg: (symbol, chrom)}"""
    panel = PanelApp(
        genes={ensg: PanelDetail(symbol=symbol, chrom=chrom) for ensg, (symbol, chrom) in genes.items()},
    )
    pa_path = tmp_path / 'panelapp.json'
    pa_path.write_text(panel.model_dump_json())
    return str(pa_path)


def fake_config_retrieve(key, default=None):
    """Stand-in config: BCSQ-era csq_string, defaults for everything else."""
    if key == ['RunHailFiltering', 'csq_string']:
        return CSQ_STRING_CONFIG
    return default


def make_pedigree(tmp_path, sample: str = 'male'):
    ped_path = tmp_path / 'pedigree.ped'
    ped_path.write_text(f'family_1\t{sample}\t0\t0\t1\t2\n')
    return str(ped_path)


def test_normalise_chrom():
    assert normalise_chrom('chrM') == 'M'
    assert normalise_chrom('MT') == 'M'
    assert normalise_chrom('M') == 'M'
    assert normalise_chrom('chr1') == '1'
    assert normalise_chrom('X') == 'X'


def test_parse_bcsq_full_entry():
    csqs = parse_bcsq_entries('missense|MT-ND1|ENST0001|protein_coding|+|62M>62V|100A>G', BCSQ_FIELDS)
    assert len(csqs) == 1
    csq = csqs[0]
    # strand is dropped, codon is derived
    assert 'strand' not in csq
    assert csq == {
        'consequence': 'missense',
        'gene': 'MT-ND1',
        'transcript': 'ENST0001',
        'biotype': 'protein_coding',
        'amino_acid_change': '62M>62V',
        'dna_change': '100A>G',
        'codon': 62,
    }


def test_parse_bcsq_truncated_entries():
    """BCFtools csq truncates the string when trailing fields are absent - pad them back."""
    for truncated in (
        'synonymous|GENE|TX',
        'synonymous|GENE|TX|protein_coding',
        'synonymous|GENE|TX|protein_coding|+',
        'synonymous|GENE|TX|protein_coding|+|',
    ):
        csq = parse_bcsq_entries(truncated, BCSQ_FIELDS)[0]
        assert csq['amino_acid_change'] == ''
        assert csq['dna_change'] == ''
        assert csq['codon'] is None


def test_parse_bcsq_multiple_entries():
    csqs = parse_bcsq_entries(
        'missense|GENE1|TX1|protein_coding|+|10K>10E|100A>G,synonymous|GENE2|TX2|protein_coding|-|55L|200C>T',
        BCSQ_FIELDS,
    )
    assert [csq['gene'] for csq in csqs] == ['GENE1', 'GENE2']
    assert [csq['codon'] for csq in csqs] == [10, 55]


def test_parse_bcsq_non_numeric_aa_change():
    csq = parse_bcsq_entries('missense|GENE|TX|protein_coding|+|?-62V|100A>G', BCSQ_FIELDS)[0]
    assert csq['codon'] is None


def test_consequences_to_csq_string(monkeypatch):
    monkeypatch.setattr('talos.vcf_streaming.config_retrieve', lambda _key: CSQ_STRING_CONFIG)
    csq_string = consequences_to_csq_string(
        [
            {
                'consequence': 'missense',
                'gene': 'MT-ND1',
                'gene_id': 'ENSG0001',
                'transcript': 'TX1',
                'biotype': 'protein_coding',
                'amino_acid_change': '62M>62V',
                'dna_change': '100A>G',
                'codon': 62,
            },
        ],
    )
    # fields absent from the consequence dict (mane_id, mane, ensp, am_*) render as empty strings
    assert csq_string == 'missense|ENSG0001|MT-ND1|TX1|||protein_coding|100A>G|62M>62V|62|||'


def test_consequences_to_csq_string_multiple(monkeypatch):
    monkeypatch.setattr('talos.vcf_streaming.config_retrieve', lambda _key: ['consequence', 'gene'])
    csq_string = consequences_to_csq_string(
        [{'consequence': 'missense', 'gene': 'A'}, {'consequence': 'synonymous', 'gene': 'B'}],
    )
    assert csq_string == 'missense|A,synonymous|B'


@pytest.mark.parametrize('gzipped', [False, True])
def test_read_clinvar_decisions_tsv(tmp_path, gzipped):
    tsv = make_clinvar_tsv(
        tmp_path,
        [
            f'chrM\t100\tA\tG\t{PATHOGENIC}\t1\t1234',
            'chrM\t200\tC\tT\tbenign\t2\t5678',
            f'chr1\t300\tG\tA\t{PATHOGENIC}\t3\t9999',
        ],
        gzipped=gzipped,
    )
    decisions = read_clinvar_decisions_tsv(tsv)
    assert len(decisions) == 3
    assert decisions[('M', 100, 'A', 'G')] == {
        'clinical_significance': PATHOGENIC,
        'gold_stars': 1,
        'allele_id': 1234,
    }

    mito_only = read_clinvar_decisions_tsv(tsv, contigs={'M'})
    assert set(mito_only) == {('M', 100, 'A', 'G'), ('M', 200, 'C', 'T')}


def test_read_clinvar_decisions_tsv_bad_columns(tmp_path):
    tsv_path = tmp_path / 'bad.tsv'
    tsv_path.write_text('contig\tposition\tsurprise_column\nchrM\t100\tvalue\n')
    with pytest.raises(ValueError, match='surprise_column'):
        read_clinvar_decisions_tsv(str(tsv_path))


def test_split_csq_header(tmp_path):
    vcf = make_vcf(tmp_path, [])
    reader = VCF(vcf)
    assert split_csq_header(reader) == BCSQ_FIELDS
    reader.close()


class TestMitoMain:
    """End-to-end tests for the streaming mito labelling process."""

    GREEN_GENES: ClassVar[dict] = {
        'ENSG0001': ('MT-ND1', 'MT'),
        'ENSG0002': ('MT-CO1', 'MT'),
        'ENSG0099': ('NUCLEAR', '1'),
    }

    def run_mito(self, tmp_path, vcf_rows, clinvar_rows, genes=None, sample='male'):
        out_path = str(tmp_path / 'labelled.vcf.bgz')
        mito_main(
            vcf_path=make_vcf(tmp_path, vcf_rows),
            output_path=out_path,
            panelapp=make_panelapp(tmp_path, genes or self.GREEN_GENES),
            pedigree=make_pedigree(tmp_path, sample=sample),
            clinvar_path=make_clinvar_tsv(tmp_path, clinvar_rows),
        )
        return list(VCF(out_path))

    def test_pathogenic_green_variant_is_kept(self, tmp_path, monkeypatch):
        monkeypatch.setattr('talos.vcf_streaming.config_retrieve', fake_config_retrieve)
        variants = self.run_mito(
            tmp_path,
            vcf_rows=[
                'chrM\t100\t.\tA\tG\t50\tPASS\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|62M>62V|100A>G\tGT\t0/1'
            ],
            clinvar_rows=[f'chrM\t100\tA\tG\t{PATHOGENIC}\t2\t1234'],
        )
        assert len(variants) == 1
        variant = variants[0]
        assert variant.INFO['categorybooleanclinvarplp'] == 1
        assert variant.INFO['clinvar_talos'] == 1
        assert variant.INFO['clinvar_stars'] == 2
        assert variant.INFO['clinvar_allele'] == 1234
        assert variant.INFO['clinvar_significance'] == PATHOGENIC
        assert variant.INFO['gene_id'] == 'ENSG0001'
        assert variant.INFO['gnomad_AC'] == 0
        # BCSQ is replaced by the talos-formatted csq string, with the panelapp-derived gene_id
        assert variant.INFO.get('BCSQ') is None
        assert variant.INFO['csq'] == 'missense|ENSG0001|MT-ND1|TX1|||protein_coding|100A>G|62M>62V|62|||'

    def test_non_pathogenic_variants_are_dropped(self, tmp_path):
        variants = self.run_mito(
            tmp_path,
            vcf_rows=[
                # zero stars
                'chrM\t100\t.\tA\tG\t50\tPASS\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|62M>62V|100A>G\tGT\t0/1',
                # benign
                'chrM\t200\t.\tC\tT\t50\tPASS\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|63M>63V|200C>T\tGT\t0/1',
                # absent from clinvar entirely
                'chrM\t300\t.\tG\tA\t50\tPASS\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|64M>64V|300G>A\tGT\t0/1',
            ],
            clinvar_rows=[
                f'chrM\t100\tA\tG\t{PATHOGENIC}\t0\t1234',
                'chrM\t200\tC\tT\tbenign\t2\t5678',
            ],
        )
        assert not variants

    def test_filter_failed_variant_is_dropped(self, tmp_path):
        variants = self.run_mito(
            tmp_path,
            vcf_rows=[
                'chrM\t100\t.\tA\tG\t50\tFAIL\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|62M>62V|100A>G\tGT\t0/1'
            ],
            clinvar_rows=[f'chrM\t100\tA\tG\t{PATHOGENIC}\t2\t1234'],
        )
        assert not variants

    def test_non_green_gene_is_dropped(self, tmp_path):
        variants = self.run_mito(
            tmp_path,
            vcf_rows=[
                'chrM\t100\t.\tA\tG\t50\tPASS\tBCSQ=missense|MT-UNLISTED|TX1|protein_coding|+|62M>62V|100A>G\tGT\t0/1'
            ],
            clinvar_rows=[f'chrM\t100\tA\tG\t{PATHOGENIC}\t2\t1234'],
        )
        assert not variants

    def test_variant_in_two_green_genes_is_exploded(self, tmp_path, monkeypatch):
        monkeypatch.setattr('talos.vcf_streaming.config_retrieve', fake_config_retrieve)
        variants = self.run_mito(
            tmp_path,
            vcf_rows=[
                (
                    'chrM\t100\t.\tA\tG\t50\tPASS\t'
                    'BCSQ=missense|MT-ND1|TX1|protein_coding|+|62M>62V|100A>G,'
                    'missense|MT-CO1|TX2|protein_coding|+|10K>10E|100A>G\tGT\t0/1'
                ),
            ],
            clinvar_rows=[f'chrM\t100\tA\tG\t{PATHOGENIC}\t2\t1234'],
        )
        # one output row per green gene, each carrying the full csq set
        assert len(variants) == 2
        assert [v.INFO['gene_id'] for v in variants] == ['ENSG0001', 'ENSG0002']
        assert variants[0].INFO['csq'] == variants[1].INFO['csq']
        assert variants[0].INFO['csq'].count(',') == 1

    def test_no_mito_genes_writes_empty_vcf(self, tmp_path):
        with pytest.raises(SystemExit) as sys_exit:
            self.run_mito(
                tmp_path,
                vcf_rows=[
                    'chrM\t100\t.\tA\tG\t50\tPASS\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|62M>62V|100A>G\tGT\t0/1'
                ],
                clinvar_rows=[f'chrM\t100\tA\tG\t{PATHOGENIC}\t2\t1234'],
                genes={'ENSG0099': ('NUCLEAR', '1')},
            )
        assert sys_exit.value.code == 0
        assert not list(VCF(str(tmp_path / 'labelled.vcf.bgz')))

    def test_no_shared_samples_writes_empty_vcf(self, tmp_path):
        with pytest.raises(SystemExit) as sys_exit:
            self.run_mito(
                tmp_path,
                vcf_rows=[
                    'chrM\t100\t.\tA\tG\t50\tPASS\tBCSQ=missense|MT-ND1|TX1|protein_coding|+|62M>62V|100A>G\tGT\t0/1'
                ],
                clinvar_rows=[f'chrM\t100\tA\tG\t{PATHOGENIC}\t2\t1234'],
                sample='someone_else',
            )
        assert sys_exit.value.code == 0
        assert not list(VCF(str(tmp_path / 'labelled.vcf.bgz')))
