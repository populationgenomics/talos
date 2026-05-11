"""Survey the prod Ensembl GFF3 for what overlaps each test locus.

Run locally against `large_files/Homo_sapiens.GRCh38.116.MTtoM.gff3.gz`. Writes
`tools/build_test_fixtures/locus_gff_overlap.tsv` with one row per
(locus, overlapping gene). Used for sanity-checking the mini GFF and as the input
to the GFF-subset step (the gene IDs are extracted from this file).

The survey only inspects feature_type='gene' so the output stays small. The full
mini GFF (gene + mRNA + exon + CDS) is generated separately in build_mini_gff.py.
"""

import gzip
import sys
from collections import defaultdict
from pathlib import Path

# Allow running as `python tools/build_test_fixtures/survey_gff.py`
sys.path.insert(0, str(Path(__file__).resolve().parent))
from loci import ALL_LOCI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
PROD_GFF = REPO_ROOT / 'large_files' / 'Homo_sapiens.GRCh38.116.MTtoM.gff3.gz'
OUTPUT_TSV = Path(__file__).resolve().parent / 'locus_gff_overlap.tsv'


def parse_attrs(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in attr_field.rstrip(';').split(';'):
        if '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        out[k.strip()] = v.strip()
    return out


def main() -> int:
    if not PROD_GFF.exists():
        print(f'[ERROR] prod GFF not found at {PROD_GFF}', file=sys.stderr)
        return 1

    # Group test loci by the GFF contig name (Ensembl uses unprefixed chrom names).
    # The prod GFF was already chrMT->chrM renamed; nuclear contigs in Ensembl are
    # unprefixed (1, 2, ...) while ours are chr1, chr2... Handle both.
    contigs_of_interest: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for locus in ALL_LOCI:
        # Match either "chr1" or "1" form.
        ensembl_name = locus.chrom.removeprefix('chr')
        contigs_of_interest[locus.chrom].append((locus.pos, locus.note))
        contigs_of_interest[ensembl_name].append((locus.pos, locus.note))
        # MT vs M
        if ensembl_name == 'M':
            contigs_of_interest['MT'].append((locus.pos, locus.note))

    rows: list[dict[str, str]] = []
    with gzip.open(PROD_GFF, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            seqid, _src, ftype, start_s, end_s, _score, strand, _phase, attrs = fields
            # Catch gene, ncRNA_gene, pseudogene, miRNA_gene etc. - any record whose
            # ID is gene:... at the top of the gene tree.
            if not attrs.startswith('ID=gene:'):
                continue
            if seqid not in contigs_of_interest:
                continue
            start, end = int(start_s), int(end_s)
            for pos, note in contigs_of_interest[seqid]:
                if start <= pos <= end:
                    a = parse_attrs(attrs)
                    rows.append({
                        'contig': seqid,
                        'pos': str(pos),
                        'gene_id': a.get('gene_id', ''),
                        'gene_name': a.get('Name', a.get('gene_name', '')),
                        'biotype': a.get('biotype', ''),
                        'gene_start': str(start),
                        'gene_end': str(end),
                        'strand': strand,
                        'note': note,
                    })

    OUTPUT_TSV.write_text(
        '\t'.join(['contig', 'pos', 'gene_id', 'gene_name', 'biotype', 'gene_start', 'gene_end', 'strand', 'note'])
        + '\n'
        + '\n'.join('\t'.join(r[k] for k in ['contig', 'pos', 'gene_id', 'gene_name', 'biotype', 'gene_start', 'gene_end', 'strand', 'note']) for r in rows)
        + '\n',
    )
    print(f'[OK] wrote {len(rows)} overlap rows to {OUTPUT_TSV}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
