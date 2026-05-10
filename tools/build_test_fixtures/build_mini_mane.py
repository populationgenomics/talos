"""Subset the prod MANE summary to rows for the test-set genes.

Output: test/fixtures/large_files/MANE.GRCh38.v1.5.summary.txt.gz - same header and
columns as the production file, just with most rows removed. The resulting fixture
is consumed by the prep workflow's ParseManeIntoJson process and (indirectly) by
the talos workflow via the pre-baked mane.json.

The MANE column 'Ensembl_Gene' carries a version suffix (e.g. ENSG00000010704.14)
that we strip when matching - mirrors parse_mane_into_json.py's behaviour.
"""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO_ROOT = Path(__file__).resolve().parents[2]
PROD_MANE = REPO_ROOT / 'large_files' / 'MANE.GRCh38.v1.5.summary.txt.gz'
OVERLAP_TSV = Path(__file__).resolve().parent / 'locus_gff_overlap.tsv'
OUTPUT = REPO_ROOT / 'test' / 'fixtures' / 'large_files' / 'MANE.GRCh38.v1.5.summary.txt.gz'


def read_target_gene_ids() -> set[str]:
    gene_ids: set[str] = set()
    lines = OVERLAP_TSV.read_text().splitlines()
    header = lines[0].split('\t')
    gid = header.index('gene_id')
    for row in lines[1:]:
        gene_ids.add(row.split('\t')[gid])
    return gene_ids


def main() -> int:
    if not PROD_MANE.exists():
        print(f'[ERROR] prod MANE summary not found at {PROD_MANE}', file=sys.stderr)
        return 1
    if not OVERLAP_TSV.exists():
        print(f'[ERROR] run survey_gff.py first', file=sys.stderr)
        return 1

    target_genes = read_target_gene_ids()
    kept = 0

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(PROD_MANE, 'rt') as src, gzip.open(OUTPUT, 'wt') as dst:
        header = src.readline()
        dst.write(header)
        col_idx = header.rstrip('\n').split('\t').index('Ensembl_Gene')
        for line in src:
            cols = line.rstrip('\n').split('\t')
            if len(cols) <= col_idx:
                continue
            ensg = cols[col_idx].split('.')[0]
            if ensg in target_genes:
                dst.write(line)
                kept += 1

    print(f'[OK] wrote {kept} MANE rows to {OUTPUT}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
