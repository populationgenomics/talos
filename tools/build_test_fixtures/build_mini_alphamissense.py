"""Subset the prod AlphaMissense TSV to rows touching the test loci.

Output: test/fixtures/large_files/AlphaMissense_hg38.tsv.gz

The annotation workflow doesn't read this directly - it reads the pre-encoded
echtvar zip in processed_annotations. But the prep workflow's EncodeAlphaMissense
process consumes the TSV, and we want to keep the prep workflow runnable against
the same fixture set. We retain rows on test contigs within +/-3 kb of any test
locus to give EncodeAlphaMissense some content to work with.
"""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from loci import ALL_LOCI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
PROD_TSV = REPO_ROOT / 'large_files' / 'AlphaMissense_hg38.tsv.gz'
OUTPUT = REPO_ROOT / 'test' / 'fixtures' / 'large_files' / 'AlphaMissense_hg38.tsv.gz'

WINDOW_BP = 3000


def windows_by_contig() -> dict[str, list[tuple[int, int]]]:
    """Build per-contig list of (start, end) windows the AlphaMissense TSV is sampled from."""
    out: dict[str, list[tuple[int, int]]] = {}
    for locus in ALL_LOCI:
        # Skip mtDNA - AlphaMissense doesn't cover MT.
        if locus.chrom in {'chrM', 'M', 'MT'}:
            continue
        out.setdefault(locus.chrom, []).append((locus.pos - WINDOW_BP, locus.pos + WINDOW_BP))
    return out


def main() -> int:
    if not PROD_TSV.exists():
        print(f'[ERROR] missing {PROD_TSV}', file=sys.stderr)
        return 1

    windows = windows_by_contig()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with gzip.open(PROD_TSV, 'rt') as src, gzip.open(OUTPUT, 'wt') as dst:
        for line in src:
            if line.startswith('#'):
                dst.write(line)
                continue
            cols = line.split('\t', 3)
            if len(cols) < 4:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            if chrom not in windows:
                continue
            if any(start <= pos <= end for start, end in windows[chrom]):
                dst.write(line)
                kept += 1

    print(f'[OK] wrote {kept} AlphaMissense rows to {OUTPUT}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
