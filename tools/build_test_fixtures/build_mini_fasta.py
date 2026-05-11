"""Build a mini reference FASTA for CI fixtures.

Strategy:
    - Keep only the contigs touched by the test VCF (chr1, chr2, chr6, chr11, chr12,
      chr16, chrX, chrM). All other contigs are dropped.
    - Each kept contig retains its real GRCh38 length, but only positions inside
      keep-windows hold real bases. Everything else is masked to N.
    - Keep-windows:
        * the entire chrM contig
        * each test gene's full extent (so bcftools csq has all exons+introns)
        * each test locus +/- FASTA_FLANK_BP (covers splice/csq context for loci
          that fall outside any retained gene)
    - Output is bgzip'd to keep size down (long N-runs compress to almost nothing)
      and a .fai index is written via `samtools faidx`.

Output: test/fixtures/large_files/ref.fa.gz + .fai (bgzipped, indexed).
The output is bgzip-compressed: 1.2 GB of mostly-N sequence collapses to ~4 MB,
which is shippable. The CI test profile sets `params.ref_genome = ".../ref.fa.gz"`.
samtools / bcftools / pysam / Hail all accept bgzipped FASTA + .fai + .gzi.
"""

import subprocess
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from loci import ALL_LOCI, FASTA_FLANK_BP  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
PROD_FASTA = REPO_ROOT / 'large_files' / 'ref.fa'
OVERLAP_TSV = Path(__file__).resolve().parent / 'locus_gff_overlap.tsv'
OUTPUT_FASTA_GZ = REPO_ROOT / 'test' / 'fixtures' / 'large_files' / 'ref.fa.gz'

# Wrap output FASTA at 60 chars to match common conventions and keep diffs clean.
LINE_WIDTH = 60


def keep_contigs() -> set[str]:
    return {locus.chrom for locus in ALL_LOCI}


def build_keep_windows() -> dict[str, list[tuple[int, int]]]:
    """Return contig -> sorted, merged list of [start, end] inclusive 1-based windows."""
    windows: dict[str, list[tuple[int, int]]] = {c: [] for c in keep_contigs()}

    # Test loci: +/- flank around each variant.
    for locus in ALL_LOCI:
        windows[locus.chrom].append((max(1, locus.pos - FASTA_FLANK_BP), locus.pos + FASTA_FLANK_BP))

    # Whole chrM: trivially one window covering the whole contig.
    # We don't know the contig length yet; use a sentinel that the merge step caps later.
    windows['chrM'].append((1, 10**9))

    # Each retained gene's full extent.
    overlap_lines = OVERLAP_TSV.read_text().splitlines()[1:]
    for row in overlap_lines:
        cols = row.split('\t')
        contig = cols[0]
        gene_start = int(cols[5])
        gene_end = int(cols[6])
        # Survey TSV records contig as Ensembl-style names ('1' not 'chr1', 'M' not 'chrM').
        chrom = f'chr{contig}' if contig != 'M' else 'chrM'
        if chrom in windows:
            windows[chrom].append((gene_start, gene_end))

    # Sort and merge.
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, items in windows.items():
        if not items:
            continue
        items.sort()
        out: list[tuple[int, int]] = [items[0]]
        for s, e in items[1:]:
            if s <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def main() -> int:
    if not PROD_FASTA.exists():
        print(f'[ERROR] missing {PROD_FASTA}', file=sys.stderr)
        return 1
    fai = PROD_FASTA.with_suffix(PROD_FASTA.suffix + '.fai')
    if not fai.exists():
        print(f'[ERROR] missing index at {fai} - run `samtools faidx {PROD_FASTA}` first', file=sys.stderr)
        return 1

    contig_lengths: dict[str, int] = {}
    with open(fai) as handle:
        for line in handle:
            name, length, *_ = line.rstrip('\n').split('\t')
            contig_lengths[name] = int(length)

    windows = build_keep_windows()

    # Cap each window to the actual contig length.
    for chrom, items in windows.items():
        clen = contig_lengths.get(chrom)
        if clen is None:
            continue
        windows[chrom] = [(s, min(e, clen)) for s, e in items if s <= clen]

    print(f'[INFO] keeping contigs: {sorted(windows.keys())}')

    OUTPUT_FASTA_GZ.parent.mkdir(parents=True, exist_ok=True)
    # Write to a temporary uncompressed file first, then bgzip it. samtools faidx
    # produces both .fai and .gzi indices for a bgzip'd FASTA.
    tmp_fasta = OUTPUT_FASTA_GZ.with_suffix('')  # strip .gz -> ref.fa
    fasta = pysam.FastaFile(str(PROD_FASTA))
    try:
        with open(tmp_fasta, 'w') as out:
            for chrom in sorted(windows.keys(), key=lambda c: list(contig_lengths.keys()).index(c)):
                clen = contig_lengths[chrom]
                seq = bytearray(b'N' * clen)
                for start, end in windows[chrom]:
                    real = fasta.fetch(chrom, start - 1, end).upper().encode('ascii')
                    seq[start - 1:start - 1 + len(real)] = real
                kept_bp = sum(e - s + 1 for s, e in windows[chrom])
                print(f'[INFO]   {chrom}: {clen} bp total, {kept_bp} bp real ({100 * kept_bp / clen:.4f}%)')
                out.write(f'>{chrom}\n')
                for i in range(0, clen, LINE_WIDTH):
                    out.write(seq[i:i + LINE_WIDTH].decode('ascii'))
                    out.write('\n')
    finally:
        fasta.close()

    # bgzip in place: produces ref.fa.gz, removes ref.fa.
    if OUTPUT_FASTA_GZ.exists():
        OUTPUT_FASTA_GZ.unlink()
    subprocess.run(['bgzip', '--force', str(tmp_fasta)], check=True)
    # Index: produces .fai and .gzi
    subprocess.run(['samtools', 'faidx', str(OUTPUT_FASTA_GZ)], check=True)
    print(f'[OK] wrote {OUTPUT_FASTA_GZ} ({OUTPUT_FASTA_GZ.stat().st_size:,} bytes)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
