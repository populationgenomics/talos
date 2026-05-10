"""Subset the prod Ensembl GFF3 to the genes hit by the test VCF + the full MT contig.

Strategy:
    1. Read locus_gff_overlap.tsv (produced by survey_gff.py) for the gene IDs of interest.
    2. Walk the prod GFF once. Keep:
         - every header line
         - every record on the MT contig (whole chrM is small)
         - every record whose ID/Parent/derives_from chain leads back to a kept gene_id
    3. Apply the prod chrMT->chrM rename on output.
    4. Gzip the result to test/fixtures/large_files/Homo_sapiens.GRCh38.116.MTtoM.gff3.gz.

Walking the GFF only once requires that records are emitted parent-before-child.
Ensembl's GFF3 is sorted gene -> mRNA/snRNA/etc -> exon/CDS so this holds.
"""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO_ROOT = Path(__file__).resolve().parents[2]
PROD_GFF = REPO_ROOT / 'large_files' / 'Homo_sapiens.GRCh38.116.MTtoM.gff3.gz'
OVERLAP_TSV = Path(__file__).resolve().parent / 'locus_gff_overlap.tsv'
OUTPUT_GFF = REPO_ROOT / 'test' / 'fixtures' / 'large_files' / 'Homo_sapiens.GRCh38.116.MTtoM.gff3.gz'


def read_target_gene_ids() -> set[str]:
    """Return the set of ENSG IDs we want to retain (excluding chrM, which is kept whole)."""
    gene_ids: set[str] = set()
    lines = OVERLAP_TSV.read_text().splitlines()
    header = lines[0].split('\t')
    gid = header.index('gene_id')
    contig = header.index('contig')
    for row in lines[1:]:
        cols = row.split('\t')
        if cols[contig] == 'M':
            continue
        gene_ids.add(cols[gid])
    return gene_ids


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
    if not OVERLAP_TSV.exists():
        print(f'[ERROR] run survey_gff.py first to produce {OVERLAP_TSV}', file=sys.stderr)
        return 1

    target_genes = read_target_gene_ids()
    print(f'[INFO] retaining {len(target_genes)} target nuclear genes + entire chrM contig')

    # ID strings we will retain; seeded with `gene:<gene_id>` for each target.
    kept_ids: set[str] = {f'gene:{g}' for g in target_genes}
    kept_lines: list[str] = []
    headers_kept: list[str] = []

    with gzip.open(PROD_GFF, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                # Keep ##gff-version, ###, etc.; skip per-contig sequence-region directives
                # for contigs we won't emit, to keep the file small.
                if line.startswith('##sequence-region'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] not in {'M', 'MT'}:
                        # Only keep sequence-regions for nuclear contigs we'll touch.
                        # We know which contigs from the overlap TSV.
                        continue
                headers_kept.append(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            seqid = fields[0]
            attrs_field = fields[8]

            # Whole chrM kept.
            if seqid in {'M', 'MT'}:
                if seqid == 'MT':
                    fields[0] = 'M'
                kept_lines.append('\t'.join(fields))
                continue

            attrs = parse_attrs(attrs_field)
            this_id = attrs.get('ID')
            parent = attrs.get('Parent')
            derives = attrs.get('Derives_from')

            # Decide whether this record belongs to a target gene by checking ancestry.
            ancestors_to_check: list[str] = []
            if parent:
                ancestors_to_check.extend(parent.split(','))
            if derives:
                ancestors_to_check.extend(derives.split(','))

            if this_id in kept_ids or any(a in kept_ids for a in ancestors_to_check):
                kept_lines.append(line.rstrip('\n'))
                if this_id:
                    kept_ids.add(this_id)

    # Add minimal sequence-region headers for the contigs we emit (Ensembl contig names).
    nuclear_contigs = sorted({line.split('\t')[0] for line in kept_lines if line.split('\t')[0] != 'M'})
    print(f'[INFO] emitting nuclear contigs: {nuclear_contigs}; total feature lines: {len(kept_lines)}')

    OUTPUT_GFF.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUTPUT_GFF, 'wt') as out:
        for h in headers_kept:
            out.write(h)
        for line in kept_lines:
            out.write(line + '\n')

    print(f'[OK] wrote mini GFF to {OUTPUT_GFF}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
