"""Subset the prod hp.obo to test-pedigree HPO terms + their ancestor closure.

Output: test/fixtures/large_files/hp.obo

The HPO obo file is parsed by obonet into a networkx graph. unified_panelapp_parser.py
then walks dfs_successors() from each pedigree HPO term to the root and intersects with
HPOs that have associated panels. We retain enough of the graph that the dfs walk
reaches the root for every term in the pedigree.

genes_to_phenotype.txt subset: filter to gene symbols that appear in the test gene set.
"""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from loci import PEDIGREE_HPO_TERMS  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
PROD_OBO = REPO_ROOT / 'large_files' / 'hp.obo'
PROD_G2P = REPO_ROOT / 'large_files' / 'genes_to_phenotype.txt'

OUTPUT_OBO = REPO_ROOT / 'test' / 'fixtures' / 'large_files' / 'hp.obo'
OUTPUT_G2P = REPO_ROOT / 'test' / 'fixtures' / 'large_files' / 'genes_to_phenotype.txt'

OVERLAP_TSV = Path(__file__).resolve().parent / 'locus_gff_overlap.tsv'


def parse_obo_terms(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    """Return (header_lines, terms_by_id). Each term value is the raw lines (without [Term])."""
    header: list[str] = []
    terms: dict[str, list[str]] = {}
    with open(path, encoding='utf-8') as handle:
        in_term = False
        block: list[str] = []
        block_kind: str | None = None
        for line in handle:
            line_stripped = line.rstrip('\n')
            if line.startswith('['):
                if in_term and block_kind == 'Term' and block:
                    tid = next((l.split(' ', 1)[1].strip() for l in block if l.startswith('id:')), None)
                    if tid:
                        terms[tid] = block
                block = []
                in_term = True
                block_kind = line.strip()[1:-1]
                continue
            if in_term:
                block.append(line_stripped)
            else:
                header.append(line_stripped)
        if in_term and block_kind == 'Term' and block:
            tid = next((l.split(' ', 1)[1].strip() for l in block if l.startswith('id:')), None)
            if tid:
                terms[tid] = block
    return header, terms


def transitive_closure(seeds: list[str], terms: dict[str, list[str]]) -> set[str]:
    """Walk is_a edges from `seeds` until we reach the root (HP:0000001)."""
    keep: set[str] = set()
    queue = list(seeds)
    while queue:
        cur = queue.pop()
        if cur in keep or cur not in terms:
            continue
        keep.add(cur)
        for line in terms[cur]:
            if line.startswith('is_a:'):
                parent = line.split(':', 1)[1].split('!')[0].strip()
                queue.append(parent)
    return keep


def write_mini_obo() -> None:
    header, terms = parse_obo_terms(PROD_OBO)
    keep = transitive_closure(PEDIGREE_HPO_TERMS, terms)
    OUTPUT_OBO.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_OBO, 'w', encoding='utf-8') as handle:
        for line in header:
            handle.write(line + '\n')
        for tid in sorted(keep):
            handle.write('\n[Term]\n')
            for body_line in terms[tid]:
                handle.write(body_line + '\n')
    print(f'[OK] wrote {len(keep)} HPO terms to {OUTPUT_OBO}')


def collect_test_gene_symbols() -> set[str]:
    """Gene symbols hit by the test loci (extracted from the survey TSV)."""
    symbols: set[str] = set()
    lines = OVERLAP_TSV.read_text().splitlines()
    header = lines[0].split('\t')
    name_idx = header.index('gene_name')
    for row in lines[1:]:
        cols = row.split('\t')
        if cols[name_idx]:
            symbols.add(cols[name_idx])
    return symbols


def write_mini_g2p() -> None:
    target_symbols = collect_test_gene_symbols()
    kept = 0
    with _open_text(PROD_G2P) as src, open(OUTPUT_G2P, 'w', encoding='utf-8') as dst:
        header_line = src.readline()
        # The Jax g2p file's first column may or may not include a header; preserve verbatim.
        dst.write(header_line)
        for line in src:
            cols = line.split('\t')
            if len(cols) >= 2 and cols[1] in target_symbols:
                dst.write(line)
                kept += 1
    print(f'[OK] wrote {kept} g2p rows to {OUTPUT_G2P}')


def _open_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, encoding='utf-8')


def main() -> int:
    if not PROD_OBO.exists():
        print(f'[ERROR] missing {PROD_OBO}', file=sys.stderr)
        return 1
    if not PROD_G2P.exists():
        print(f'[ERROR] missing {PROD_G2P}', file=sys.stderr)
        return 1
    write_mini_obo()
    write_mini_g2p()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
