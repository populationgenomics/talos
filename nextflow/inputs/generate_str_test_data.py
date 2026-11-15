"""
script used to create an artificial joint-called STR VCF for testing the STR path through ValidateMOI

STR VCFs aren't produced by the annotation workflow - they come from STRipy JSON reports, merged by
talos.scripts.stripy_json_to_vcf. Rather than simulating STRipy JSON, this builds the per-sample locus dicts that
stripy_json_to_vcf.load_sample returns, and hands them to the real write_multisample_vcf, so the header, GT encoding
(1 = pathogenic-range allele, '.' = absent 2nd allele) and DISEASE_DETAILS format all match production output.

Loci are real PanelApp STR entities (panelapp_2026-09.json), all Green on the Mendeliome (137) apart from ATXN2.
Repeat coordinates are approximate GRCh38 positions inside each gene span, not exact STRipy catalogue coordinates -
nothing downstream checks them.

Every locus except DMPK sits on a contig present in nextflow/inputs/joint.vcf.bgz. ValidateMOI only iterates the
contigs declared in the small-variant VCF header, so STRs elsewhere are never read.

Expected behaviour against the premerged test cohort (proband: affected male, mother, father: unaffected), using
nextflow/inputs/config.toml. That sets pheno_match_strs = false - with the default (true), no STR survives, as none of
the proband's phenotype-matched panels carry one. "Gene MOI" is flagged where ValidateMOI currently applies the
PanelApp *gene* MOI rather than the STR MOI:

| Locus  | STR MOI     | Genotypes                         | Expected                                               |
|--------|-------------|-----------------------------------|--------------------------------------------------------|
| HOXD13 | Mono        | proband het, parents normal       | reported, de novo dominant                             |
| XYLT1  | Biallelic   | proband hom, parents het          | reported, homozygous                                   |
| FMR1   | X-linked    | proband hemi, mother het          | reported, hemizygous                                   |
| ARX_1  | X-linked    | proband hemi, parents normal      | reported, unless ARX_1 is in ValidateMOI.noisy_strs    |
| SOX3   | X-linked    | proband hemi, STR_FILTER=LowDepth | reported - STR_FILTER is not read at present           |
| VWA1   | Biallelic   | proband + mother het              | not reported, single het                               |
| GLS    | Biallelic   | proband + father het              | not reported (gene MOI Mono_And_Biallelic, father het) |
| AFF2   | X-linked    | mother het only                   | not reported, only an unaffected female carrier        |
| ABCD3  | Mono        | proband no call, parents normal   | row dropped, no carriers                               |
| BCLAF3 | X-linked    | all normal                        | row dropped, no carriers                               |
| ATXN2  | not on 137  | proband het, parents normal       | dropped unless an HPO-matched panel carries it         |
| STARD7 | Mono        | proband het, parents normal       | should report; dropped today - STR-only, no gene entry |
| DMPK   | Mono        | proband het, parents normal       | never read - chr19 absent from small VCF header        |

Run from the repository root:
    uv run python nextflow/inputs/generate_str_test_data.py
"""

import json
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

# write_multisample_vcf looks up an optional locus allow-list via config_retrieve, which needs a config path
os.environ.setdefault('TALOS_CONFIG', str(Path(__file__).parent / 'config.toml'))

from talos.scripts.stripy_json_to_vcf import write_multisample_vcf

OUTPUT = Path(__file__).parent / 'joint_str.vcf'

# matches the sample column order in nextflow/inputs/joint.vcf.bgz; proband and father male, mother female
SAMPLES = ['proband', 'mother', 'father']
MALES = {'proband', 'father'}


@dataclass
class StrLocus:
    locus_id: str
    ensg: str
    chrom: str
    start: int
    end: int
    motif: str
    disease: str
    inheritance: str
    normal_max: int
    pathogenic: int
    # sample -> repeat counts, one entry per allele. A sample absent here had no genotyping result for this locus
    calls: dict[str, tuple[int, ...]] = field(default_factory=dict)
    # sample -> STRipy Filter value, anything absent is PASS
    filters: dict[str, str] = field(default_factory=dict)

    def allele_range(self, repeats: int) -> str:
        if repeats >= self.pathogenic:
            return 'Pathogenic'
        if repeats > self.normal_max:
            return 'Intermediate'
        return 'Normal'

    @property
    def disease_details(self) -> str:
        inter = f'min{self.normal_max + 1}max{self.pathogenic - 1}' if self.pathogenic - self.normal_max > 1 else '.'
        return f'{self.disease}__{self.inheritance}__min1max{self.normal_max}__{inter}__{self.pathogenic}'

    def sample_record(self, sample: str) -> dict | None:
        """Build the per-sample locus dict in the shape stripy_json_to_vcf.load_sample produces."""
        if (repeats := self.calls.get(sample)) is None:
            return None
        a1 = repeats[0]
        a2 = repeats[1] if len(repeats) > 1 else None
        return {
            'chrom': self.chrom,
            'pos': self.start,
            'end': self.end,
            'id': self.locus_id,
            'motif': self.motif,
            'period': len(self.motif),
            'a1_rep': float(a1),
            'a2_rep': float(a2) if a2 is not None else None,
            'a1_ci': (max(a1 - 2, 0), a1 + 2),
            'a2_ci': (max(a2 - 2, 0), a2 + 2) if a2 is not None else (None, None),
            'a1_range': self.allele_range(a1),
            'a2_range': self.allele_range(a2) if a2 is not None else None,
            'a1_out': int(a1 >= self.pathogenic),
            'a2_out': int(a2 >= self.pathogenic) if a2 is not None else None,
            'a1_z': 8.0 if a1 >= self.pathogenic else 0.1,
            'a2_z': (8.0 if a2 >= self.pathogenic else 0.1) if a2 is not None else None,
            'coverage': 12 if sample in self.filters else 40,
            'filter': self.filters.get(sample, 'PASS'),
            'diseases': self.disease,
            'disease_details': self.disease_details,
        }


LOCI = [
    StrLocus(
        'HOXD13', 'ENSG00000128714', 'chr2', 176093059, 176093103, 'GCG', 'SPD1', 'AD', 15, 22,
        calls={'proband': (15, 29), 'mother': (15, 15), 'father': (15, 15)},
    ),
    StrLocus(
        'XYLT1', 'ENSG00000103489', 'chr16', 17470908, 17470922, 'GCC', 'DBQD2', 'AR', 20, 120,
        calls={'proband': (350, 380), 'mother': (9, 350), 'father': (9, 380)},
    ),
    StrLocus(
        'FMR1', 'ENSG00000102081', 'chrX', 147912051, 147912110, 'CGG', 'FXS', 'XLD', 44, 200,
        calls={'proband': (450,), 'mother': (30, 120), 'father': (30,)},
    ),
    StrLocus(
        'ARX_1', 'ENSG00000004848', 'chrX', 25013650, 25013697, 'GCN', 'EIEE1', 'XLR', 16, 23,
        calls={'proband': (27,), 'mother': (16, 16), 'father': (16,)},
    ),
    StrLocus(
        'SOX3', 'ENSG00000134595', 'chrX', 140503360, 140503405, 'GCN', 'PHPX', 'XLR', 15, 22,
        calls={'proband': (26,), 'mother': (15, 15), 'father': (15,)},
        filters={'proband': 'LowDepth'},
    ),
    StrLocus(
        'VWA1', 'ENSG00000179403', 'chr1', 1435798, 1435817, 'GCGCGGAGCG', 'HMNMYO', 'AR', 2, 3,
        calls={'proband': (2, 3), 'mother': (2, 3), 'father': (2, 2)},
    ),
    StrLocus(
        'GLS', 'ENSG00000115419', 'chr2', 190880873, 190880920, 'GCA', 'GDPAG', 'AR', 16, 400,
        calls={'proband': (8, 680), 'mother': (8, 8), 'father': (8, 680)},
    ),
    StrLocus(
        'AFF2', 'ENSG00000155966', 'chrX', 148500638, 148500682, 'GCC', 'FRAXE', 'XLR', 44, 200,
        calls={'proband': (15,), 'mother': (15, 250), 'father': (15,)},
    ),
    StrLocus(
        'ABCD3', 'ENSG00000117528', 'chr1', 94418422, 94418442, 'GCC', 'OPDM', 'AD', 50, 118,
        calls={'mother': (8, 8), 'father': (8, 8)},
    ),
    StrLocus(
        'BCLAF3', 'ENSG00000173681', 'chrX', 19990900, 19990930, 'CCG', 'FRAXG', 'XLR', 57, 117,
        calls={'proband': (12,), 'mother': (12, 14), 'father': (14,)},
    ),
    StrLocus(
        'ATXN2', 'ENSG00000204842', 'chr12', 111598951, 111599019, 'CAG', 'SCA2', 'AD', 31, 35,
        calls={'proband': (22, 41), 'mother': (22, 22), 'father': (22, 23)},
    ),
    StrLocus(
        'STARD7', 'ENSG00000084090', 'chr2', 96197067, 96197121, 'ATTTC', 'FAME2', 'AD', 20, 150,
        calls={'proband': (11, 400), 'mother': (11, 11), 'father': (11, 11)},
    ),
    StrLocus(
        'DMPK', 'ENSG00000104936', 'chr19', 45770205, 45770264, 'CTG', 'DM1', 'AD', 34, 50,
        calls={'proband': (5, 800), 'mother': (5, 13), 'father': (13, 13)},
    ),
]  # fmt: skip


def main():
    """write the VCF via the production STRipy writer, then bgzip and tabix it in place"""
    per_sample = [
        (sample, {locus.locus_id: record for locus in LOCI if (record := locus.sample_record(sample)) is not None})
        for sample in SAMPLES
    ]

    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as mapping:
        # write_multisample_vcf maps the symbol, i.e. the locus ID up to the first '_', not the full locus ID
        json.dump({locus.locus_id.split('_')[0]: locus.ensg for locus in LOCI}, mapping)

    try:
        write_multisample_vcf(samples=per_sample, out_path=str(OUTPUT), gene_map=mapping.name)
    finally:
        Path(mapping.name).unlink()

    bgzipped = OUTPUT.with_suffix('.vcf.bgz')
    with bgzipped.open('wb') as handle:
        subprocess.run(['bgzip', '-c', str(OUTPUT)], stdout=handle, check=True)  # noqa: S603, S607
    subprocess.run(['tabix', '-f', '-p', 'vcf', str(bgzipped)], check=True)  # noqa: S603, S607
    OUTPUT.unlink()
    print(f'wrote {len(LOCI)} STR loci to {bgzipped}')


if __name__ == '__main__':
    main()
