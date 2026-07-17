"""
Tests for the optional super-logging mode that records the reason each variant is
excluded by the MOI checks.
"""

import gzip
import json
from pathlib import Path

from mendelbrot.pedigree_parser import PedigreeParser

from talos import exclusion_log
from talos.models import Coordinates, SmallVariant
from talos.moi_tests import DominantAutosomal, RecessiveAutosomalCH

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
TEST_COORDS_PARTNER = Coordinates(chrom='1', pos=2, ref='A', alt='C')


def _read_records(path: Path) -> list[dict]:
    if str(path).endswith('.gz'):
        with gzip.open(path, 'rt') as fh:
            return [json.loads(line) for line in fh]
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_frequency_filter_records_reason(
    tmp_path,
    pedigree_path,
    moi_config,
):
    """A variant rejected for high gnomAD AF emits one frequency_filter record."""
    out = tmp_path / 'exclusions_filter.jsonl'

    _logger = exclusion_log.reset_exclusion_logger(str(out))

    info = {'gnomad_af': 0.5, 'ac': 0, 'af': 0.0, 'gene_id': 'GENEY'}
    variant = SmallVariant(
        info=info,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 100},
        alt_depths={'male': 50},
        transcript_consequences=[],
    )
    runner = DominantAutosomal(pedigree=PedigreeParser(pedigree_path), config=moi_config)
    assert runner.run(principal=variant) == []

    _logger.close()

    records = _read_records(out)
    # one frequency_filter record, no per-sample records
    assert len(records) == 1
    rec = records[0]
    assert rec['stage'] == 'frequency_filter'
    assert rec['reason'] == 'gnomad_af_too_high'
    assert rec['sample'] is None
    assert rec['applied_moi'] == 'Autosomal Dominant'
    assert rec['details']['value'] == 0.5


def test_first_failure_short_circuit(tmp_path, pedigree_path, moi_config):
    """A sample failing both depth and category checks records only the first failure."""
    out = tmp_path / 'exclusions_short.jsonl'

    _logger = exclusion_log.reset_exclusion_logger(str(out))

    # passes the population frequency, fails per-sample (no category, low depth, low alt)
    info = {'gnomad_af': 0.0, 'ac': 0, 'af': 0.0, 'gene_id': 'GENEZ'}
    variant = SmallVariant(
        info=info,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 1},
        alt_depths={'male': 0},
        transcript_consequences=[],
    )
    runner = DominantAutosomal(pedigree=PedigreeParser(pedigree_path), config=moi_config)
    runner.run(principal=variant)
    _logger.close()
    records = _read_records(out)

    # exactly one per-sample record for the male; the category check fires first
    sample_records = [r for r in records if r['sample'] == 'male']
    assert len(sample_records) == 1
    assert sample_records[0]['reason'] == 'sample_not_categorised'


def test_gz_path_writes_gzip(tmp_path, pedigree_path, moi_config):
    """Output path ending in .gz should produce a valid gzip stream."""
    out = tmp_path / 'exclusions_compressed.jsonl.gz'

    _logger = exclusion_log.reset_exclusion_logger(str(out))

    info = {'gnomad_af': 0.5, 'ac': 0, 'af': 0.0, 'gene_id': 'GENEY'}
    variant = SmallVariant(
        info=info,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 100},
        alt_depths={'male': 50},
        transcript_consequences=[],
    )
    DominantAutosomal(pedigree=PedigreeParser(pedigree_path), config=moi_config).run(principal=variant)
    _logger.close()

    records = _read_records(out)
    assert len(records) == 1
    assert records[0]['reason'] == 'gnomad_af_too_high'


def test_comp_het_partner_rejection_logged(tmp_path, pedigree_path, moi_config):
    """A comp-het partner with too-high gnomAD AF should be recorded as a comp_het rejection."""
    out = tmp_path / 'exclusions_rej.jsonl'

    _logger = exclusion_log.reset_exclusion_logger(str(out))

    boolean_categories = ['categorybooleansample']
    principal_info = {
        'gnomad_af': 0.0,
        'ac': 0,
        'af': 0.0,
        'gene_id': 'GENEY',
        'categorybooleansample': True,
    }
    principal = SmallVariant(
        info=principal_info,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 100},
        alt_depths={'male': 50},
        transcript_consequences=[],
        boolean_categories=boolean_categories,
    )
    partner_info = {
        'gnomad_af': 0.99,  # forces partner_too_common (after partner depth check passes)
        'ac': 0,
        'af': 0.0,
        'gene_id': 'GENEY',
        'categorybooleansample': True,
    }
    partner = SmallVariant(
        info=partner_info,
        het_samples={'male'},
        coordinates=TEST_COORDS_PARTNER,
        depths={'male': 100},
        alt_depths={'male': 50},
        transcript_consequences=[],
        boolean_categories=boolean_categories,
    )
    comp_het = {'male': {principal.coordinates.string_format: [partner]}}
    RecessiveAutosomalCH(pedigree=PedigreeParser(pedigree_path), config=moi_config).run(
        principal=principal, comp_het=comp_het
    )
    _logger.close()

    records = _read_records(out)
    # the partner-side too_common emits a frequency_filter record (variant-level) AND
    # the comp-het site emits a partner_too_common record (sample-level, stage=comp_het)
    comp_het_records = [r for r in records if r['stage'] == 'comp_het']
    assert any(r['reason'] == 'partner_too_common' for r in comp_het_records)
