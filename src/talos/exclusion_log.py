"""
Optional 'super logging' for variant exclusions during MOI validation.

When enabled (via [ValidateMOI] super_logging), records one JSON line per
(variant, sample, applied_moi) exclusion to a file at [ValidateMOI]
super_logging_path. If the path ends with '.gz', output is gzip-compressed.

When disabled, every record() call is a single attribute check and returns
immediately, so the hot loops in moi_tests.py pay essentially nothing.

First-failure-wins: callers preserve their existing short-circuit semantics
and emit a single record per (variant, sample, MOI) iteration.
"""

import gzip
import json
from typing import Any

from loguru import logger

from talos.config import config_retrieve
from talos.models import VARIANT_MODELS

_logger: 'ExclusionLogger | None' = None


class ExclusionLogger:
    """Stream JSONL exclusion records. Disabled by default; configure to enable."""

    def __init__(self) -> None:
        self.enabled: bool = bool(config_retrieve(['ValidateMOI', 'super_logging'], False))
        self.path: str | None = config_retrieve(['ValidateMOI', 'super_logging_path'], None)
        self._handle: Any = None

        if self.enabled and not self.path:
            logger.warning('super_logging enabled but super_logging_path not set; disabling super logging')
            self.enabled = False

    def _open(self) -> None:
        if self._handle is not None or not self.enabled:
            return
        assert self.path is not None
        opener = gzip.open if self.path.endswith('.gz') else open
        self._handle = opener(self.path, 'wt')
        logger.info(f'Super logging enabled, writing exclusions to {self.path}')

    def record(
        self,
        *,
        variant: VARIANT_MODELS | None,
        gene: str | None,
        sample: str | None,
        applied_moi: str | None,
        stage: str,
        reason: str,
        details: dict[str, Any] | None = None,
    ) -> None:
        if not self.enabled:
            return
        self._open()
        payload = {
            'variant': variant.coordinates.string_format if variant is not None else None,
            'gene': gene,
            'sample': sample,
            'applied_moi': applied_moi,
            'stage': stage,
            'reason': reason,
            'details': details or {},
        }
        assert self._handle is not None
        self._handle.write(json.dumps(payload) + '\n')

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None


def get_exclusion_logger() -> ExclusionLogger:
    """Return the process-wide exclusion logger, initialising it on first use."""
    global _logger
    if _logger is None:
        _logger = ExclusionLogger()
    return _logger


def reset_exclusion_logger() -> None:
    """Close and reset the singleton. Intended for tests."""
    global _logger
    if _logger is not None:
        _logger.close()
    _logger = None
