import pytest

from talos.static_values import get_evidence_date, get_granular_date


def test_evidence_date_defaults_to_today(monkeypatch):
    monkeypatch.delenv('TALOS_EVIDENCE_DATE', raising=False)
    assert get_evidence_date() == get_granular_date()


def test_evidence_date_can_be_historical(monkeypatch):
    monkeypatch.setenv('TALOS_EVIDENCE_DATE', '2025-10-07')
    assert get_evidence_date() == '2025-10-07'


def test_evidence_date_rejects_invalid_date(monkeypatch):
    monkeypatch.setenv('TALOS_EVIDENCE_DATE', '2025-02-30')
    with pytest.raises(ValueError, match='TALOS_EVIDENCE_DATE'):
        get_evidence_date()
