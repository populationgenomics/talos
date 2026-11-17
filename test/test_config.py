"""
config_retrieve when no config can be loaded
"""

import pytest

from talos import config
from talos.config import ConfigError, config_retrieve


@pytest.fixture
def no_config(monkeypatch):
    """TALOS_CONFIG unset, and no cached config"""
    monkeypatch.delenv(config.ENV_VAR, raising=False)
    monkeypatch.setattr(config, '_config', None)


def test_no_config_no_default_raises(no_config):  # noqa: ARG001
    with pytest.raises(ConfigError):
        config_retrieve(['a', 'b'])


def test_no_config_default_none_returns_none(no_config):  # noqa: ARG001
    assert config_retrieve(['a', 'b'], None) is None


def test_no_config_default_returned(no_config):  # noqa: ARG001
    assert config_retrieve(['a', 'b'], 'fallback') == 'fallback'
