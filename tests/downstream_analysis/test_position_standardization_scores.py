"""
Tests for score transformation semantics in position_standardization.py.
"""

import importlib.util
from pathlib import Path
import sys


_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "position_standardization.py"
_spec = importlib.util.spec_from_file_location("position_standardization", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)


def test_g_stat_score_is_passthrough():
    assert _mod._score_metric_from_raw("g_stat", 42.5) == 42.5


def test_neglog10_fdr_score_is_passthrough():
    assert _mod._score_metric_from_raw("neglog10_fdr", 97.2) == 97.2


def test_fdr_score_is_transformed_to_neglog10():
    value = _mod._score_metric_from_raw("fdr", 1e-5)
    assert value > 4.9

