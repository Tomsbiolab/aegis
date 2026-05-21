"""
Tests for aegis.equivalence — utility functions.
"""

import math
import pytest

from aegis.utils.evalue import parse_evalue, round_evalue


# ============================================================
# round_evalue
# ============================================================

class TestRoundEvalue:
    def test_small_evalue(self):
        result = round_evalue("1e-50")
        assert result == "1.00e-50"

    def test_zero(self):
        result = round_evalue("0")
        assert result == "0.00e+00"

    def test_large_value(self):
        result = round_evalue("10")
        assert result == "1.00e+01"

    def test_decimal(self):
        result = round_evalue("0.005")
        assert result == "5.00e-03"


# ============================================================
# parse_evalue
# ============================================================

class TestParseEvalue:
    def test_normal_evalue(self):
        result = parse_evalue("1e-50")
        assert result == pytest.approx(1e-50)

    def test_greater_than_prefix(self):
        result = parse_evalue(">10")
        assert result == pytest.approx(10.0)

    def test_na_string(self):
        result = parse_evalue("NA")
        assert math.isnan(result)

    def test_nan_string(self):
        result = parse_evalue("nan")
        assert math.isnan(result)

    def test_empty_string(self):
        result = parse_evalue("")
        assert math.isnan(result)

    def test_whitespace_stripped(self):
        result = parse_evalue("  1e-10  ")
        assert result == pytest.approx(1e-10)

