"""
Tests for aegis.equivalence — utility functions.
"""

import math
import pytest
from aegis.equivalence import parse_evalue, round_evalue, clean_annotation_tag


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


# ============================================================
# clean_annotation_tag
# ============================================================

class TestCleanAnnotationTag:
    def test_no_change_needed(self):
        assert clean_annotation_tag("myannot_v1") == "myannot_v1"

    def test_remove_lifton_prefix(self):
        assert clean_annotation_tag("Lifton_myannot") == "myannot"

    def test_remove_liftoff_prefix(self):
        assert clean_annotation_tag("Liftoff_myannot") == "myannot"

    def test_remove_lifton_word(self):
        assert clean_annotation_tag("myLiftonannot") == "myannot"

    def test_split_on_from(self):
        assert clean_annotation_tag("annot_from_genome1") == "annot"

    def test_split_on_on(self):
        assert clean_annotation_tag("annot_on_genome1") == "annot"

    def test_split_on_to(self):
        assert clean_annotation_tag("annot_to_genome1") == "annot"

    def test_combined_liftoff_and_from(self):
        result = clean_annotation_tag("Liftoff_myannot_from_genome")
        assert result == "myannot"
