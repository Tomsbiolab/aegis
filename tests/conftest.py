"""
Shared Pytest fixtures for the Aegis test suite.
"""

import os
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA_DIR = REPO_ROOT / "tests" / "test_data"


@pytest.fixture
def test_data_dir():
    """Return the path to the real test_data/ directory."""
    return TEST_DATA_DIR


# ---------------------------------------------------------------------------
# GFF3 fixtures — loaded from test_data/
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_gff3_file():
    """Return path to the minimal GFF3 test file."""
    return str(TEST_DATA_DIR / "minimal.gff3")


@pytest.fixture
def multi_gene_gff3_file():
    """Return path to the multi-gene GFF3 test file."""
    return str(TEST_DATA_DIR / "multi_gene.gff3")


@pytest.fixture
def rich_gff3_file():
    """Return path to the rich GFF3 test file (3 genes, 2 chromosomes, alt transcripts)."""
    return str(TEST_DATA_DIR / "rich.gff3")


# ---------------------------------------------------------------------------
# FASTA fixture — loaded from test_data/
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_fasta_file():
    """Return path to the minimal FASTA test file."""
    return str(TEST_DATA_DIR / "minimal.fasta")


# ---------------------------------------------------------------------------
# Ready-made class instances
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_feature():
    """A simple Feature instance for quick tests."""
    from aegis.feature import Feature
    return Feature(
        feature_id="feat001",
        ch="chr1",
        source="aegis",
        feature="gene",
        strand="+",
        start=100,
        end=500,
        score=".",
        phase=".",
        attributes="ID=feat001;Name=TestFeature;Alias=TF1"
    )


@pytest.fixture
def sample_gene():
    """A simple Gene instance for quick tests."""
    from aegis.gene import Gene
    return Gene(
        pseudogene=False,
        transposable=False,
        feature_id="gene001",
        ch="chr1",
        source="aegis",
        feature="gene",
        strand="+",
        start=100,
        end=2000,
        score=".",
        phase=".",
        attributes="ID=gene001;Name=TestGene"
    )
