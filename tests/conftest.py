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
    return str(TEST_DATA_DIR / "input/annotation/minimal.gff3")


@pytest.fixture
def multi_gene_gff3_file():
    """Return path to the multi-gene GFF3 test file."""
    return str(TEST_DATA_DIR / "input/annotation/multi_gene.gff3")


@pytest.fixture
def rich_gff3_file():
    """Return path to the rich GFF3 test file (3 genes, 2 chromosomes, alt transcripts)."""
    return str(TEST_DATA_DIR / "input/annotation/rich.gff3")


# ---------------------------------------------------------------------------
# Edge cases GFF3 files
# ---------------------------------------------------------------------------

@pytest.fixture
def exon_only_gff3_file():
    """Gene with only exons (no CDS/UTR) — non-coding mRNA."""
    return str(TEST_DATA_DIR / "input/annotation/exon_only.gff3")


@pytest.fixture
def cds_only_gff3_file():
    """Gene with only CDS segments (no explicit exons or UTRs)."""
    return str(TEST_DATA_DIR / "input/annotation/cds_only.gff3")


@pytest.fixture
def no_subfeatures_gff3_file():
    """Gene with mRNA but zero subfeatures."""
    return str(TEST_DATA_DIR / "input/annotation/no_subfeatures.gff3")


@pytest.fixture
def noncoding_transcripts_gff3_file():
    """Several non-coding transcript types (lnc_RNA, tRNA, rRNA, snoRNA)."""
    return str(TEST_DATA_DIR / "input/annotation/noncoding_transcripts.gff3")


@pytest.fixture
def pseudogene_gff3_file():
    """Pseudogene with exons pointing directly at the gene (no transcript)."""
    return str(TEST_DATA_DIR / "input/annotation/pseudogene.gff3")


@pytest.fixture
def overlapping_exons_gff3_file():
    """Gene with overlapping exons that need collapsing."""
    return str(TEST_DATA_DIR / "input/annotation/overlapping_exons.gff3")


@pytest.fixture
def multi_cds_ids_gff3_file():
    """Transcript with CDS segments having different IDs (polycistronic)."""
    return str(TEST_DATA_DIR / "input/annotation/multi_cds_ids.gff3")


@pytest.fixture
def transcript_no_parent_gff3_file():
    """Transcript without a Parent attribute (gene inferred)."""
    return str(TEST_DATA_DIR / "input/annotation/transcript_no_parent.gff3")


@pytest.fixture
def cds_no_parent_gff3_file():
    """CDS subfeatures with no Parent attribute (gene+transcript inferred)."""
    return str(TEST_DATA_DIR / "input/annotation/cds_no_parent.gff3")


@pytest.fixture
def multiple_isoforms_gff3_file():
    """Gene with 3 mRNA isoforms sharing exon regions but different CDS."""
    return str(TEST_DATA_DIR / "input/annotation/multiple_isoforms.gff3")


@pytest.fixture
def subfeature_parent_is_gene_gff3_file():
    """Exons and CDSs reference a gene ID as Parent (no mRNA/transcript line)."""
    return str(TEST_DATA_DIR / "input/annotation/subfeature_parent_is_gene.gff3")


@pytest.fixture
def geneID_attribute_as_parent_gff3_file():
    """GFF3 without gene features and with geneID rather than Parent."""
    return str(TEST_DATA_DIR / "input/annotation/geneID_attribute_as_parent.gff3")


# ---------------------------------------------------------------------------
# FASTA fixture — loaded from test_data/
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_fasta_file():
    """Return path to the minimal FASTA test file."""
    return str(TEST_DATA_DIR / "input/fasta/minimal.fasta")


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
        attributes={"Name": ["TestFeature"], "Alias": ["TF1"]}
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
        attributes={"Name": ["TestGene"]}
    )
