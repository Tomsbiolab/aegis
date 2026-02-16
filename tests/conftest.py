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
# Minimal GFF3 fixture (1 gene, 1 mRNA, 2 exons, 1 CDS segment)
# ---------------------------------------------------------------------------

MINIMAL_GFF3 = """\
##gff-version 3
chr1\taegis\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=TestGene
chr1\taegis\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=TestTranscript
chr1\taegis\texon\t1000\t2000\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\taegis\texon\t3000\t5000\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\taegis\tCDS\t1200\t2000\t.\t+\t0\tID=CDS1;Parent=mRNA1
chr1\taegis\tCDS\t3000\t4500\t.\t+\t0\tID=CDS1;Parent=mRNA1
chr1\taegis\tfive_prime_UTR\t1000\t1199\t.\t+\t.\tID=utr5_1;Parent=mRNA1
chr1\taegis\tthree_prime_UTR\t4501\t5000\t.\t+\t.\tID=utr3_1;Parent=mRNA1
"""

# Multi-gene GFF3: 2 genes on different chromosomes, plus a minus-strand gene
MULTI_GENE_GFF3 = """\
##gff-version 3
chr1\taegis\tgene\t1000\t5000\t.\t+\t.\tID=geneA;Name=GeneA
chr1\taegis\tmRNA\t1000\t5000\t.\t+\t.\tID=mRNA_A;Parent=geneA
chr1\taegis\texon\t1000\t2000\t.\t+\t.\tID=exonA1;Parent=mRNA_A
chr1\taegis\texon\t3000\t5000\t.\t+\t.\tID=exonA2;Parent=mRNA_A
chr1\taegis\tCDS\t1200\t2000\t.\t+\t0\tID=CDS_A;Parent=mRNA_A
chr2\taegis\tgene\t500\t3000\t.\t-\t.\tID=geneB;Name=GeneB
chr2\taegis\tmRNA\t500\t3000\t.\t-\t.\tID=mRNA_B;Parent=geneB
chr2\taegis\texon\t500\t1500\t.\t-\t.\tID=exonB1;Parent=mRNA_B
chr2\taegis\texon\t2000\t3000\t.\t-\t.\tID=exonB2;Parent=mRNA_B
chr2\taegis\tCDS\t600\t1500\t.\t-\t0\tID=CDS_B;Parent=mRNA_B
"""


@pytest.fixture
def sample_gff3_file(tmp_path):
    """Write minimal GFF3 to a temp file and return its path."""
    gff = tmp_path / "sample.gff3"
    gff.write_text(MINIMAL_GFF3)
    return str(gff)


@pytest.fixture
def multi_gene_gff3_file(tmp_path):
    """Write multi-gene GFF3 to a temp file and return its path."""
    gff = tmp_path / "multi.gff3"
    gff.write_text(MULTI_GENE_GFF3)
    return str(gff)


# ---------------------------------------------------------------------------
# Minimal FASTA fixture (2 chromosomes)
# ---------------------------------------------------------------------------

MINIMAL_FASTA = """\
>chr1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGAT
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
"""


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Write minimal FASTA to a temp file and return its path."""
    fa = tmp_path / "sample.fasta"
    fa.write_text(MINIMAL_FASTA)
    return str(fa)


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
