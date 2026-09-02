"""
Tests for aegis.genome — Scaffold and Genome classes.
"""

import pytest

from aegis.genome import Scaffold, Genome


# ============================================================
# Scaffold
# ============================================================

class TestScaffold:
    def test_basic_init(self):
        s = Scaffold("chr1", "ATGCGATCGATCGATCGATCGATCGATCG")
        assert s.name == "chr1"
        assert s.size == 29
        assert s.original_name == "chr1"

    def test_custom_original_name(self):
        s = Scaffold("chr1_renamed", "ATGCGATCG", original_name="scaffold_001")
        assert s.original_name == "scaffold_001"
        assert s.name == "chr1_renamed"

    def test_chromosome_detection(self):
        """Names starting with 'ch' set chromosome=True."""
        s = Scaffold("chr1", "ATGCGATCG")
        assert s.chromosome is True

    def test_non_chromosome(self):
        """Names not starting with 'ch' are not chromosomes."""
        s = Scaffold("scaffold_99", "ATGCGATCG")
        assert s.chromosome is False

    def test_mitochondria_detection(self):
        s = Scaffold("chrM", "ATGCGATCG")
        assert s.mitochondria is True

    def test_chloroplast_detection(self):
        s = Scaffold("chrC", "ATGCGATCG")
        assert s.chloroplast is True

    def test_unknown_chromosome(self):
        """The 'chrUn' name is detected as unknown_chromosome."""
        s = Scaffold("chrUn", "ATGCGATCG")
        assert s.unknown_chromosome is True

    def test_update_with_new_name(self):
        s = Scaffold("scaffold_1", "ATGCGATCG")
        s.update(new_name="chr5")
        assert s.name == "chr5"
        assert s.chromosome is True

    def test_copy(self):
        s = Scaffold("chr1", "ATGCGATCG")
        s2 = s.copy()
        s2.name = "changed"
        assert s.name == "chr1"

    def test_dapfit_numbered_chromosome(self):
        s = Scaffold("chr1", "ATGCGATCG")
        assert s.dapfit is True

    def test_dapfit_non_numbered(self):
        s = Scaffold("chrM", "ATGCGATCG")
        assert s.dapfit is False


# ============================================================
# Genome
# ============================================================

class TestGenome:
    def test_load_from_fasta(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        assert g.name == "test_genome"
        assert len(g.scaffolds) == 2
        assert "chr1" in g.scaffolds
        assert "chr2" in g.scaffolds

    def test_scaffold_sizes(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        assert g.scaffolds["chr1"].size > 0
        assert g.scaffolds["chr2"].size > 0

    def test_update(self, sample_fasta_file):
        """Genome.update sets self.size (total genome size)."""
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        g.update()
        assert g.size > 0
        assert len(g.scaffolds) == 2

    def test_remove_features(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        g.remove_features({"chr2"})
        assert "chr2" not in g.scaffolds
        assert "chr1" in g.scaffolds

    def test_subset(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        g.subset(chosen_features={"chr1"})
        assert "chr1" in g.scaffolds

        # Test list input
        g_list = Genome("test_genome", sample_fasta_file, quiet=True)
        g_list.subset(chosen_features=["chr1"])
        assert "chr1" in g_list.scaffolds

        # Test tuple input
        g_tuple = Genome("test_genome", sample_fasta_file, quiet=True)
        g_tuple.subset(chosen_features=("chr1",))
        assert "chr1" in g_tuple.scaffolds

    def test_copy(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        g2 = g.copy()
        g2.name = "changed"
        assert g.name == "test_genome"
