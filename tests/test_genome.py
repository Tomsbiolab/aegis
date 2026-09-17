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

    def test_get_stats(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        stats = g.get_stats()
        assert stats["total_size"] == 5800
        assert stats["chromosome_size"] == 5800
        assert stats["num_sequences"] == 2
        assert stats["num_chromosomes"] == 2
        assert stats["num_scaffolds"] == 0
        assert stats["n50"] == 5500
        assert stats["l50"] == 1
        assert stats["n90"] == 5500
        assert stats["l90"] == 1
        assert stats["gc_content"] == 50.0
        assert stats["gap_content"] == 0.0
        assert stats["longest_scaffold"] == ("chr1", 5500)
        assert stats["shortest_scaffold"] == ("chr2", 300)

    def test_stats_property(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        assert g.stats == g.get_stats()

    def test_get_sorted_features(self, sample_fasta_file):
        g = Genome("test_genome", sample_fasta_file, quiet=True)
        sorted_feats = g.get_sorted_features(sort_by="name")
        assert sorted_feats == ["chr1", "chr2"]

        sorted_by_size = g.get_sorted_features(sort_by="size")
        assert sorted_by_size == ["chr1", "chr2"]

    def test_assembly_stats_and_sorting_logic(self, tmp_path):
        """
        Verify underlying biological assembly metrics (N50, L50, GC content, Gap content)
        and natural chromosome ordering on a known synthetic genome.
        """
        fasta_file = tmp_path / "synthetic_genome.fasta"
        # chr1: 400 bp (all G/C)
        # chr10: 300 bp (all A/T)
        # chr2: 250 bp (balanced)
        # chrM: 150 bp (organelle)
        # scaffold_1: 100 bp (50 Ns)
        fasta_content = (
            ">chr1\n" + "GC" * 200 + "\n" +
            ">chr10\n" + "AT" * 150 + "\n" +
            ">chr2\n" + "ATGC" * 62 + "AT\n" +
            ">chrM\n" + "ATGC" * 37 + "AT\n" +
            ">scaffold_1\n" + "N" * 50 + "GC" * 25 + "\n"
        )
        fasta_file.write_text(fasta_content, encoding="utf-8")

        g = Genome("synthetic", str(fasta_file), quiet=True)
        stats = g.get_stats()

        # Assembly sizes
        assert stats["total_size"] == 1200
        # Nuclear chr (400 + 300 + 250 = 950) + organelle (150) = 1100
        assert stats["chromosome_size"] == 1100
        assert stats["nuclear_chromosome_size"] == 950
        assert stats["scaffold_size"] == 100
        assert stats["num_sequences"] == 5
        assert stats["num_chromosomes"] == 4
        assert stats["num_scaffolds"] == 1

        # N50 and L50 (sorted lengths: 400, 300, 250, 150, 100; half = 600; 400+300=700 >= 600)
        assert stats["n50"] == 300
        assert stats["l50"] == 2

        # N90 and L90 (90% = 1080; 400+300+250+150 = 1100 >= 1080)
        assert stats["n90"] == 150
        assert stats["l90"] == 4

        # auN (length-weighted mean): sum(L^2) / total = 345,000 / 1200 = 287.5 -> 288
        assert stats["auN"] == 288

        # With estimated genome size = 1500 (threshold = 750; 400+300+250=950 >= 750)
        stats_ng = g.get_stats(estimated_genome_size=1500)
        assert stats_ng["ng50"] == 250
        assert stats_ng["lg50"] == 3
        assert stats_ng["auNG"] == 230

        # With large estimated size = 3000 (threshold = 1500; assembly never reaches 50%)
        stats_large = g.get_stats(estimated_genome_size=3000)
        assert stats_large["ng50"] == 0
        assert stats_large["lg50"] == 0
        assert stats_large["auNG"] == 115

        # Gap content: 50 Ns out of 1200 bp = 4.17%
        assert stats["gap_content"] == round((50 / 1200) * 100, 2)

        # Natural sorting: nuclear (chr1, chr2, chr10) -> organelle (chrM) -> scaffold (scaffold_1)
        # Crucial: chr2 MUST precede chr10 (natural numeric sort, not lexicographical)
        sorted_all = g.get_sorted_features(chromosomes_only=False, sort_by="name")
        assert sorted_all == ["chr1", "chr2", "chr10", "chrM", "scaffold_1"]

        sorted_chrs = g.get_sorted_features(chromosomes_only=True, sort_by="name")
        assert sorted_chrs == ["chr1", "chr2", "chr10", "chrM"]

    def test_preserve_case_and_soft_masking(self, tmp_path):
        fasta_path = tmp_path / "soft_masked.fa"
        fasta_path.write_text(">chr1\nATGCatgcNNNN\n>chr2\natgcatgc\n")

        # Default: preserve_case=True
        g = Genome("soft", str(fasta_path), quiet=True)
        assert g.preserve_case is True
        assert g["chr1"].seq == "ATGCatgcNNNN"
        assert g["chr1"].upper_seq == "ATGCATGCNNNN"
        assert g["chr1"].soft_masked_bp == 4
        assert g["chr1"].soft_masked_fraction == round(4 / 12, 4)

        assert g["chr2"].seq == "atgcatgc"
        assert g["chr2"].soft_masked_bp == 8
        assert g["chr2"].soft_masked_fraction == 1.0

        stats = g.get_stats()
        assert stats["soft_masked_bp"] == 12
        assert stats["soft_masked_pct"] == round(12 / 20 * 100, 2)
        assert stats["gap_content"] == round(4 / 20 * 100, 2)

        # preserve_case=False forces uppercase
        g_upper = Genome("upper", str(fasta_path), quiet=True, preserve_case=False)
        assert g_upper.preserve_case is False
        assert g_upper["chr1"].seq == "ATGCATGCNNNN"
        assert g_upper["chr1"].soft_masked_bp == 0
        assert g_upper["chr1"].soft_masked_fraction == 0.0

        # Sequence checksums match regardless of casing
        assert g["chr1"].seq_hash == g_upper["chr1"].seq_hash
        assert g["chr2"].seq_hash == g_upper["chr2"].seq_hash



