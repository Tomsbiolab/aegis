"""
Tests for aegis.genome — Scaffold and Genome classes.
"""

import pytest

from aegis.genome import Scaffold, Genome, classify_compartment


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

    @pytest.mark.parametrize("name, description, expected", [
        ("chrM", "", "mitochondria"),
        ("chrMT", "", "mitochondria"),
        ("MT", "", "mitochondria"),
        ("Mt", "", "mitochondria"),
        ("Ch_M", "", "mitochondria"),
        ("mitochondrion", "", "mitochondria"),
        ("ChrC", "", "chloroplast"),
        ("chrPt", "", "chloroplast"),
        ("Pt", "", "chloroplast"),
        ("chloroplast", "", "chloroplast"),
        ("NC_012920.1", "NC_012920.1 Homo sapiens mitochondrion, complete genome", "mitochondria"),
        ("NC_000932.1", "NC_000932.1 Arabidopsis thaliana chloroplast, complete genome", "chloroplast"),
        ("chr1", "", None),
        ("chromosome_Sym", "", None),
        ("ChrSyntenic_c", "", None),
        ("scaffold_9", "", None),
        ("M", "", None),
        ("C", "", None),
    ])
    def test_organelle_classification(self, name, description, expected):
        s = Scaffold(name, "ATGCGATCG", description=description)
        assert Scaffold.organelle_type(name, description) == expected
        assert s.mitochondria is (expected == "mitochondria")
        assert s.chloroplast is (expected == "chloroplast")
        assert s.organelle is (expected is not None)
        if expected is not None:
            assert s.chromosome is True

    def test_organelle_flags_survive_renaming(self):
        """DAP renaming turns organelles into chr0001...; they must stay organelles."""
        s = Scaffold("chrM", "ATGCGATCG")
        s.update(new_name="chr0001")
        assert s.mitochondria is True
        assert s.organelle is True

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

    @pytest.mark.parametrize("name", ["chrUn", "ChrUn", "Chr00", "chr0", "chrUn_KI270302v1", "chrUn_NW_123.1", "chr00_scaffold_1", "Chr0_ctg5"])
    def test_unknown_chromosome_names(self, name):
        s = Scaffold(name, "ATGC")
        assert s.unknown_chromosome is True
        assert s.chromosome is False

    @pytest.mark.parametrize("name", ["chr1_random", "chr4_KI270790v1_alt", "chr1_KN196472v1_fix", "chr1_KI270706v1_random", "chr2_NW_0001.1"])
    def test_chromosome_scaffold_names(self, name):
        s = Scaffold(name, "ATGC", description=f"{name} Homo sapiens chromosome 1")
        assert s.chromosome is False
        assert s.unknown_chromosome is False

    @pytest.mark.parametrize("name", ["chr1", "Chr01", "chr1A", "chr1_hap1", "chrX"])
    def test_chromosome_names_stay_chromosomes(self, name):
        assert Scaffold(name, "ATGC").chromosome is True

    def test_renaming_keeps_input_id_in_description(self):
        s = Scaffold("NC_012007.3", "ATGC", description="NC_012007.3 Malus domestica chromosome 1")
        s.update(new_name="chr01")
        assert s.description == "chr01 NC_012007.3 Malus domestica chromosome 1"
        # renaming again keeps the input id once
        s.update(new_name="chr02")
        assert s.description == "chr02 NC_012007.3 Malus domestica chromosome 1"
        # renaming back to the input id drops it
        s.update(new_name="NC_012007.3")
        assert s.description == "NC_012007.3 Malus domestica chromosome 1"

    def test_renaming_without_description(self):
        s = Scaffold("scf1", "ATGC")
        s.update(new_name="chr00000001")
        assert s.description == "chr00000001 scf1"

    def test_renaming_only_replaces_a_whole_id(self):
        """chr1 must not be read as the start of chr10."""
        s = Scaffold("chr1", "ATGC", description="chr10 something")
        s.update(new_name="chr01")
        assert s.description == "chr10 something"


# ============================================================
# Genome
# ============================================================


class TestClassifyCompartment:
    def test_user_lists_come_first(self):
        assert classify_compartment("chrM", scaffold=Scaffold("chrM", "A"), chloroplast_chroms=("chrM",)) == "chloroplast"
        assert classify_compartment("contig_9", mitochondria_chroms=("contig_9",)) == "mitochondria"

    def test_user_lists_apply_without_auto_detection(self):
        assert classify_compartment("contig_9", chloroplast_chroms=("contig_9",), auto_detection=False) == "chloroplast"
        assert classify_compartment("chrM", scaffold=Scaffold("chrM", "A"), auto_detection=False) == "nuclear"
        assert classify_compartment("chrM", auto_detection=False) == "nuclear"

    def test_scaffold_classification_is_used_when_given(self):
        renamed = Scaffold("chrM", "A")
        renamed.update(new_name="chr0001")
        assert classify_compartment("chr0001", scaffold=renamed) == "mitochondria"
        described = Scaffold("NC_000932.1", "A", description="NC_000932.1 Arabidopsis thaliana chloroplast, complete genome")
        assert classify_compartment("NC_000932.1", scaffold=described) == "chloroplast"
        assert classify_compartment("chr1", scaffold=Scaffold("chr1", "A")) == "nuclear"

    @pytest.mark.parametrize("name, expected", [
        ("chrM", "mitochondria"), ("Mt", "mitochondria"), ("ChrC", "chloroplast"), ("Pt", "chloroplast"),
        ("chr1", "nuclear"), ("scaffold_9", "nuclear"),
    ])
    def test_name_is_classified_without_scaffold(self, name, expected):
        assert classify_compartment(name) == expected

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





# ============================================================
# DAP renaming
# ============================================================

class TestDapRenaming:
    @staticmethod
    def _write_fasta(path, names):
        path.write_text("".join(f">{name}\n{'ACGT' * (10 - i)}\n" for i, name in enumerate(names)))
        return str(path)

    @staticmethod
    def _headers(path):
        return [line[1:].strip() for line in path.read_text().splitlines() if line.startswith(">")]

    def test_scaffolds_stay_scaffolds_after_renaming(self, tmp_path):
        """DAP scaffold names (chr00000001) start with "chr" but are not chromosomes."""
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "chr02", "NW_0001.1", "scf2", "chrM"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap()
        assert g.scaffold_names == {"chr00000001", "chr00000002"}
        assert g.chromosome_names == {"chr01", "chr02"}
        assert g.accessory_chromosome_names == {"chr0001"}
        assert g.suffix == "_dapmod"

    def test_scffree_export_drops_scaffolds(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "chr02", "NW_0001.1", "scf2", "chrM"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path))
        g.remove_organelles(export=True, output_dir=str(tmp_path))
        g.remove_scaffolds(export=True, output_dir=str(tmp_path))
        assert self._headers(tmp_path / "toy_dapmod_organellefree.fasta") == ["chr01", "chr02", "chr00000001 NW_0001.1", "chr00000002 scf2"]
        assert self._headers(tmp_path / "toy_dapmod_scffree_organellefree_no_description.fasta") == ["chr01", "chr02"]

    def test_dap_export_can_drop_descriptions(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "NW_0001.1", "chrM"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path), keep_description=False)
        assert self._headers(tmp_path / "toy_dapmod_no_description.fasta") == ["chr01", "chr00000001", "chr0001"]

    def test_reloaded_dap_genome_keeps_input_ids(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "NW_0001.1", "chrM"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path))
        reloaded = Genome("toy", str(tmp_path / "toy_dapmod.fasta"), quiet=True)
        reloaded.remove_organelles(export=True, output_dir=str(tmp_path))
        # the unknown_dapmod placeholder is never written
        assert self._headers(tmp_path / "toy_dapmod_organellefree.fasta") == ["chr01", "chr00000001 NW_0001.1"]

    def test_dapmod_survives_removing_renamed_organelles(self, tmp_path):
        """When only the organelles were renamed, removing them must not turn the genome into "_dapfit"."""
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "chr02", "chrM", "chrC"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap()
        g.remove_organelles()
        assert g.suffix == "_dapmod_scffree_organellefree"

    def test_reloaded_dap_genome_keeps_classification(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "chr02", "NW_0001.1", "chrM"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path))
        reloaded = Genome("toy", str(tmp_path / "toy_dapmod.fasta"), quiet=True)
        assert reloaded.suffix == "_dapmod"
        assert reloaded.chromosome_names == {"chr01", "chr02"}
        assert reloaded.scaffold_names == {"chr00000001"}
        assert reloaded.accessory_chromosome_names == {"chr0001"}
        assert reloaded.scaffolds["chr0001"].organelle is True

        reloaded.remove_organelles()
        reloaded.remove_scaffolds()
        assert list(reloaded.scaffolds) == ["chr01", "chr02"]
        assert reloaded.suffix == "_dapmod_scffree_organellefree"

    def test_dap_name_patterns_need_dapmod(self):
        """Outside DAP genomes, chr00000001 and chr0001 keep the generic classification."""
        assert Scaffold("chr00000001", "ACGT").chromosome is True
        assert Scaffold("chr0001", "ACGT").organelle is False
        assert Scaffold("chr00000001", "ACGT", dapmod=True).chromosome is False
        assert Scaffold("chr0001", "ACGT", dapmod=True).organelle is True

    @staticmethod
    def _mixed_fasta(tmp_path):
        path = tmp_path / "mixed.fasta"
        path.write_text(">chr1 desc one\nAAAA\n>chr2\nCCCC\n>chr00\nGGGG\n>scf1 some scaffold\nTTTT\n>ChrC\nACGT\n")
        return str(path)

    @staticmethod
    def _flags(genome):
        return {name: (s.chromosome, s.unknown_chromosome, s.organelle, s.mitochondria, s.chloroplast) for name, s in genome.scaffolds.items()}

    def test_single_unknown_chromosome_is_named_chr00(self, tmp_path):
        g = Genome("mixed", self._mixed_fasta(tmp_path), quiet=True)
        equivalences = g.rename_features_dap(return_equivalences=True)
        assert {k: equivalences[k] for k in ("chr1", "chr2", "chr00", "scf1", "ChrC")} == {
            "chr1": "chr01", "chr2": "chr02", "chr00": "chr00", "scf1": "chr00000001", "ChrC": "chr0001"}
        assert g.scaffolds["chr00"].unknown_chromosome is True
        assert g.suffix == "_dapmod_chr00"

    def test_single_unknown_chromosome_name_is_kept_for_chrun(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "un.fa", ["chr1", "chrUn_KI270302v1", "scf1"])
        g = Genome("un", fasta, quiet=True)
        equivalences = g.rename_features_dap(return_equivalences=True)
        assert equivalences["chrUn_KI270302v1"] == "chr00"
        assert equivalences["scf1"] == "chr00000001"

    def test_several_unknown_chromosomes_get_scaffold_names(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "multi.fa", ["chr1", "chrUn", "chr00", "scf1"])
        g = Genome("multi", fasta, quiet=True)
        equivalences = g.rename_features_dap(return_equivalences=True)
        assert [equivalences[k] for k in ("chrUn", "chr00", "scf1")] == ["chr00000001", "chr00000002", "chr00000003"]
        assert g.scaffolds["chr00000001"].unknown_chromosome is True
        assert g.scaffolds["chr00000002"].unknown_chromosome is True
        assert g.scaffolds["chr00000003"].unknown_chromosome is False

    def test_chromosome_numbered_0_keeps_chr00(self, tmp_path):
        """An unknown chromosome must not take the chr00 name of a chromosome numbered 0."""
        fasta = self._write_fasta(tmp_path / "zero.fa", ["Ch0", "Ch1", "chrUn"])
        g = Genome("zero", fasta, quiet=True)
        equivalences = g.rename_features_dap(return_equivalences=True)
        assert equivalences["Ch0"] == "chr00"
        assert equivalences["chrUn"] == "chr00000001"
        assert len(g.scaffolds) == 3

    def test_dapfit_genome_with_chr00_is_unchanged(self, tmp_path, capsys):
        fasta = self._write_fasta(tmp_path / "fit.fa", ["chr01", "chr02", "chr00"])
        g = Genome("fit", fasta, quiet=True)
        assert g.rename_features_dap(return_equivalences=True) is None
        assert "already fit for DAP-Seq analysis" in capsys.readouterr().out
        assert list(g.scaffolds) == ["chr01", "chr02", "chr00"]
        assert g.suffix == "_scffree_dapfit_organellefree_chr00"

    @pytest.mark.parametrize("fasta_text", [
        ">chr1 desc one\nAAAA\n>chr2\nCCCC\n>chr00\nGGGG\n>scf1 some scaffold\nTTTT\n>ChrC\nACGT\n",
        # several unknown chromosomes get scaffold names, recognised from their input ids
        ">chr1\nAAAA\n>chrUn\nCCCC\n>chr00\nGGGG\n>scf1\nTTTT\n>NC_1 mitochondrion genome\nACGT\n>ChrC\nAC\n",
    ])
    @pytest.mark.parametrize("remove_00", [None, False, True])
    def test_reloaded_dap_exports_keep_classification(self, tmp_path, fasta_text, remove_00):
        (tmp_path / "mixed.fasta").write_text(fasta_text)
        g = Genome("mixed", str(tmp_path / "mixed.fasta"), quiet=True)
        out_dir = tmp_path / "out"
        g.rename_features_dap()
        if remove_00 is not None:
            g.remove_scaffolds(remove_00=remove_00)
        g.export(output_dir=str(out_dir), quiet=True)
        reloaded = Genome("mixed", str(out_dir / f"mixed{g.suffix}.fasta"), quiet=True)
        assert reloaded.suffix == g.suffix
        assert self._flags(reloaded) == self._flags(g)
        for genome in (g, reloaded):
            genome.remove_scaffolds(remove_00=True)
        assert reloaded.suffix == g.suffix
        assert self._flags(reloaded) == self._flags(g)

    def test_reloaded_dap_genome_reads_input_ids_of_dap_names_only(self, tmp_path):
        g = Genome("mixed", self._mixed_fasta(tmp_path), quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path / "out"))
        reloaded = Genome("mixed", str(tmp_path / "out" / "mixed_dapmod_chr00.fasta"), quiet=True)
        # chr01 was renamed too, but its second header word is only known to be an id for DAP scaffold/organelle names
        assert {name: s.original_name for name, s in reloaded.scaffolds.items()} == {
            "chr01": "unknown_dapmod_1", "chr02": "unknown_dapmod_2", "chr00": "unknown_dapmod_3",
            "chr00000001": "scf1", "chr0001": "ChrC"}
        reloaded.export(output_dir=str(tmp_path / "again"), quiet=True)
        assert self._headers(tmp_path / "again" / "mixed_dapmod_chr00.fasta") == [
            "chr01 chr1 desc one", "chr02 chr2", "chr00", "chr00000001 scf1 some scaffold", "chr0001 ChrC"]

    def test_dapmod_is_read_from_the_file_name_only(self, tmp_path):
        folder = tmp_path / "results_dapmod"
        folder.mkdir()
        (folder / "toy.fa").write_text(">chr01\nAAAA\n>chr00000001 chrUn\nCCCC\n>chr0001 ChrM\nGGGG\n")
        g = Genome("toy", str(folder / "toy.fa"), quiet=True)
        assert g.dapmod is False
        assert "_dapmod" not in g.suffix
        # outside DAP genomes, the second header word is not an input id
        assert g.scaffolds["chr00000001"].original_name == "chr00000001"
        assert g.scaffolds["chr00000001"].unknown_chromosome is False

    def test_reloaded_dap_genome_without_input_ids(self, tmp_path):
        g = Genome("mixed", self._mixed_fasta(tmp_path), quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path / "out"), keep_description=False)
        reloaded = Genome("mixed", str(tmp_path / "out" / "mixed_dapmod_chr00_no_description.fasta"), quiet=True)
        assert reloaded.scaffolds["chr00000001"].original_name == "unknown_dapmod_4"
        assert reloaded.scaffolds["chr0001"].original_name == "unknown_dapmod_5"
        # the organelle type is lost with the input id, but it is still an organelle
        assert reloaded.scaffolds["chr0001"].organelle is True
        assert reloaded.scaffolds["chr0001"].chloroplast is False


# ============================================================
# Chromosomes named by accession, warnings
# ============================================================

class TestDescriptionChromosomes:
    @pytest.mark.parametrize("description, expected", [
        ("NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly", "1"),
        ("NC_000023.11 Homo sapiens chromosome X, GRCh38.p14 Primary Assembly", "X"),
        ("NT_187361.1 Homo sapiens chromosome 1 unlocalized genomic scaffold", None),
        ("NT_167214.1 Homo sapiens chromosome Un unplaced genomic scaffold", None),
        ("NW_025791756.1 Homo sapiens chromosome 1 genomic patch of type FIX", None),
        ("NC_012920.1 Homo sapiens mitochondrion, complete genome", None),
        ("chromosome_Sym", None),
    ])
    def test_description_chromosome(self, description, expected):
        assert Scaffold.description_chromosome(description) == expected

    def test_accession_named_chromosomes_are_dap_renamed_as_chromosomes(self, tmp_path):
        fasta = tmp_path / "grch.fa"
        fasta.write_text(
            ">NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly\nAAAAAAAA\n"
            ">NT_187361.1 Homo sapiens chromosome 1 unlocalized genomic scaffold\nCCCC\n"
            ">NC_000002.12 Homo sapiens chromosome 2, GRCh38.p14 Primary Assembly\nGGGGGG\n"
            ">NC_012920.1 Homo sapiens mitochondrion, complete genome\nATATAT\n"
        )
        g = Genome("grch", str(fasta), quiet=True)
        assert g.chromosome_names == {"NC_000001.11", "NC_000002.12"}
        assert g.scaffold_names == {"NT_187361.1"}
        equivalences = g.rename_features_dap(return_equivalences=True)
        assert equivalences["NC_000001.11"] == "chr01"
        assert equivalences["NC_000002.12"] == "chr02"
        assert equivalences["NT_187361.1"] == "chr00000001"
        assert equivalences["NC_012920.1"] == "chr0001"


class TestGenomeWarnings:
    def test_warns_when_no_chromosomes_are_left(self, tmp_path):
        fasta = tmp_path / "acc.fa"
        fasta.write_text(">NC_012007.3\nAAAA\n>NC_012008.3\nCCCC\n>MT\nGGGG\n")
        g = Genome("acc", str(fasta), quiet=True)
        with pytest.warns(UserWarning, match="No recognised nuclear chromosomes in acc genome.*only organelles"):
            g.remove_scaffolds()

    def test_organelle_only_scffree_genome_is_not_exported(self, tmp_path):
        fasta = tmp_path / "acc.fa"
        fasta.write_text(">NC_012007.3\nAAAA\n>NC_012008.3\nCCCC\n>MT\nGGGG\n")
        g = Genome("acc", str(fasta), quiet=True)
        out_dir = tmp_path / "out"
        with pytest.warns(UserWarning) as record:
            g.remove_scaffolds(export=True, output_dir=str(out_dir))
        messages = [str(w.message) for w in record]
        assert any("only organelles (MT) are left" in m for m in messages)
        assert any("nothing to export for acc genome" in m for m in messages)
        assert not out_dir.exists() or not any(out_dir.iterdir())

    def test_dap_renaming_warns_when_no_chromosomes_are_found(self, tmp_path):
        fasta = tmp_path / "acc.fa"
        fasta.write_text(">NC_012007.3\nAAAA\n>NC_012008.3\nCCCC\n>MT\nGGGG\n")
        g = Genome("acc", str(fasta), quiet=True)
        with pytest.warns(UserWarning, match="No recognised nuclear chromosomes in acc genome.*DAP renaming"):
            equivalences = g.rename_features_dap(return_equivalences=True)
        assert equivalences["NC_012007.3"] == "chr00000001"
        assert equivalences["MT"] == "chr0001"

    def _mixed_genome(self, tmp_path):
        fasta = tmp_path / "mixed.fa"
        fasta.write_text(">Chr1\nAAAA\n>Chr2\nCCCC\n>chr00\nGGGG\n>scaffold_1\nTTTT\n>ChrM\nAA\n>ChrC\nCC\n")
        return Genome("mixed", str(fasta), quiet=True)

    def test_scffree_names_with_chr00(self, tmp_path):
        g = self._mixed_genome(tmp_path)
        assert g.suffix == "_chr00"
        g.remove_scaffolds(remove_00=False)
        assert g.suffix == "_scffree_chr00"
        assert set(g.scaffolds) == {"Chr1", "Chr2", "chr00", "ChrM", "ChrC"}
        g.remove_scaffolds(remove_00=True)
        assert g.suffix == "_scffree"

    def test_scffree_export_keeping_chr00_does_not_overwrite_full_export(self, tmp_path):
        g = self._mixed_genome(tmp_path)
        g.export(output_dir=str(tmp_path), quiet=True)
        g.remove_scaffolds(export=True, output_dir=str(tmp_path), remove_00=False)
        assert (tmp_path / "mixed_chr00_no_description.fasta").exists()
        assert (tmp_path / "mixed_scffree_chr00_no_description.fasta").exists()

    def test_unchanged_scffree_genome_is_not_exported(self, tmp_path):
        fasta = tmp_path / "only00.fa"
        fasta.write_text(">chr1\nAAAA\n>chr00\nGGGG\n")
        g = Genome("only00", str(fasta), quiet=True)
        assert g.suffix == "_scffree_dapfit_organellefree_chr00"
        with pytest.warns(UserWarning, match="No scaffolds were removed from only00 genome"):
            g.remove_scaffolds(export=True, output_dir=str(tmp_path / "out"), remove_00=False)
        assert not (tmp_path / "out").exists()

    def test_ucsc_unplaced_and_alt_scaffolds_are_removed(self, tmp_path):
        fasta = tmp_path / "ucsc.fa"
        fasta.write_text(">chr1\nAAAA\n>chr1_random\nCC\n>chr4_KI270790v1_alt\nGG\n>chrUn_KI270302v1\nTT\n>chrM\nAA\n")
        g = Genome("ucsc", str(fasta), quiet=True)
        assert g.chromosome_names == {"chr1"}
        g.remove_scaffolds(remove_00=False)
        assert set(g.scaffolds) == {"chr1", "chrUn_KI270302v1", "chrM"}

    def test_dap_renaming_does_not_warn_with_chromosomes(self, tmp_path, recwarn):
        fasta = tmp_path / "toy.fa"
        fasta.write_text(">chr01\nAAAA\n>scf1\nCCCC\n>chrM\nGGGG\n")
        g = Genome("toy", str(fasta), quiet=True)
        g.rename_features_dap()
        assert not [w for w in recwarn if "No recognised nuclear chromosomes" in str(w.message)]

    def test_warns_when_export_overwrites_earlier_export(self, tmp_path):
        fasta = tmp_path / "toy.fa"
        fasta.write_text(">chr1\nAAAA\n")
        g = Genome("toy", str(fasta), quiet=True)
        g.export(output_dir=str(tmp_path), quiet=True)
        with pytest.warns(UserWarning, match="already written earlier in this run"):
            g.export(output_dir=str(tmp_path), quiet=True)

    def test_scaffold_and_organelle_removal_passes_output_arguments(self, tmp_path):
        g = self._mixed_genome(tmp_path)
        g.remove_scaffolds(export=True, remove_organelles=True, use_genome_dir=True, subfolder=True, subfolder_name="sub", filename="clean", extension="fa")
        out = tmp_path / "sub" / "clean.fa"
        assert out.exists()
        assert [line[1:] for line in out.read_text().splitlines() if line.startswith(">")] == ["Chr1", "Chr2"]

        g = self._mixed_genome(tmp_path)
        g.remove_scaffolds(export=True, remove_organelles=True, filepath=str(tmp_path / "exact" / "clean.fasta"))
        assert (tmp_path / "exact" / "clean.fasta").exists()

    def test_unchanged_organelle_free_genome_is_not_exported(self, tmp_path):
        fasta = tmp_path / "noorg.fa"
        fasta.write_text(">chr1\nAAAA\n>scf1\nCCCC\n")
        g = Genome("noorg", str(fasta), quiet=True)
        with pytest.warns(UserWarning, match="No organelles were removed from noorg genome"):
            g.remove_organelles(export=True, output_dir=str(tmp_path / "out"))
        assert not (tmp_path / "out").exists()

    def test_scaffold_and_organelle_removal_exports_without_organelles(self, tmp_path):
        """Only scaffolds are removed, which is still worth exporting."""
        fasta = tmp_path / "noorg.fa"
        fasta.write_text(">chr1\nAAAA\n>scf1\nCCCC\n")
        g = Genome("noorg", str(fasta), quiet=True)
        g.remove_scaffolds(export=True, remove_organelles=True, output_dir=str(tmp_path / "out"))
        assert (tmp_path / "out" / "noorg_scffree_dapfit_organellefree_no_description.fasta").exists()

    def test_scaffold_and_organelle_removal_exports_without_scaffolds(self, tmp_path):
        fasta = tmp_path / "noscf.fa"
        fasta.write_text(">chr1\nAAAA\n>chrM\nCCCC\n")
        g = Genome("noscf", str(fasta), quiet=True)
        g.remove_scaffolds(export=True, remove_organelles=True, output_dir=str(tmp_path / "out"))
        assert (tmp_path / "out" / "noscf_scffree_dapfit_organellefree_no_description.fasta").exists()

    @pytest.mark.parametrize("remove_organelles, message", [
        (False, "No scaffolds were removed from clean genome, so its scaffold-free version"),
        (True, "No scaffolds or organelles were removed from clean genome, so its scaffold- and organelle-free version"),
    ])
    def test_genome_without_scaffolds_is_not_exported(self, tmp_path, remove_organelles, message):
        fasta = tmp_path / "clean.fa"
        fasta.write_text(">chr1\nAAAA\n>chr2\nCCCC\n")
        g = Genome("clean", str(fasta), quiet=True)
        with pytest.warns(UserWarning, match=message):
            g.remove_scaffolds(export=True, remove_organelles=remove_organelles, output_dir=str(tmp_path / "out"))
        assert not (tmp_path / "out").exists()


class TestNoDescriptionTag:
    @staticmethod
    def _export_names(tmp_path, fasta_text, **export_kwargs):
        (tmp_path / "toy.fa").write_text(fasta_text)
        g = Genome("toy", str(tmp_path / "toy.fa"), quiet=True)
        g.export(output_dir=str(tmp_path / "out"), quiet=True, **export_kwargs)
        return sorted(p.name for p in (tmp_path / "out").iterdir())

    def test_dropped_descriptions(self, tmp_path):
        names = self._export_names(tmp_path, ">chr1 some text\nAAAA\n>chr2 more\nCCCC\n", keep_description=False)
        assert names == ["toy_scffree_dapfit_organellefree_no_description.fasta"]

    def test_input_without_descriptions(self, tmp_path):
        names = self._export_names(tmp_path, ">chr1\nAAAA\n>chr2\nCCCC\n")
        assert names == ["toy_scffree_dapfit_organellefree_no_description.fasta"]

    def test_kept_descriptions(self, tmp_path):
        names = self._export_names(tmp_path, ">chr1 some text\nAAAA\n>chr2\nCCCC\n")
        assert names == ["toy_scffree_dapfit_organellefree.fasta"]

    def test_explicit_filename_has_no_tag(self, tmp_path):
        names = self._export_names(tmp_path, ">chr1\nAAAA\n", filename="plain.fasta")
        assert names == ["plain.fasta"]

    def test_renamed_input_ids_count_as_description(self, tmp_path):
        """DAP renaming keeps each renamed sequence's input id in its header."""
        (tmp_path / "toy.fa").write_text(">chr1\nAAAA\n>scf1\nCCCC\n")
        g = Genome("toy", str(tmp_path / "toy.fa"), quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path / "out"))
        assert (tmp_path / "out" / "toy_dapmod_organellefree.fasta").exists()

    def test_tag_is_not_in_genome_suffix(self, tmp_path):
        """Other outputs (e.g. feature sizes) do not depend on FASTA headers."""
        (tmp_path / "toy.fa").write_text(">chr1\nAAAA\n")
        g = Genome("toy", str(tmp_path / "toy.fa"), quiet=True)
        g.export_feature_sizes(output_dir=str(tmp_path / "out"))
        assert "_no_description" not in g.suffix
        assert (tmp_path / "out" / "toy_genome_feature_sizes_scffree_dapfit_organellefree.tsv").exists()


class TestExportBothDescriptions:
    @staticmethod
    def _headers(path):
        return [line[1:] for line in path.read_text().splitlines() if line.startswith(">")]

    @staticmethod
    def _genome(tmp_path, fasta_text):
        (tmp_path / "toy.fa").write_text(fasta_text)
        return Genome("toy", str(tmp_path / "toy.fa"), quiet=True)

    def test_writes_both_versions(self, tmp_path):
        g = self._genome(tmp_path, ">chr1 some text\nAAAA\n>chr2\nCCCC\n")
        g.export(output_dir=str(tmp_path / "out"), quiet=True, keep_description="both")
        out = tmp_path / "out"
        assert sorted(p.name for p in out.iterdir()) == [
            "toy_scffree_dapfit_organellefree.fasta", "toy_scffree_dapfit_organellefree_no_description.fasta"]
        assert self._headers(out / "toy_scffree_dapfit_organellefree.fasta") == ["chr1 some text", "chr2"]
        assert self._headers(out / "toy_scffree_dapfit_organellefree_no_description.fasta") == ["chr1", "chr2"]

    def test_id_only_headers_are_written_once(self, tmp_path, recwarn):
        g = self._genome(tmp_path, ">chr1\nAAAA\n>chr2\nCCCC\n")
        g.export(output_dir=str(tmp_path / "out"), quiet=True, keep_description="both")
        assert [p.name for p in (tmp_path / "out").iterdir()] == ["toy_scffree_dapfit_organellefree_no_description.fasta"]
        assert not [w for w in recwarn if "already written earlier" in str(w.message)]

    @pytest.mark.parametrize("method", ["export", "rename_features_dap", "remove_scaffolds", "remove_organelles"])
    @pytest.mark.parametrize("name_kwargs", [{"filename": "clean.fa"}, {"filepath": "out/clean.fa"}])
    def test_explicit_names_are_rejected(self, tmp_path, method, name_kwargs):
        g = self._genome(tmp_path, ">chr1 some text\nAAAA\n>scf1\nCCCC\n>chrM\nGGGG\n")
        name_kwargs = {k: str(tmp_path / v) if k == "filepath" else v for k, v in name_kwargs.items()}
        export_kwargs = {} if method == "export" else {"export": True}
        with pytest.raises(ValueError, match="cannot be used with filepath or filename"):
            getattr(g, method)(output_dir=str(tmp_path / "out"), keep_description="both", **export_kwargs, **name_kwargs)
        # rejected before the genome is changed or anything is written
        assert list(g.scaffolds) == ["chr1", "scf1", "chrM"]
        assert not (tmp_path / "out").exists()

    def test_explicit_names_are_fine_without_export(self, tmp_path):
        """Output arguments are ignored when nothing is exported."""
        g = self._genome(tmp_path, ">chr1 some text\nAAAA\n>scf1\nCCCC\n")
        g.remove_scaffolds(filename="clean.fa", keep_description="both")
        assert list(g.scaffolds) == ["chr1"]

    def test_empty_genome_warns_once(self, tmp_path):
        g = self._genome(tmp_path, ">chr1 some text\nAAAA\n")
        g.remove_features({"chr1"})
        with pytest.warns(UserWarning) as record:
            g.export(output_dir=str(tmp_path / "out"), quiet=True, keep_description="both")
        assert sum("nothing to export" in str(w.message) for w in record) == 1

    def test_invalid_value(self, tmp_path):
        g = self._genome(tmp_path, ">chr1\nAAAA\n")
        with pytest.raises(ValueError, match="keep_description must be True, False or 'both'"):
            g.export(output_dir=str(tmp_path / "out"), quiet=True, keep_description="yes")

    def test_dap_renaming(self, tmp_path):
        g = self._genome(tmp_path, ">chr1\nAAAA\n>scf1 some scaffold\nCCCC\n")
        g.rename_features_dap(export=True, output_dir=str(tmp_path / "out"), keep_description="both")
        out = tmp_path / "out"
        assert self._headers(out / "toy_dapmod_organellefree.fasta") == ["chr01 chr1", "chr00000001 scf1 some scaffold"]
        assert self._headers(out / "toy_dapmod_organellefree_no_description.fasta") == ["chr01", "chr00000001"]

    def test_scaffold_and_organelle_removal(self, tmp_path):
        fasta_text = ">chr1 one\nAAAA\n>scf1\nCCCC\n>chrM mito\nGGGG\n"
        out = tmp_path / "out"
        g = self._genome(tmp_path, fasta_text)
        g.remove_scaffolds(export=True, output_dir=str(out / "scf"), keep_description="both")
        # chrM is not DAP-fit, so _dapfit only comes once organelles are removed
        assert sorted(p.name for p in (out / "scf").iterdir()) == ["toy_scffree.fasta", "toy_scffree_no_description.fasta"]

        g = self._genome(tmp_path, fasta_text)
        g.remove_organelles(export=True, output_dir=str(out / "org"), keep_description="both")
        assert sorted(p.name for p in (out / "org").iterdir()) == ["toy_organellefree.fasta", "toy_organellefree_no_description.fasta"]

        g = self._genome(tmp_path, fasta_text)
        g.remove_scaffolds(export=True, remove_organelles=True, output_dir=str(out / "both"), keep_description="both")
        assert sorted(p.name for p in (out / "both").iterdir()) == [
            "toy_scffree_dapfit_organellefree.fasta", "toy_scffree_dapfit_organellefree_no_description.fasta"]


class TestExportPaths:
    @staticmethod
    def _genome(tmp_path, fasta_text, name="toy"):
        (tmp_path / f"{name}.fa").write_text(fasta_text)
        return Genome(name, str(tmp_path / f"{name}.fa"), quiet=True)

    @pytest.mark.parametrize("keep_description, expected", [
        (True, {"description": "toy_scffree_dapfit_organellefree.fasta"}),
        (False, {"no_description": "toy_scffree_dapfit_organellefree_no_description.fasta"}),
        ("both", {"description": "toy_scffree_dapfit_organellefree.fasta", "no_description": "toy_scffree_dapfit_organellefree_no_description.fasta"}),
    ])
    def test_export_returns_absolute_paths_by_header_type(self, tmp_path, monkeypatch, keep_description, expected):
        monkeypatch.chdir(tmp_path)
        g = self._genome(tmp_path, ">chr1 some text\nAAAA\n")
        paths = g.export(output_dir="out", quiet=True, keep_description=keep_description)
        assert paths == {header_type: tmp_path / "out" / name for header_type, name in expected.items()}
        assert all(path.is_absolute() and path.exists() for path in paths.values())

    @pytest.mark.parametrize("keep_description, expected_keys", [(True, ["description"]), ("both", ["description", "no_description"])])
    def test_id_only_headers_point_every_type_to_one_file(self, tmp_path, keep_description, expected_keys):
        g = self._genome(tmp_path, ">chr1\nAAAA\n")
        paths = g.export(output_dir=str(tmp_path / "out"), quiet=True, keep_description=keep_description)
        assert list(paths) == expected_keys
        assert set(paths.values()) == {tmp_path / "out" / "toy_scffree_dapfit_organellefree_no_description.fasta"}
        assert len(list((tmp_path / "out").iterdir())) == 1

    def test_nothing_to_export(self, tmp_path):
        g = self._genome(tmp_path, ">chr1\nAAAA\n")
        g.remove_features({"chr1"})
        with pytest.warns(UserWarning, match="nothing to export"):
            assert g.export(output_dir=str(tmp_path / "out"), quiet=True) == {}

    def test_symlinked_paths_are_not_rewritten(self, tmp_path):
        (tmp_path / "real").mkdir()
        (tmp_path / "link").symlink_to(tmp_path / "real")
        g = self._genome(tmp_path, ">chr1 some text\nAAAA\n")
        paths = g.export(output_dir=str(tmp_path / "link"), quiet=True)
        assert paths == {"description": tmp_path / "link" / "toy_scffree_dapfit_organellefree.fasta"}
        # the overwrite check still sees through the link
        with pytest.warns(UserWarning, match="already written earlier in this run"):
            g.export(output_dir=str(tmp_path / "real"), quiet=True)

    @pytest.mark.parametrize("method", ["remove_scaffolds", "remove_organelles"])
    def test_removal_returns_exported_paths(self, tmp_path, method):
        g = self._genome(tmp_path, ">chr1 one\nAAAA\n>scf1\nCCCC\n>chrM\nGGGG\n")
        paths = getattr(g, method)(export=True, output_dir=str(tmp_path / "out"), keep_description="both")
        assert set(paths) == {"description", "no_description"}
        assert paths["description"] != paths["no_description"]
        assert all(path.exists() for path in paths.values())

    def test_scaffold_and_organelle_removal_returns_exported_paths(self, tmp_path):
        g = self._genome(tmp_path, ">chr1 one\nAAAA\n>scf1\nCCCC\n>chrM\nGGGG\n")
        paths = g.remove_scaffolds(export=True, remove_organelles=True, output_dir=str(tmp_path / "out"))
        assert paths == {"description": tmp_path / "out" / "toy_scffree_dapfit_organellefree.fasta"}

    @pytest.mark.parametrize("method", ["remove_scaffolds", "remove_organelles"])
    def test_removal_without_export_returns_nothing(self, tmp_path, method):
        g = self._genome(tmp_path, ">chr1 one\nAAAA\n>scf1\nCCCC\n>chrM\nGGGG\n")
        assert getattr(g, method)(output_dir=str(tmp_path / "out")) == {}

    def test_refused_scaffold_free_export_returns_nothing(self, tmp_path):
        g = self._genome(tmp_path, ">NC_012007.3\nAAAA\n>MT\nGGGG\n")
        with pytest.warns(UserWarning, match="nothing to export"):
            assert g.remove_scaffolds(export=True, output_dir=str(tmp_path / "out")) == {}

    @pytest.mark.parametrize("method, fasta_text", [
        ("remove_scaffolds", ">chr1 one\nAAAA\n>chrM\nGGGG\n"),
        ("remove_organelles", ">chr1 one\nAAAA\n>scf1\nCCCC\n"),
    ])
    def test_skipped_export_returns_paths_written_this_run(self, tmp_path, method, fasta_text):
        out = str(tmp_path / "out")
        g = self._genome(tmp_path, fasta_text)
        exported = g.export(output_dir=out, quiet=True, keep_description="both")
        mtimes = {path: path.stat().st_mtime_ns for path in exported.values()}
        with pytest.warns(UserWarning, match="was not exported"):
            paths = getattr(g, method)(export=True, output_dir=out, keep_description="both")
        assert paths == exported
        # not written again by the skipped export
        assert {path: path.stat().st_mtime_ns for path in paths.values()} == mtimes

    @pytest.mark.parametrize("method, fasta_text", [
        ("remove_scaffolds", ">chr1 one\nAAAA\n>chrM\nGGGG\n"),
        ("remove_organelles", ">chr1 one\nAAAA\n>scf1\nCCCC\n"),
    ])
    def test_skipped_export_returns_nothing_not_written_this_run(self, tmp_path, method, fasta_text):
        # a stale file from an earlier run is not reported
        (tmp_path / "out").mkdir()
        g = self._genome(tmp_path, fasta_text, name=f"stale_{method}")
        (tmp_path / "out" / f"stale_{method}{g.suffix}.fasta").write_text(">old\nA\n")
        with pytest.warns(UserWarning, match="was not written earlier in this run, so it is not returned"):
            assert getattr(g, method)(export=True, output_dir=str(tmp_path / "out")) == {}

    def test_skipped_export_checks_each_path(self, tmp_path):
        out = str(tmp_path / "out")
        g = self._genome(tmp_path, ">chr1 one\nAAAA\n>chrM\nGGGG\n")
        exported = g.export(output_dir=out, quiet=True, keep_description=True)
        with pytest.warns(UserWarning, match="_no_description.fasta was not written earlier in this run"):
            paths = g.remove_scaffolds(export=True, output_dir=out, keep_description="both")
        assert paths == exported


class TestGenomeAttributes:
    def test_other_scaffolds(self, tmp_path):
        fasta = tmp_path / "toy.fa"
        fasta.write_text(">chr1\nAAAA\n>chr00\nCCCC\n>scf1\nGGGG\n")
        g = Genome("toy", str(fasta), quiet=True)
        assert g.non_chromosomal_scaffolds is True
        assert g.other_scaffolds is True
        g.remove_scaffolds(remove_00=False)
        assert g.non_chromosomal_scaffolds is True
        assert g.other_scaffolds is False

    def test_export_replaces_tabs_in_descriptions(self, tmp_path):
        fasta = tmp_path / "gwh.fa"
        fasta.write_text(">GWHAAAA00000001\tOriSeqID=chr1\tLen=4\nAAAA\n>GWHAAAA00000002\tOriSeqID=chr2\tLen=4\nCCCC\n")
        g = Genome("gwh", str(fasta), quiet=True, gwh=True)
        assert "\t" in g.scaffolds["chr1"].description
        g.export(output_dir=str(tmp_path / "out"), quiet=True)
        exported = tmp_path / "out" / f"gwh{g.suffix}.fasta"
        headers = [line[1:] for line in exported.read_text().splitlines() if line.startswith(">")]
        assert headers == ["chr1 GWHAAAA00000001 OriSeqID=chr1 Len=4", "chr2 GWHAAAA00000002 OriSeqID=chr2 Len=4"]
        reloaded = Genome("gwh", str(exported), quiet=True, gwh=True)
        assert list(reloaded.scaffolds) == ["chr1", "chr2"]
