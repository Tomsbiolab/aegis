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
        assert self._headers(tmp_path / "toy_dapmod_scffree_organellefree.fasta") == ["chr01", "chr02"]

    def test_dap_export_can_drop_descriptions(self, tmp_path):
        fasta = self._write_fasta(tmp_path / "toy.fa", ["chr01", "NW_0001.1", "chrM"])
        g = Genome("toy", fasta, quiet=True)
        g.rename_features_dap(export=True, output_dir=str(tmp_path), keep_description=False)
        assert self._headers(tmp_path / "toy_dapmod.fasta") == ["chr01", "chr00000001", "chr0001"]

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
        assert (tmp_path / "mixed_chr00.fasta").exists()
        assert (tmp_path / "mixed_scffree_chr00.fasta").exists()

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
