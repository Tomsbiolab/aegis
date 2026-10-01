"""
Tests for aegis.gene — the Gene class.
"""

import pytest
from aegis.feature import Feature

# ============================================================
# __init__
# ============================================================

class TestGeneInit:
    def test_inherits_feature_properties(self, make_gene):
        g = make_gene()
        assert g.id == "gene001"
        assert g.ch == "chr1"
        assert g.size == 1901  # (2000 - 100) + 1

    def test_pseudogene_flag(self, make_gene):
        g = make_gene(pseudogene=True)
        assert g.pseudogene is True

    def test_transposable_flag(self, make_gene):
        g = make_gene(transposable=True)
        assert g.transposable is True

    def test_transcripts_dict_initially_empty(self, make_gene):
        g = make_gene()
        assert isinstance(g.transcripts, dict)
        assert len(g.transcripts) == 0

    def test_coding_defaults(self, make_gene):
        g = make_gene()
        assert g.coding == 0
        assert g.noncoding == 0


# ============================================================
# rename
# ============================================================

class TestGeneRename:
    def test_rename_with_prefix(self, make_gene):
        """Gene.rename with a prefix generates: {prefix}{ch}g{count:0{digits}d}"""
        g = make_gene()
        g.rename(count=5, prefix="VIT")
        # Format: VIT + chr1 + g + 00005
        assert g.id == "VITchr1g00005"
        assert g.renamed is True

    def test_rename_with_prefix_and_suffix(self, make_gene):
        """Gene.rename with prefix+suffix generates: {prefix}{ch}g{count:0{digits}d}{sep}{suffix}"""
        g = make_gene()
        g.rename(count=3, prefix="VIT", suffix="v1")
        assert g.id == "VITchr1g00003_v1"
        assert g.renamed is True

    def test_rename_custom_digits(self, make_gene):
        g = make_gene()
        g.rename(count=1, digits=3, prefix="AT")
        assert g.id == "ATchr1g001"

    def test_rename_no_prefix_does_not_change_id(self, make_gene):
        """Without a prefix, rename doesn't change the id formatting."""
        g = make_gene()
        original_id = g.id
        g.rename(count=1)
        # Without prefix, the if-prefix branch is not entered
        assert g.id == original_id

    def test_rename_remove_point_suffix(self, make_gene):
        g = make_gene(feature_id="gene001.v2")
        g.rename(count=1, remove_point_suffix=True)
        assert "." not in g.id.split("g")[0]  # point suffix removed


# ============================================================
# obtain_base_id
# ============================================================

class TestGeneObtainBaseId:
    def test_base_id_strips_gene_prefix(self, make_gene):
        """obtain_base_id sets self.base_id by stripping 'gene' prefix."""
        g = make_gene(feature_id="gene001")
        # "gene001" starts with "gene" → base_id = id[4:] = "001"
        assert g.base_id == "001"

    def test_base_id_with_gene_suffix(self, make_gene):
        g = make_gene(feature_id="Vitvi01g00123_gene")
        assert g.base_id == "Vitvi01g00123"

    def test_base_id_no_gene_in_name(self, make_gene):
        g = make_gene(feature_id="AT1G00100")
        assert g.base_id == "AT1G00100"

    def test_base_id_original_flag(self, make_gene):
        g = make_gene(feature_id="gene001")
        assert g.base_id == "001"
        assert g.original_base_id == "001"


# ============================================================
# sort_transcripts
# ============================================================

class TestGeneSortTranscripts:
    def test_sort_by_start(self, make_gene, make_transcript):
        g = make_gene(start=100, end=5000)
        g.transcripts["t2"] = make_transcript("t2", start=500, end=2000)
        g.transcripts["t1"] = make_transcript("t1", start=100, end=1500)
        g.sort_transcripts()
        ids = list(g.transcripts.keys())
        starts = [g.transcripts[k].start for k in ids]
        assert starts == sorted(starts)


# ============================================================
# clear_UTRs
# ============================================================

class TestGeneClearUTRs:
    def test_clear_utrs_delegates_to_transcripts(self, make_gene, make_transcript):
        g = make_gene()
        t = make_transcript()
        g.transcripts["t1"] = t
        g.clear_UTRs()
        assert t.temp_UTRs is None


# ============================================================
# __str__
# ============================================================

class TestGeneStr:
    def test_str_with_names(self, make_gene):
        """Gene.__str__ includes names/symbols if present."""
        g = make_gene()
        result = str(g)
        # Gene has Name=TestGene → names=['TestGene'], no symbols
        assert "gene001" in result
        assert "TestGene" in result

    def test_str_no_names_or_symbols(self, make_gene):
        g = make_gene(feature_id="gene001", attributes={"Name":[], "Symbol":[]})
        result = str(g)
        assert result == "gene001"


# ============================================================
# longer_CDS
# ============================================================

class TestGeneLongerCDS:
    def test_longer_cds_returns_none_without_main_transcripts(self, make_gene, make_transcript):
        """longer_CDS returns None when no main transcripts exist."""
        g1 = make_gene(feature_id="g1")
        g2 = make_gene(feature_id="g2")
        t1 = make_transcript("t1", start=100, end=1000)
        t2 = make_transcript("t2", start=100, end=2000)
        g1.transcripts["t1"] = t1
        g2.transcripts["t2"] = t2
        # No transcript has main=True, so it returns None
        result = g1.longer_CDS(g2)
        assert result is None


# ============================================================
# compare_protein_blast_hits
# ============================================================

class TestGeneCompareProteinBlastHits:
    def test_compare_protein_blast_hits_quiet(self, make_gene, make_transcript, capsys):
        g1 = make_gene(feature_id="g1")
        g2 = make_gene(feature_id="g2")
        t1 = make_transcript("t1", start=100, end=1000)
        t2 = make_transcript("t2", start=100, end=2000)
        g1.transcripts["t1"] = t1
        g2.transcripts["t2"] = t2

        res = g1.compare_protein_blast_hits(g2, source_priority=["sp"], quiet=True)
        assert res is None
        captured = capsys.readouterr()
        assert "Warning" not in captured.out

    def test_compare_protein_blast_hits_not_quiet(self, make_gene, make_transcript, capsys):
        g1 = make_gene(feature_id="g1")
        g2 = make_gene(feature_id="g2")
        t1 = make_transcript("t1", start=100, end=1000)
        t2 = make_transcript("t2", start=100, end=2000)
        g1.transcripts["t1"] = t1
        g2.transcripts["t2"] = t2

        res = g1.compare_protein_blast_hits(g2, source_priority=["sp"], quiet=False)
        assert res is None
        captured = capsys.readouterr()
        assert "Warning: g1 and g2 genes have no associated proteins." in captured.out


# ============================================================
# combine_transcripts
# ============================================================

class MockScaffold:
    def __init__(self, seq: str):
        self.seq = seq


class MockGenome:
    def __init__(self, seq_dict: dict[str, str]):
        self.name = "mock_genome"
        self.scaffolds = {k: MockScaffold(v) for k, v in seq_dict.items()}


class TestGeneCombineTranscripts:
    @pytest.fixture
    def setup_mock_genome(self):
        saved = Feature._ACTIVE_GENOME
        def _activate(seq: str, chrom: str = "chr1"):
            Feature._ACTIVE_GENOME = MockGenome({chrom: seq})
        yield _activate
        Feature._ACTIVE_GENOME = saved

    def test_combine_transcripts_orf_or_end_fallback(self, setup_mock_genome, make_gene, make_transcript, make_exon):
        # 30bp sequence without start/stop codons
        seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        setup_mock_genome(seq)

        g = make_gene(feature_id="g1", start=1, end=30, strand="+")
        t1 = make_transcript(feature_id="t1", start=1, end=20, strand="+")
        e1 = make_exon(feature_id="e1", start=1, end=20, strand="+")
        t1.exons = [e1]
        t1.main = True
        t1.coding = True

        t2 = make_transcript(feature_id="t2", start=10, end=30, strand="+")
        e2 = make_exon(feature_id="e2", start=10, end=30, strand="+")
        t2.exons = [e2]

        g.transcripts = {"t1": t1, "t2": t2}
        g.coding = True

        g.combine_transcripts(redetect_CDS=True, respect_non_coding=True, quiet=True)

        assert len(g.transcripts) == 1
        combined_t = list(g.transcripts.values())[0]
        assert combined_t.id == "g1_t001"
        assert combined_t.coding is True
        assert "g1_t001_CDS1" in combined_t.CDSs
        assert combined_t.CDSs["g1_t001_CDS1"].protein is not None
        assert combined_t.CDS_size == 30
