"""
Tests for aegis.annotation — parsing helpers and Annotation class.
"""

import os
import pytest
from pathlib import Path

from aegis.transcript import Transcript
from aegis.annotation import Annotation, default_features
from aegis.utils.gtf_gff import parse_gtf_attributes, format_gff3_attributes, convert_gtf_to_gff3, detect_file_format
from aegis.utils.misc import read_file_with_fallback
from aegis.utils.genefunctions import sort_and_update_genes


# ============================================================
# read_file_with_fallback
# ============================================================

class TestReadFileWithFallback:
    def test_utf8_file(self, tmp_path):
        f = tmp_path / "test.gff3"
        f.write_text("##gff-version 3\nchr1\taegis\tgene\t1\t100\t.\t+\t.\tID=g1\n", encoding="utf-8")
        enc = read_file_with_fallback(str(f))
        assert enc in ("utf-8", "latin-1")

    def test_latin1_file(self, tmp_path):
        f = tmp_path / "test_latin.gff3"
        f.write_bytes("##gff-version 3\nchr1\taegis\tgene\t1\t100\t.\t+\t.\tID=gé1\n".encode("latin-1"))
        enc = read_file_with_fallback(str(f))
        assert enc in ("utf-8", "latin-1")

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(Exception):
            read_file_with_fallback(str(tmp_path / "nonexistent.gff3"))


# ============================================================
# detect_file_format
# ============================================================

class TestDetectFileFormat:
    def test_gff3_with_header(self, tmp_path):
        f = tmp_path / "test.gff3"
        f.write_text("##gff-version 3\nchr1\taegis\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        fmt = detect_file_format(str(f), "utf-8")
        assert fmt == "gff3"

    def test_gff3_without_header(self, tmp_path):
        f = tmp_path / "test.gff3"
        f.write_text("chr1\taegis\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        fmt = detect_file_format(str(f), "utf-8")
        assert fmt == "gff3"

    def test_gtf_format(self, tmp_path):
        f = tmp_path / "test.gtf"
        f.write_text('chr1\taegis\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
        fmt = detect_file_format(str(f), "utf-8")
        assert fmt == "gtf"


# ============================================================
# parse_gtf_attributes
# ============================================================

class TestParseGtfAttributes:
    def test_standard_gtf(self):
        attrs = parse_gtf_attributes('gene_id "GENE1"; transcript_id "T1";')
        assert attrs["gene_id"] == "GENE1"
        assert attrs["transcript_id"] == "T1"

    def test_empty_string(self):
        attrs = parse_gtf_attributes("")
        assert attrs == {}

    def test_single_attribute(self):
        attrs = parse_gtf_attributes('gene_id "ABC123";')
        assert attrs["gene_id"] == "ABC123"


# ============================================================
# format_gff3_attributes
# ============================================================

class TestFormatGff3Attributes:
    def test_gene_format(self):
        attrs = {"gene_id": "g1", "gene_name": "MyGene"}
        result = format_gff3_attributes(attrs, "gene")
        assert "ID=g1" in result
        assert "Symbol=MyGene" in result

    def test_transcript_format(self):
        attrs = {"transcript_id": "t1", "gene_id": "g1"}
        result = format_gff3_attributes(attrs, "mRNA")
        assert "ID=t1" in result
        assert "Parent=g1" in result

    def test_cds_format(self):
        attrs = {"transcript_id": "t1"}
        result = format_gff3_attributes(attrs, "CDS")
        assert "Parent=t1" in result

    def test_exon_with_number(self):
        attrs = {"transcript_id": "t1", "exon_number": "3"}
        result = format_gff3_attributes(attrs, "exon")
        assert "Parent=t1" in result
        assert "ID=t1_e3" in result


# ============================================================
# sort_and_update_genes
# ============================================================

class TestSortAndUpdateGenes:
    def test_sort_genes(self):
        from aegis.gene import Gene
        g1 = Gene(False, False, "gB", "chr1", "aegis", "gene", "+", 500, 1000, ".", ".", "ID=gB")
        g2 = Gene(False, False, "gA", "chr1", "aegis", "gene", "+", 100, 400, ".", ".", "ID=gA")
        genes_dict = {"gB": g1, "gA": g2}
        ch, sorted_dict = sort_and_update_genes("chr1", genes_dict)
        assert ch == "chr1"
        ids = list(sorted_dict.keys())
        starts = [sorted_dict[k].start for k in ids]
        assert starts == sorted(starts)


# ============================================================
# Annotation — integration with small GFF3
# ============================================================

class TestAnnotationSmallGFF3:
    def test_load_minimal_gff3(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        assert annot.name == "sample"
        assert len(annot.all_gene_ids) == 1
        assert "gene1" in annot.all_gene_ids
        assert len(annot.all_transcript_ids) >= 1

    def test_chromosomes_loaded(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        assert "chr1" in annot.chrs

    def test_features_counted(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        assert "gene" in annot.features
        assert "mRNA" in annot.features

    def test_transcript_has_parent(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        assert len(gene.transcripts) >= 1
        t = list(gene.transcripts.values())[0]
        assert "gene1" in t.parents

    def test_multi_gene_gff3(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 2
        assert "geneA" in annot.all_gene_ids
        assert "geneB" in annot.all_gene_ids
        assert "chr1" in annot.chrs
        assert "chr2" in annot.chrs

    def test_copy(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot2 = annot.copy()
        annot2.name = "changed"
        assert annot.name == "sample"


# ============================================================
# Annotation — slow integration test with real test_data
# ============================================================

class TestAnnotationRealData:
    def test_load_grapevine(self, test_data_dir):
        gff_path = str(test_data_dir / "grapevine_v5.1.gff3")
        if not os.path.exists(gff_path):
            pytest.skip("test_data not available")
        annot = Annotation(gff_path, quiet=True)
        assert len(annot.all_gene_ids) == 1820
        assert len(annot.chrs) == 1


# ============================================================
# convert_gtf_to_gff3
# ============================================================

class TestConvertGtfToGff3:
    def test_convert_gtf_basic(self, tmp_path):
        gtf_file = tmp_path / "test.gtf"
        gff_file = tmp_path / "test.gff3"
        gtf_content = (
            "##sequence-region chr1 1 1000\n"
            "chr1\taegis\tgene\t100\t500\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MyGene\";\n"
            "chr1\taegis\ttranscript\t100\t500\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"T1\"; transcript_biotype \"mRNA\";\n"
            "chr1\taegis\texon\t100\t200\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"T1\";\n"
            "chr1\taegis\texon\t300\t500\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"T1\";\n"
        )
        gtf_file.write_text(gtf_content)

        convert_gtf_to_gff3(str(gtf_file), str(gff_file), "utf-8", quiet=True)

        gff_content = gff_file.read_text()
        # Verify headers and format
        assert "##gff-version 3" in gff_content
        # GTF parser drops ## lines so sequence-region won't be in output, except for ##gff-version 3
        assert "ID=GENE1" in gff_content 
        assert "Parent=GENE1" in gff_content
        assert "ID=T1" in gff_content


# ============================================================
# Annotation unique IDs
# ============================================================

class TestAnnotationUniqueIDs:
    def test_get_unique_gene_id(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        # Assuming sample has "gene1"
        annot.all_gene_ids = {"gene1", "gene1_1"}
        new_id = annot._get_unique_gene_id("gene1")
        assert new_id == "gene1_2"

        new_id2 = annot._get_unique_gene_id("gene2")
        assert new_id2 == "gene2"

    def test_get_unique_transcript_id(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.all_transcript_ids = {"t1", "t1_1"}
        new_id = annot._get_unique_transcript_id("t1")
        assert new_id == "t1_2"


# ============================================================
# Annotation marking functions
# ============================================================

class TestAnnotationMarkingFunctions:
    def test_mark_transposable_element_genes(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        te_file = tmp_path / "te_genes.txt"
        te_file.write_text("geneA\n")
        
        annot.mark_transposable_element_genes(str(te_file))
        assert annot.chrs["chr1"]["geneA"].transposable is True
        assert annot.chrs["chr2"]["geneB"].transposable is False

    def test_mark_rRNA_transcripts(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        tx_A = list(annot.chrs["chr1"]["geneA"].transcripts.values())[0]
        tx_A_id = tx_A.id

        rrna_file = tmp_path / "rrna.txt"
        rrna_file.write_text(f"{tx_A_id}\n")

        annot.mark_rRNA_transcripts(str(rrna_file), clean=False)
        assert tx_A.feature == "rRNA"
        tx_B = list(annot.chrs["chr2"]["geneB"].transcripts.values())[0]
        assert not tx_B.feature == "rRNA"

    def test_remove_other_mRNA_transcripts_from_rRNA_genes(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        gene1 = annot.chrs["chr1"]["gene1"]
        gene1.feature = "rRNA_gene"
        rrna_t = Transcript("rRNA1", "chr1", "aegis", "rRNA", "+", 10, 50, ".", ".", "ID=rRNA1;Parent=gene1")
        gene1.transcripts["rRNA1"] = rrna_t

        assert len(gene1.transcripts) == 2
        annot.remove_other_mRNA_transcripts_from_rRNA_genes()
        assert "mRNA1" not in gene1.transcripts
        assert "rRNA1" in gene1.transcripts
