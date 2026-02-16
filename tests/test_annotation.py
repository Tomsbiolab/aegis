"""
Tests for aegis.annotation — parsing helpers and Annotation class.
"""

import os
import pytest
from pathlib import Path

from aegis.annotation import (
    read_file_with_fallback,
    detect_file_format,
    parse_gtf_attributes,
    format_gff3_attributes,
    sort_and_update_genes,
    Annotation,
    default_features,
)


# ============================================================
# read_file_with_fallback
# ============================================================

class TestReadFileWithFallback:
    def test_utf8_file(self, tmp_path):
        f = tmp_path / "test.gff3"
        f.write_text("##gff-version 3\nchr1\taegis\tgene\t1\t100\t.\t+\t.\tID=g1\n", encoding="utf-8")
        enc = read_file_with_fallback(str(f))
        assert enc == "utf-8"

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
        attrs = parse_gtf_attributes('gene_id "GENE1"; transcript_id "TX1";')
        assert attrs["gene_id"] == "GENE1"
        assert attrs["transcript_id"] == "TX1"

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
        assert len(annot.all_gene_ids) > 0
        assert len(annot.chrs) > 0
