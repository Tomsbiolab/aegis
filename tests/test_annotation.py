"""
Tests for aegis.annotation — parsing helpers and Annotation class.
"""

import os
import warnings
import pytest
from pathlib import Path
import pytest
from aegis.gene import Gene

from aegis.transcript import Transcript

from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.utils.gtf_gff import parse_gtf_attributes, format_gff3_attributes, convert_gtf_to_gff3, detect_file_format
from aegis.utils.misc import read_file_with_fallback
from aegis.utils.genefunctions import sort_and_update_genes

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"

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

    def test_ascii_file(self, tmp_path):
        f = tmp_path / "ascii.gff3"
        f.write_text("chr1\taegis\tgene\t1\t100\t.\t+\t.\tID=g1\n", encoding="ascii")
        enc = read_file_with_fallback(str(f))
        assert enc in ("utf-8", "latin-1", "ascii")

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.gff3"
        f.write_text("")
        enc = read_file_with_fallback(str(f))
        assert enc is not None


# ============================================================
# detect_file_format
# ============================================================

class TestDetectFileFormat:
    def test_gff3_with_header(self):
        f = str(TEST_DATA_DIR / "input/annotation/detect_gff3_with_header.gff3")
        fmt = detect_file_format(f, "utf-8")
        assert fmt == "gff3"

    def test_gff3_without_header(self):
        f = str(TEST_DATA_DIR / "input/annotation/detect_gff3_without_header.gff3")
        fmt = detect_file_format(f, "utf-8")
        assert fmt == "gff3"

    def test_gtf_format(self):
        f = str(TEST_DATA_DIR / "input/annotation/detect_gtf.gtf")
        fmt = detect_file_format(f, "utf-8")
        assert fmt == "gtf"

    def test_gff3_with_comment_and_blank_lines(self):
        f = str(TEST_DATA_DIR / "input/annotation/detect_gff3_commented.gff3")
        fmt = detect_file_format(f, "utf-8")
        assert fmt == "gff3"

    def test_gtf_without_header(self):
        f = str(TEST_DATA_DIR / "input/annotation/detect_gtf_no_header.gtf")
        fmt = detect_file_format(f, "utf-8")
        assert fmt == "gtf"

    def test_gtf_keyless_attributes(self, keyless_attributes_gtf_file):
        fmt = detect_file_format(keyless_attributes_gtf_file, "utf-8")
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

    def test_multiple_attributes_with_spaces(self):
        attrs = parse_gtf_attributes(
            'gene_id "GENE1"; transcript_id "T1"; gene_name "My Gene"; exon_number "2";'
        )
        assert attrs["gene_id"] == "GENE1"
        assert attrs["gene_name"] == "My Gene"
        assert attrs["exon_number"] == "2"

    def test_trailing_whitespace(self):
        attrs = parse_gtf_attributes('gene_id "GENE1";  ')
        assert attrs["gene_id"] == "GENE1"

    def test_keyless_attribute_unquoted(self):
        attrs = parse_gtf_attributes("PITA_19277")
        assert attrs["id"] == "PITA_19277"

    def test_keyless_attribute_quoted(self):
        attrs = parse_gtf_attributes('"PITA_19277"')
        assert attrs["id"] == "PITA_19277"

    def test_keyless_attribute_custom_default_key(self):
        attrs = parse_gtf_attributes("PITA_19277", default_key="gene_id")
        assert attrs["gene_id"] == "PITA_19277"


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

    def test_five_prime_utr_format(self):
        attrs = {"transcript_id": "t1"}
        result = format_gff3_attributes(attrs, "five_prime_utr")
        assert "Parent=t1" in result

    def test_three_prime_utr_format(self):
        attrs = {"transcript_id": "t1"}
        result = format_gff3_attributes(attrs, "three_prime_utr")
        assert "Parent=t1" in result


# ============================================================
# sort_and_update_genes
# ============================================================

class TestSortAndUpdateGenes:
    def test_sort_genes(self):
        g1 = Gene(False, False, "gB", "chr1", "aegis", "gene", "+", 500, 1000, ".")
        g2 = Gene(False, False, "gA", "chr1", "aegis", "gene", "+", 100, 400, ".")
        genes_dict = {"gB": g1, "gA": g2}
        ch, sorted_dict = sort_and_update_genes("chr1", genes_dict)
        assert ch == "chr1"
        ids = list(sorted_dict.keys())
        starts = [sorted_dict[k].start for k in ids]
        assert starts == sorted(starts)

    def test_sort_single_gene(self):
        g = Gene(False, False, "gX", "chr1", "aegis", "gene", "+", 100, 500, ".")
        ch, sorted_dict = sort_and_update_genes("chr1", {"gX": g})
        assert ch == "chr1"
        assert list(sorted_dict.keys()) == ["gX"]

    def test_sort_already_sorted(self):
        g1 = Gene(False, False, "g1", "chr1", "aegis", "gene", "+", 100, 300, ".")
        g2 = Gene(False, False, "g2", "chr1", "aegis", "gene", "+", 400, 600, ".")
        ch, sorted_dict = sort_and_update_genes("chr1", {"g1": g1, "g2": g2})
        assert list(sorted_dict.keys()) == ["g1", "g2"]


# ============================================================
# Annotation — integration with small GFF3
# ============================================================

class TestAnnotationSmallGFF3:
    def test_load_minimal_gff3(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        assert annot.name == "minimal"
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
        assert "gene1" in t.parents # type: ignore

    def test_multi_gene_gff3(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 2
        assert "geneA" in annot.all_gene_ids
        assert "geneB" in annot.all_gene_ids
        assert "chr1" in annot.chrs
        assert "chr2" in annot.chrs

    def test_create_gene_id(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 2
        assert "geneA" in annot.all_gene_ids
        assert "geneB" in annot.all_gene_ids
        assert "chr1" in annot.chrs
        assert "chr2" in annot.chrs
        assert annot.chrs["chr2"]["geneB"].gene_id == None

        annot.create_featurecounts_ids()

        assert annot.chrs["chr2"]["geneB"].gene_id == "geneB"
        assert annot.chrs["chr2"]["geneB"].transcripts["mRNA_B"].gene_id == "geneB"
        assert annot.chrs["chr2"]["geneB"].transcripts["mRNA_B"].exons[0].gene_id == "geneB"

    def test_print_gff(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        output = annot.chrs["chr1"]["geneA"].print_gff()
        assert "geneA" in output
        assert "TF1" not in output
        output = annot.chrs["chr1"]["geneA"].print_gff(aliases=True, names=True, print_empty_attributes=True)
        assert "TF1" in output
        assert "TF2" in output

    def test_copy(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot2 = annot.copy()
        annot2.name = "changed"
        assert annot.name == "minimal"

    def test_feature_counts_include_subfeatures(self, sample_gff3_file):
        """Verify features dict counts exons, CDS, UTR types."""
        annot = Annotation(sample_gff3_file, quiet=True)
        assert "exon" in annot.features
        assert "CDS" in annot.features

        assert annot.features["exon"] >= 1
        assert annot.features["CDS"] >= 1

    def test_transcript_has_exons(self, sample_gff3_file):
        """Verify transcript objects have exon sub-features."""
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        t = list(gene.transcripts.values())[0]
        assert len(t.exons) >= 2

    def test_transcript_has_cds(self, sample_gff3_file):
        """Verify transcript objects have CDS sub-features."""
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        t = list(gene.transcripts.values())[0]
        assert len(t.CDSs) >= 1


# ============================================================
# Annotation — rich GFF3 (3 genes, alt transcripts, TE, NC)
# ============================================================

class TestAnnotationRichGFF3:
    def test_load_three_genes(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 3
        assert "geneR1" in annot.all_gene_ids
        assert "geneR2" in annot.all_gene_ids
        assert "geneR3" in annot.all_gene_ids

    def test_two_chromosomes(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        assert "chr1" in annot.chrs
        assert "chr2" in annot.chrs

    def test_alternative_transcripts(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["geneR1"]
        assert len(gene.transcripts) == 2

    def test_noncoding_transcript_type(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        gene = annot.chrs["chr2"]["geneR3"]
        t = list(gene.transcripts.values())[0]
        assert t.feature == "lnc_RNA"


# ============================================================
# Annotation — slow integration test with real test_data
# ============================================================

class TestAnnotationRealData:
    def test_load_grapevine(self, test_data_dir):
        gff_path = str(test_data_dir / "input/annotation/grapevine_v5.1.gff3")
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
        gtf_file = str(TEST_DATA_DIR / "input/annotation/convert_basic.gtf")
        gff_file = tmp_path / "test.gff3"

        convert_gtf_to_gff3(gtf_file, str(gff_file), "utf-8", quiet=True)

        gff_content = gff_file.read_text()
        # Verify headers and format
        assert "##gff-version 3" in gff_content
        assert "ID=GENE1" in gff_content
        assert "Parent=GENE1" in gff_content
        assert "ID=T1" in gff_content

    def test_convert_gtf_with_cds_and_exon(self, tmp_path):
        """convert_gtf_to_gff3 emits gene, transcript, AND subfeature lines"""
        gtf_file = str(TEST_DATA_DIR / "input/annotation/convert_cds.gtf")
        gff_file = tmp_path / "cds.gff3"
        convert_gtf_to_gff3(gtf_file, str(gff_file), "utf-8", quiet=True)
        gff_content = gff_file.read_text()
        assert "ID=G1" in gff_content
        assert "ID=T1" in gff_content
        # Subfeature lines should now be present
        lines = [l for l in gff_content.strip().split("\n") if not l.startswith("#")]
        features = [l.split("\t")[2] for l in lines]
        assert "exon" in features
        assert "CDS" in features

    def test_convert_gtf_exon_cds_only(self, tmp_path):
        """GTF with only exon/CDS rows should infer gene and transcript lines"""
        gtf_file = str(TEST_DATA_DIR / "input/annotation/convert_exon_cds_only.gtf")
        gff_file = tmp_path / "inferred.gff3"
        convert_gtf_to_gff3(gtf_file, str(gff_file), "utf-8", quiet=True)
        gff_content = gff_file.read_text()

        assert "##gff-version 3" in gff_content
        lines = [l for l in gff_content.strip().split("\n") if not l.startswith("#")]
        features = [l.split("\t")[2] for l in lines]

        assert "gene" in features
        assert "ID=g1" in gff_content

        assert "mRNA" in features
        assert "ID=t1" in gff_content
        assert "ID=t2" in gff_content

        assert "exon" in features
        assert "CDS" in features

    def test_convert_gtf_exon_cds_only_gene_boundaries(self, tmp_path):
        """Inferred gene boundaries should span all subfeatures"""
        gtf_file = str(TEST_DATA_DIR / "input/annotation/convert_exon_cds_only.gtf")
        gff_file = tmp_path / "bounds.gff3"
        convert_gtf_to_gff3(gtf_file, str(gff_file), "utf-8", quiet=True)
        lines = [l for l in gff_file.read_text().strip().split("\n") if not l.startswith("#")]
        gene_lines = [l for l in lines if l.split("\t")[2] == "gene"]
        assert len(gene_lines) == 1
        parts = gene_lines[0].split("\t")
        # Gene should span from earliest exon start to latest exon end
        assert int(parts[3]) == 100
        assert int(parts[4]) == 7400

    def test_annotation_from_exon_only_gtf(self, tmp_path):
        """Full integration: Annotation should load exon/CDS-only GTF correctly"""
        gtf_file = str(TEST_DATA_DIR / "input/annotation/convert_exon_cds_only.gtf")
        annot = Annotation(gtf_file, quiet=True)

        assert "g1" in annot.all_gene_ids

        assert "t1" in annot.all_transcript_ids
        assert "t2" in annot.all_transcript_ids
        
        gene = annot.chrs["chr1"]["g1"]
        assert len(gene.transcripts) == 2
        t_re = gene.transcripts["t1"]
        assert len(t_re.exons) == 5
        t_rf = gene.transcripts["t2"]
        assert len(t_rf.exons) == 3

    def test_convert_gtf_keyless_attributes(self, keyless_attributes_gtf_file, tmp_path):
        """GTF with bare/keyless attributes should convert to valid GFF3"""
        gff_file = tmp_path / "keyless.gff3"
        convert_gtf_to_gff3(keyless_attributes_gtf_file, str(gff_file), "utf-8", quiet=True)
        gff_content = gff_file.read_text()

        assert "##gff-version 3" in gff_content
        assert "ID=PITA_19277" in gff_content
        assert "Parent=PITA_19277" in gff_content
        assert "ID=PITA_36893" in gff_content
        assert "ID=PITA_43071" in gff_content


# ============================================================
# Testing Keyless Attributes GTF (e.g. GFACS output)
# ============================================================

class TestKeylessAttributesGtf:
    def test_load_keyless_attributes_gtf(self, keyless_attributes_gtf_file):
        """Full integration: Annotation should load GTF with keyless attributes directly"""
        annot = Annotation(keyless_attributes_gtf_file, quiet=True)

        expected_genes = {
            "PITA_19277", "PITA_36893", "PITA_43071",
            "PITA_31484", "PITA_22913", "PITA_43120", "PITA_19387"
        }
        assert set(annot.all_gene_ids.keys()) == expected_genes
        assert len(annot.chrs) == 7

        # Verify gene features and structure on super32
        gene_19277 = annot.chrs["super32"]["PITA_19277"]
        assert gene_19277.start == 267921
        assert gene_19277.end == 811217
        assert gene_19277.strand == "+"
        assert len(gene_19277.transcripts) == 1

        transcript = list(gene_19277.transcripts.values())[0]
        assert len(transcript.exons) == 9
        assert len(transcript.CDSs) == 1

    def test_keyless_attributes_roundtrip_export(self, keyless_attributes_gtf_file, tmp_path):
        """Verify exporting to GFF3 and reloading preserves all genes and structure"""
        annot = Annotation(keyless_attributes_gtf_file, quiet=True)
        annot.export.gff(output_dir=tmp_path, filename="keyless_exported.gff3", subfolder=False, quiet=True)

        reloaded = Annotation(str(tmp_path / "keyless_exported.gff3"), quiet=True)
        assert set(reloaded.all_gene_ids.keys()) == set(annot.all_gene_ids.keys())
        assert set(reloaded.chrs.keys()) == set(annot.chrs.keys())

    def test_keyless_attributes_subset(self, keyless_attributes_gtf_file):
        """Verify subsetting works as expected on keyless attributes annotations"""
        annot = Annotation(keyless_attributes_gtf_file, quiet=True)
        annot.subset(chosen_features=["super32"], no_gene_cap=True, quiet=True)
        assert len(annot.chrs) == 1
        assert "super32" in annot.chrs
        assert "PITA_19277" in annot.all_gene_ids


# ============================================================
# Annotation unique IDs
# ============================================================

class TestAnnotationUniqueIDs:
    def test_get_unique_gene_id(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        # Assuming sample has "gene1"
        annot.all_gene_ids = {"gene1": "chr1", "gene1_1": "chr2"}
        new_id = annot._get_unique_gene_id("gene1")
        assert new_id == "gene1_2"

        new_id2 = annot._get_unique_gene_id("gene2")
        assert new_id2 == "gene2"

    def test_get_unique_transcript_id(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.all_transcript_ids = {"t1": ("chr1", "gene1"), "t1_1": ("chr2", "gene1_1")}
        new_id = annot._get_unique_transcript_id("t1")
        assert new_id == "t1_2"

    def test_get_unique_gene_id_no_collision(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.all_gene_ids = {"geneX": "chr1"}
        # No collision for "geneY"
        assert annot._get_unique_gene_id("geneY") == "geneY"

    def test_get_unique_transcript_id_no_collision(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.all_transcript_ids = {"tX": ("chr1", "gene1")}
        assert annot._get_unique_transcript_id("tY") == "tY"


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
        rrna_t = Transcript("rRNA1", "chr1", "aegis", "rRNA", "+", 10, 50, ".", ["gene1"])
        gene1.transcripts["rRNA1"] = rrna_t

        assert len(gene1.transcripts) == 2
        annot.remove_other_mRNA_transcripts_from_rRNA_genes()
        assert "mRNA1" not in gene1.transcripts
        assert "rRNA1" in gene1.transcripts

    def test_mark_te_file_with_missing_genes(self, multi_gene_gff3_file, tmp_path):
        """TE file with genes not in annotation should not crash."""
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        te_file = tmp_path / "te_genes.txt"
        te_file.write_text("geneA\nnonexistent_gene\n")
        annot.mark_transposable_element_genes(str(te_file))
        assert annot.chrs["chr1"]["geneA"].transposable is True


# ============================================================
# Annotation — update, update_features
# ============================================================

class TestAnnotationUpdate:
    def test_update_refreshes_features(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        old_features = dict(annot.features)
        annot.update(quiet=True)
        # features dict should be non-empty and consistent after update
        assert "gene" in annot.features
        assert "mRNA" in annot.features
        assert annot.features["gene"] == old_features["gene"]

    def test_update_features_counts(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        annot.update_features(standardise=False, quiet=True)
        # Should count multiple feature types
        assert annot.features["gene"] == 3
        assert "exon" in annot.features


# ============================================================
# Annotation — sort_genes
# ============================================================

class TestAnnotationSortGenes:
    def test_sort_sets_sorted_flag(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.sorted = False
        annot.sort_genes(processes=1, quiet=True)
        assert annot.sorted is True

    def test_sort_genes_order(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        annot.sort_genes(processes=1, quiet=True)
        # chr1 has geneR1 (start=1000) and geneR2 (start=6000)
        gene_ids = list(annot.chrs["chr1"].keys())
        starts = [annot.chrs["chr1"][g].start for g in gene_ids]
        assert starts == sorted(starts)


# ============================================================
# Annotation — gene_count
# ============================================================

class TestAnnotationGeneCount:
    def test_gene_count_prints(self, multi_gene_gff3_file, capsys):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.stats.gene_count()
        captured = capsys.readouterr()
        assert "2 gene objects" in captured.out
        assert "2 genes in all gene ids" in captured.out


# ============================================================
# Annotation — return_random_gene_ids
# ============================================================

class TestAnnotationReturnRandomGeneIds:
    def test_returns_correct_number(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        ids = annot.return_random_gene_ids(number=1, coding=False)
        assert len(ids) == 1
        assert ids[0] in annot.all_gene_ids

    def test_avoids_specified_ids(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        # Both genes should be in all_gene_ids; avoid one
        ids = annot.return_random_gene_ids(number=1, to_avoid=["geneA"], coding=False)
        assert "geneA" not in ids

    def test_returns_unique_ids(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        ids = annot.return_random_gene_ids(number=2, coding=False)
        assert len(ids) == 2
        assert len(set(ids)) == 2


# ============================================================
# Annotation — remove_genes
# ============================================================

class TestAnnotationRemoveGenes:
    def test_remove_single_gene(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        assert "geneA" in annot.all_gene_ids
        annot.remove_genes(to_remove={"geneA"}, quiet=True)
        assert "geneA" not in annot.all_gene_ids
        assert "geneA" not in annot.chrs.get("chr1", {})

    def test_remove_nonexistent_gene_warns(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            annot.remove_genes(to_remove={"fake_gene"}, quiet=True)
            assert any("fake_gene" in str(warning.message) for warning in w)

    def test_remove_all_genes(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.remove_genes(to_remove={"geneA", "geneB"}, quiet=True)
        assert len(annot.all_gene_ids) == 0

    def test_update_gene_and_transcript_list(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        # Manual manipulation then update
        annot.all_gene_ids = {}
        annot.all_transcript_ids = {}
        annot.update_gene_and_transcript_list(quiet=True)
        assert len(annot.all_gene_ids) == 2
        assert len(annot.all_transcript_ids) == 2


# ============================================================
# Annotation — remove_transcripts
# ============================================================

class TestAnnotationRemoveTranscripts:
    def test_remove_transcript_by_id(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        assert "mRNA_A" in annot.all_transcript_ids
        annot.remove_transcripts(to_remove={"mRNA_A"}, remove_genes_accordingly=True, quiet=True)
        # Removing the only transcript should also remove the parent gene
        assert "geneA" not in annot.all_gene_ids

    def test_remove_nonexistent_transcript_warns(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            annot.remove_transcripts(to_remove={"fake_transcript"}, remove_genes_accordingly=True, quiet=False)
            assert any("fake_transcript" in str(warning.message) for warning in w)


# ============================================================
# Annotation — remove_transcripts_with_no_exons
# ============================================================

class TestAnnotationRemoveTranscriptsWithNoExons:
    def test_removes_exon_less_transcript(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        # Add a bogus transcript with no exons
        empty_t = Transcript("emptyT", "chr1", "aegis", "mRNA", "+", 100, 200, ".", ["gene1"])
        empty_t.exons = []
        gene.transcripts["emptyT"] = empty_t
        original_t_count = len(gene.transcripts)
        assert original_t_count >= 2

        annot.detect_transcripts_with_no_exons(remove_transcripts=True, remove_genes_accordingly=True, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        assert "emptyT" not in gene.transcripts
        assert "transcript_with_no_exons" in annot.warnings
        assert "emptyT" in annot.warnings["transcript_with_no_exons"]


# ============================================================
# Annotation — remove_exons_with_unmatched_strand
# ============================================================

class TestAnnotationRemoveExonsWithUnmatchedStrand:
    def test_removes_mismatched_strand_exon(self, sample_gff3_file):
        from aegis.subfeatures import Exon
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        t = list(gene.transcripts.values())[0]
        assert t.strand == "+"

        # Add an exon on the wrong strand
        wrong_exon = Exon("wrongE", "chr1", "aegis", "exon", "-", 100, 200, ".", ["mRNA1"])
        t.exons.append(wrong_exon)
        original_count = len(t.exons)

        annot.remove_exons_with_unmatched_strand(quiet=True)
        t_after = list(annot.chrs["chr1"]["gene1"].transcripts.values())[0]
        assert len(t_after.exons) < original_count
        for e in t_after.exons:
            assert e.strand == t_after.strand


# ============================================================
# Annotation — remove_chromosomes
# ============================================================

class TestAnnotationRemoveChromosomes:
    def test_remove_chromosome(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        assert "chr2" in annot.chrs
        annot.remove_chromosomes(features_to_remove={"chr2"}, update=True, quiet=True)
        assert "chr2" not in annot.chrs
        assert "geneB" not in annot.all_gene_ids

    def test_remove_empty_set(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        orig_chrs = set(annot.chrs.keys())
        annot.remove_chromosomes(features_to_remove=set(), update=False, quiet=True)
        assert set(annot.chrs.keys()) == orig_chrs


# ============================================================
# Annotation — filter_by_rna_class
# ============================================================

class TestAnnotationFilterByRnaClass:
    def test_keep_only_mrna(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        # geneR1 has mRNA + lnc_RNA; geneR3 has lnc_RNA only
        annot.filter_by_rna_class(rna_classes=["mRNA"], remove_genes_accordingly=True, quiet=True)
        # geneR3 (lnc_RNA only) should have been removed
        assert "geneR3" not in annot.all_gene_ids
        # geneR1 should keep only the mRNA transcript
        gene = annot.chrs["chr1"]["geneR1"]
        for t in gene.transcripts.values():
            assert t.feature == "mRNA"

    def test_keep_lnc_rna(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        annot.filter_by_rna_class(rna_classes=["lnc_RNA"], remove_genes_accordingly=True, quiet=True)
        # geneR2 (mRNA only) should have been removed
        assert "geneR2" not in annot.all_gene_ids
        # geneR3 should still exist
        assert "geneR3" in annot.all_gene_ids


# ============================================================
# Annotation — correct_gene_transcript_and_subfeature_coordinates
# ============================================================

class TestAnnotationCorrectCoordinates:
    def test_fix_transcript_extending_beyond_gene(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        t = list(gene.transcripts.values())[0]

        # Make gene coordinates narrower than the transcript
        gene.start = 1500
        gene.end = 4000
        annot.correct_gene_transcript_and_subfeature_coordinates(quiet=True)
        # Gene should now match transcript boundaries
        assert gene.start <= t.start
        assert gene.end >= t.end

    def test_fix_gene_too_long(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        gene = annot.chrs["chr1"]["gene1"]
        # Make gene coordinates wider than all transcripts
        gene.start = 1
        gene.end = 50000
        annot.correct_gene_transcript_and_subfeature_coordinates(quiet=True)
        # Gene should have been trimmed to match transcript extent
        t = list(gene.transcripts.values())[0]
        assert gene.start == t.start
        assert gene.end == t.end


# ============================================================
# Annotation — gene_list
# ============================================================

class TestAnnotationGeneList:
    def test_gene_list_basic(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.gene_list(output_dir=str(tmp_path))
        output_file = tmp_path / "multi_gene_genes.txt"
        assert output_file.exists()
        content = output_file.read_text()
        assert "gene_id" in content  # header
        assert "geneA" in content
        assert "geneB" in content

    def test_gene_list_with_coordinates(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.gene_list(output_dir=str(tmp_path), coordinates=True)
        output_file = tmp_path / "multi_gene_genes.txt"
        content = output_file.read_text()
        assert "chromosome" in content
        assert "gene_start" in content
        assert "gene_end" in content

    def test_gene_list_with_lengths(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.gene_list(output_dir=str(tmp_path), lengths=True)
        output_file = tmp_path / "multi_gene_genes.txt"
        content = output_file.read_text()
        assert "gene_length" in content

    def test_gene_list_custom_output_file(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.gene_list(output_dir=str(tmp_path), filename="custom_list.txt")
        output_file = tmp_path / "custom_list.txt"
        assert output_file.exists()

    def test_gene_list_skip_both_coding_and_noncoding(self, multi_gene_gff3_file, tmp_path, capsys):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.gene_list(output_dir=str(tmp_path), skip_coding=True, skip_non_coding=True)
        captured = capsys.readouterr()
        assert "Warning" in captured.out


# ============================================================
# Annotation — transcript_list
# ============================================================

class TestAnnotationTranscriptList:
    def test_transcript_list_basic(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.transcript_list(output_dir=str(tmp_path))
        output_file = tmp_path / "multi_gene_transcripts.txt"
        assert output_file.exists()
        content = output_file.read_text()
        assert "transcript_id" in content
        assert "mRNA_A" in content
        assert "mRNA_B" in content

    def test_transcript_list_with_coordinates(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.transcript_list(output_dir=str(tmp_path), coordinates=True)
        output_file = tmp_path / "multi_gene_transcripts.txt"
        content = output_file.read_text()
        assert "transcript_start" in content
        assert "transcript_end" in content

    def test_transcript_list_skip_both(self, multi_gene_gff3_file, tmp_path, capsys):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.export.transcript_list(output_dir=str(tmp_path), skip_coding=True, skip_non_coding=True)
        captured = capsys.readouterr()
        assert "Warning" in captured.out


# ============================================================
# Annotation — export_gff
# ============================================================

class TestAnnotationExportGff:
    def test_export_gff3(self, sample_gff3_file, tmp_path):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.export.gff(output_dir=str(tmp_path), subfolder=True, quiet=True)
        out_dir = tmp_path / "out_gffs"
        assert out_dir.exists()
        gff_files = list(out_dir.glob("*.gff3"))
        assert len(gff_files) >= 1
        content = gff_files[0].read_text()
        assert "##gff-version 3" in content
        assert "gene1" in content

    def test_export_gff3_just_genes(self, sample_gff3_file, tmp_path):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.export.gff(output_dir=str(tmp_path), just_genes=True, subfolder=True, quiet=True)
        out_dir = tmp_path / "out_gffs"
        gff_files = list(out_dir.glob("*.gff3"))
        content = gff_files[0].read_text()
        assert "gene1" in content
        # Should not have transcript-level features when just_genes=True
        lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 3:
                assert parts[2] == "gene"

    def test_export_gff3_no_subfolder(self, sample_gff3_file, tmp_path):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.export.gff(output_dir=str(tmp_path), subfolder=False, quiet=True)
        gff_files = list(tmp_path.glob("*.gff3"))
        assert len(gff_files) >= 1


# ============================================================
# Annotation — export_gtf
# ============================================================

class TestAnnotationExportGtf:
    def test_export_gtf(self, sample_gff3_file, tmp_path):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.export.gtf(output_dir=str(tmp_path), subfolder=True, quiet=True)
        out_dir = tmp_path / "out_gtfs"
        assert out_dir.exists()
        gtf_files = list(out_dir.glob("*.gtf"))
        assert len(gtf_files) >= 1
        content = gtf_files[0].read_text()
        assert "#gtf-version 2.2" in content
        assert "gene1" in content

    def test_export_gtf_just_genes(self, sample_gff3_file, tmp_path):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.export.gtf(output_dir=str(tmp_path), just_genes=True, subfolder=True, quiet=True)
        out_dir = tmp_path / "out_gtfs"
        gtf_files = list(out_dir.glob("*.gtf"))
        content = gtf_files[0].read_text()
        lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 3:
                assert parts[2] == "gene"


# ============================================================
# Annotation — rename_chromosomes
# ============================================================

class TestAnnotationRenameChromosomes:
    def test_rename_chromosome(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        equivalences = {"chr1": "scaffold_1", "chr2": "scaffold_2"}
        annot.rename_chromosomes(equivalences)
        assert "scaffold_1" in annot.chrs
        assert "scaffold_2" in annot.chrs
        assert "chr1" not in annot.chrs
        assert "chr2" not in annot.chrs
        # Verify nested gene chromosome updated
        gene = annot.chrs["scaffold_1"]["geneA"]
        assert gene.ch == "scaffold_1"

    def test_rename_updates_transcript_chromosome(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        equivalences = {"chr1": "sc1"}
        annot.rename_chromosomes(equivalences)
        gene = annot.chrs["sc1"]["geneA"]
        t = list(gene.transcripts.values())[0]
        assert t.ch == "sc1"

    def test_rename_nonexistent_chrom_is_noop(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        equivalences = {"chrX": "scaffoldX"}
        annot.rename_chromosomes(equivalences)
        # Original chromosomes should be untouched
        assert "chr1" in annot.chrs
        assert "chr2" in annot.chrs

    def test_rename_sets_confrenamed(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        equivalences = {"chr1": "scaffold_1"}
        annot.rename_chromosomes(equivalences)
        assert "confrenamed" in annot.tags

# ============================================================
# Annotation — clear_gene_names_and_symbols
# ============================================================

class TestAnnotationClearGeneNamesAndSymbols:
    def test_clear_names_and_symbols(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        # Pre-populate
        gene = annot.chrs["chr1"]["geneA"]
        gene.names = ["MyName"]
        gene.symbols = ["SYM"]
        gene.synonyms = ["Syn1"]

        annot.clear_gene_names_and_symbols(quiet=True)
        assert gene.names is None
        assert gene.symbols is None
        assert gene.synonyms is None

        assert gene.aliases == ["TF1", "TF2"]

        annot.clear_aliases()

        assert gene.aliases is None        


# ============================================================
# Annotation — remove_TE_genes
# ============================================================

class TestAnnotationRemoveTEGenes:
    def test_remove_te_genes(self, multi_gene_gff3_file, tmp_path):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        # Mark one gene as TE
        annot.chrs["chr1"]["geneA"].transposable = True
        annot.remove_TE_genes(quiet=True)
        assert "geneA" not in annot.all_gene_ids
        assert "geneB" in annot.all_gene_ids
        assert "minus_TE" in annot.feature_tags

    def test_remove_te_genes_none_marked(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.remove_TE_genes(quiet=True)
        # Nothing marked, so nothing removed
        assert len(annot.all_gene_ids) == 2
        assert "minus_TE" not in annot.feature_tags


# ============================================================
# Annotation — remove_non_coding_genes_and_transcripts
# ============================================================

class TestAnnotationRemoveNonCodingGenes:
    def test_remove_noncoding_gene(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        # geneR3 is noncoding (lnc_RNA only, no CDS)
        annot.remove_non_coding_genes_and_transcripts(quiet=True)
        assert "geneR3" not in annot.all_gene_ids
        # geneR1 and geneR2 should remain (they have coding transcripts)
        assert "geneR1" in annot.all_gene_ids
        assert "geneR2" in annot.all_gene_ids


# ============================================================
# Annotation — remove_coding_genes_and_transcripts
# ============================================================

class TestAnnotationRemoveCodingGenes:
    def test_remove_coding_gene(self, rich_gff3_file):
        annot = Annotation(rich_gff3_file, quiet=True)
        annot.remove_coding_genes_and_transcripts(quiet=True)
        # Pure coding genes should be removed
        assert "geneR2" not in annot.all_gene_ids
        # geneR3 (noncoding) should remain
        assert "geneR3" in annot.all_gene_ids


# ============================================================
# Annotation — clear_sequences
# ============================================================

class TestAnnotationClearSequences:
    def test_contains_no_proteins_initially(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)

        assert annot.contains_protein_sequences is False

    def test_clear_protein_flag(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.contains_protein_sequences = True
        annot.clear_proteins()

        assert annot.contains_protein_sequences is False

# ============================================================
# Annotation — copy (deeper tests)
# ============================================================

class TestAnnotationCopyDeep:
    def test_copy_independence(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot2 = annot.copy()
        # Mutation in copy must not affect original
        annot2.chrs["chr1"]["gene1"].start = 9999
        assert annot.chrs["chr1"]["gene1"].start != 9999

    def test_copy_preserves_gene_ids(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot2 = annot.copy()
        assert set(annot2.all_gene_ids.keys()) == set(annot.all_gene_ids.keys())

    def test_copy_preserves_chromosomes(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot2 = annot.copy()
        assert set(annot2.chrs.keys()) == set(annot.chrs.keys())


# ============================================================
# Annotation — remove_genes_with_no_transcripts
# ============================================================

class TestAnnotationRemoveGenesWithNoTranscripts:
    def test_removes_empty_gene(self, sample_gff3_file):
        from aegis.gene import Gene
        annot = Annotation(sample_gff3_file, quiet=True)
        # Add a gene with no transcripts
        empty_gene = Gene(False, False, "emptyG", "chr1", "aegis", "gene", "+", 100, 200, ".")
        annot.chrs["chr1"]["emptyG"] = empty_gene
        annot.all_gene_ids["emptyG"] = "chr1"

        annot.detect_genes_with_no_transcripts(remove=True, quiet=True)
        assert "emptyG" not in annot.chrs["chr1"]
        assert "gene_with_no_transcripts" in annot.warnings
        assert "emptyG" in annot.warnings["gene_with_no_transcripts"]

    def test_keeps_gene_with_transcripts(self, sample_gff3_file):
        annot = Annotation(sample_gff3_file, quiet=True)
        annot.detect_genes_with_no_transcripts(remove=True, quiet=True)
        # gene1 has transcripts, so it should remain
        assert "gene1" in annot.chrs["chr1"]


# ============================================================
# Annotation — remove_genes_without_symbols
# ============================================================

class TestAnnotationRemoveGenesWithoutSymbols:
    def test_removes_genes_with_empty_symbols(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        # All genes start without symbols, so all should be removed
        annot.remove_genes_without_symbols(quiet=True)
        assert len(annot.all_gene_ids) == 0

    def test_keeps_genes_with_symbols(self, multi_gene_gff3_file):
        annot = Annotation(multi_gene_gff3_file, quiet=True)
        annot.chrs["chr1"]["geneA"].symbols = ["SYM1"]
        annot.remove_genes_without_symbols(quiet=True)
        assert "geneA" in annot.all_gene_ids
        assert "geneB" not in annot.all_gene_ids


# ============================================================
# Annotation — rework_CDSs (generate CDS from exon-only GFF3)
# ============================================================

class TestAnnotationReworkCDSs:
    """Load genome and exon-only Arabidopsis mRNA GFF3, run rework_CDSs,
    and verify the generated CDS segments match expected values."""

    @pytest.fixture
    def reworked_annotation(self):
        from aegis.genome import Genome
        genome = Genome(
            "ara_genome",
            str(TEST_DATA_DIR / "input/fasta/arabidopsis_tair10.fasta"),
            quiet=True,
        )
        annot = Annotation(
            str(TEST_DATA_DIR / "input/annotation/arabidopsis_generate_cds_mRNA.gff3"),
            "ara_annotation",
            quiet=True,
            genome=genome,
        )
        annot.rework_CDSs(quiet=True)
        return annot

    def test_all_transcripts_have_cds(self, reworked_annotation):
        """After rework_CDSs every mRNA transcript should be coding."""
        annot = reworked_annotation
        gene = annot.chrs["Chr4"]["AT4G00050"]
        for t in gene.transcripts.values():
            assert len(t.CDSs) == 1

    def test_generated_cds_segments(self, reworked_annotation):
        """Verify every CDS segment matches the expected coordinates,
        phase, ID and parent — grouped by transcript."""
        annot = reworked_annotation
        gene = annot.chrs["Chr4"]["AT4G00050"]

        # Expected CDS segments keyed by transcript ID
        expected = {
            "AT4G00050.1": [
                ("Chr4", "Araport11", 17863, 17954, "+", 0, "AT4G00050.1_CDS001", "AT4G00050.1"),
                ("Chr4", "Araport11", 18030, 18513, "+", 1, "AT4G00050.1_CDS001", "AT4G00050.1"),
                ("Chr4", "Araport11", 18600, 18692, "+", 0, "AT4G00050.1_CDS001", "AT4G00050.1"),
                ("Chr4", "Araport11", 18805, 18870, "+", 0, "AT4G00050.1_CDS001", "AT4G00050.1"),
                ("Chr4", "Araport11", 19296, 19673, "+", 0, "AT4G00050.1_CDS001", "AT4G00050.1"),
                ("Chr4", "Araport11", 19762, 19848, "+", 0, "AT4G00050.1_CDS001", "AT4G00050.1"),
            ],
            "AT4G00050.3": [
                ("Chr4", "Araport11", 17863, 17954, "+", 0, "AT4G00050.3_CDS001", "AT4G00050.3"),
                ("Chr4", "Araport11", 18030, 18513, "+", 1, "AT4G00050.3_CDS001", "AT4G00050.3"),
                ("Chr4", "Araport11", 18600, 18692, "+", 0, "AT4G00050.3_CDS001", "AT4G00050.3"),
                ("Chr4", "Araport11", 18805, 18870, "+", 0, "AT4G00050.3_CDS001", "AT4G00050.3"),
                ("Chr4", "Araport11", 19296, 19715, "+", 0, "AT4G00050.3_CDS001", "AT4G00050.3"),
            ],
            "AT4G00050.2": [
                ("Chr4", "Araport11", 18244, 18513, "+", 0, "AT4G00050.2_CDS001", "AT4G00050.2"),
                ("Chr4", "Araport11", 18600, 18692, "+", 0, "AT4G00050.2_CDS001", "AT4G00050.2"),
                ("Chr4", "Araport11", 18805, 18870, "+", 0, "AT4G00050.2_CDS001", "AT4G00050.2"),
                ("Chr4", "Araport11", 19296, 19673, "+", 0, "AT4G00050.2_CDS001", "AT4G00050.2"),
                ("Chr4", "Araport11", 19762, 19848, "+", 0, "AT4G00050.2_CDS001", "AT4G00050.2"),
            ],
        }

        for t_id, expected_segments in expected.items():

            t = gene.transcripts[t_id]

            # Collect all CDS segments across all CDS objects in the transcript
            actual_segments = []
            for cds in t.CDSs.values():
                for seg in cds.CDS_segments:
                    actual_segments.append((
                        seg.ch,
                        seg.source,
                        seg.start,
                        seg.end,
                        seg.strand,
                        seg.phase,
                        seg.id,
                        seg.parents[0],
                    ))

            print(actual_segments)

            print(expected_segments)

            assert len(actual_segments) == len(expected_segments)

            for i, (act, exp) in enumerate(zip(actual_segments, expected_segments)):
                assert act == exp


# ============================================================
# Edge-case GFF3 tests
# ============================================================

# ---- 1. Exon-only (no CDS/UTR) ----

class TestExonOnlyGFF3:
    """exons with no gene or mRNA"""

    def test_gene_and_transcript_loaded(self, exon_only_gff3_file):
        annot = Annotation(exon_only_gff3_file, quiet=True)
        assert "mRNA_eo1_gene" in annot.all_gene_ids
        assert "mRNA_eo1" in annot.all_transcript_ids

    def test_gene_coordinates(self, exon_only_gff3_file):
        annot = Annotation(exon_only_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["mRNA_eo1_gene"]
        coords = (g.start, g.end)
        assert coords == (1000, 5000)


# ---- 2. CDS-only ----

class TestCDSOnlyGFF3:
    """Gene→mRNA→CDS with no explicit exons"""

    def test_transcript_is_coding(self, cds_only_gff3_file):
        annot = Annotation(cds_only_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_co1"].transcripts["mRNA_co1"]
        assert t.coding is True

    def test_one_cds_with_three_segments(self, cds_only_gff3_file):
        annot = Annotation(cds_only_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_co1"].transcripts["mRNA_co1"]
        assert len(t.CDSs) == 1
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 3

    def test_exons_auto_generated_from_cds(self, cds_only_gff3_file):
        """When no exons are provided, aegis reconstructs them from CDS segments"""
        annot = Annotation(cds_only_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_co1"].transcripts["mRNA_co1"]
        assert len(t.exons) == 3

    def test_auto_exon_coordinates_match_cds(self, cds_only_gff3_file):
        annot = Annotation(cds_only_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_co1"].transcripts["mRNA_co1"]
        coords = [(e.start, e.end) for e in t.exons]
        assert coords == [(1200, 2000), (3000, 3800), (4200, 4800)]

    def test_two_introns_generated(self, cds_only_gff3_file):
        annot = Annotation(cds_only_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_co1"].transcripts["mRNA_co1"]
        annot.generate_introns()
        assert len(t.introns) == 2 # type: ignore

    def test_gene_coordinates_corrected(self, cds_only_gff3_file):
        """Gene coordinates are corrected to match the actual subfeature span"""
        annot = Annotation(cds_only_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_co1"]
        assert g.start == 1200
        assert g.end == 4800


# ---- 3. No subfeatures (single-exon fallback) ----

class TestNoSubfeaturesGFF3:
    """Gene with mRNA but zero subfeatures — single exon generated"""

    def test_two_genes_loaded(self, no_subfeatures_gff3_file):
        annot = Annotation(no_subfeatures_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 2

    def test_single_exon_spanning_transcript(self, no_subfeatures_gff3_file):
        """When no subfeatures exist, a single exon spanning the full transcript is created"""
        annot = Annotation(no_subfeatures_gff3_file, quiet=True)
        t1 = annot.chrs["chr1"]["gene_ns1"].transcripts["mRNA_ns1"]
        assert len(t1.exons) == 1
        assert t1.exons[0].start == 500
        assert t1.exons[0].end == 3000

    def test_not_coding(self, no_subfeatures_gff3_file):
        annot = Annotation(no_subfeatures_gff3_file, quiet=True)
        t1 = annot.chrs["chr1"]["gene_ns1"].transcripts["mRNA_ns1"]
        t2 = annot.chrs["chr2"]["gene_ns2"].transcripts["mRNA_ns2"]
        assert t1.coding is False
        assert t2.coding is False

    def test_no_introns(self, no_subfeatures_gff3_file):
        annot = Annotation(no_subfeatures_gff3_file, quiet=True)
        t1 = annot.chrs["chr1"]["gene_ns1"].transcripts["mRNA_ns1"]
        annot.generate_introns()
        assert len(t1.introns) == 0 # type: ignore

    def test_minus_strand_transcript(self, no_subfeatures_gff3_file):
        annot = Annotation(no_subfeatures_gff3_file, quiet=True)
        t2 = annot.chrs["chr2"]["gene_ns2"].transcripts["mRNA_ns2"]
        assert t2.strand == "-"
        assert len(t2.exons) == 1
        assert t2.exons[0].start == 100
        assert t2.exons[0].end == 2500


# ---- 4. Non-coding transcript types ----

class TestNoncodingTranscriptsGFF3:
    """lnc_RNA, tRNA, rRNA, snoRNA transcript types from conf.py"""

    def test_four_genes_loaded(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 4

    def test_all_transcripts_noncoding(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        for genes in annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    assert t.coding is False

    def test_transcript_feature_types_preserved(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        features = set()
        for genes in annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    features.add(t.feature)
        assert features == {"lnc_RNA", "tRNA", "rRNA", "snoRNA"}

    def test_lnc_rna_has_two_exons(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_nc1"].transcripts["lncRNA_nc1"]
        assert len(t.exons) == 2

    def test_trna_has_one_exon(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_nc2"].transcripts["tRNA_nc2"]
        assert len(t.exons) == 1

    def test_rrna_minus_strand(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        t = annot.chrs["chr2"]["gene_nc3"].transcripts["rRNA_nc3"]
        assert t.strand == "-"
        assert len(t.exons) == 2

    def test_no_cds_on_any_transcript(self, noncoding_transcripts_gff3_file):
        annot = Annotation(noncoding_transcripts_gff3_file, quiet=True)
        for genes in annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    assert len(t.CDSs) == 0


# ---- 5. Pseudogene ----

class TestPseudogeneGFF3:
    """Pseudogene with exons parented directly to the gene (no transcript).
    Aegis auto-creates a pseudotranscript"""

    def test_gene_is_pseudogene(self, pseudogene_gff3_file):
        annot = Annotation(pseudogene_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_ps1"]
        assert g.pseudogene is True

    def test_pseudotranscript_auto_created(self, pseudogene_gff3_file):
        annot = Annotation(pseudogene_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_ps1"]
        assert len(g.transcripts) == 1
        assert "pseudo_t_ps1" in g.transcripts

    def test_pseudotranscript_has_three_exons(self, pseudogene_gff3_file):
        annot = Annotation(pseudogene_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ps1"].transcripts["pseudo_t_ps1"]
        assert len(t.exons) == 3

    def test_pseudotranscript_not_coding(self, pseudogene_gff3_file):
        annot = Annotation(pseudogene_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ps1"].transcripts["pseudo_t_ps1"]
        assert t.coding is False

    def test_pseudotranscript_exon_coordinates(self, pseudogene_gff3_file):
        annot = Annotation(pseudogene_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ps1"].transcripts["pseudo_t_ps1"]
        coords = [(e.start, e.end) for e in t.exons]
        assert coords == [(1000, 2000), (3000, 4000), (4500, 5000)]


# ---- 6. Overlapping, and adjacent exons (collapse) ----

class TestExonsToCollapseGFF3:
    """Gene with overlapping and adjacent exons is collapsed into fewer exons"""

    def test_exons_collapsed(self, exons_to_collapse_gff3_file):
        """4 overlapping/adjacent input exons -> 2 collapsed exons"""
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        assert len(t.exons) == 2

    def test_collapsed_exon_coordinates(self, exons_to_collapse_gff3_file):
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        coords = [(e.start, e.end) for e in t.exons]
        # exons 1000-2500 + 2000-3500 + 3500-4000 overlap → merged into 1000-4000
        # exon 5000-6000 remains separate
        assert coords == [(1000, 4000), (5000, 6000)]

    def test_transcript_is_coding(self, exons_to_collapse_gff3_file):
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        assert t.coding is True

    def test_one_intron_after_collapse(self, exons_to_collapse_gff3_file):
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        annot.generate_introns()
        assert len(t.introns) == 1 # type: ignore

    def test_cds_segments_collapsed(self, exons_to_collapse_gff3_file):
        """3 CDS input segments (2 overlapping) -> 2 collapsed CDS segments"""
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 2

    def test_collapsed_cds_coordinates(self, exons_to_collapse_gff3_file):
        """CDS 1200-2500 + 2000-3500 -> merged 1200-3500; CDS 5000-5800 stays"""
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        cds = list(t.CDSs.values())[0]
        coords = [(s.start, s.end) for s in cds.CDS_segments]
        assert coords == [(1200, 3500), (5000, 5800)]

    def test_no_collapse_cds_flag(self, exons_to_collapse_gff3_file):
        """With collapse_CDSs=False the 3 original CDS segments are kept"""
        annot = Annotation(exons_to_collapse_gff3_file, quiet=True, collapse_CDSs=False)
        t = annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"]
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 3

# ---- 7. Multi CDS IDs ----

class TestMultiCDSIdsGFF3:
    """CDS segments with different IDs"""

    def test_cds_segments_combined(self, multi_cds_ids_gff3_file):
        """CDS segments with different IDs but sorted into CDS objects"""
        annot = Annotation(multi_cds_ids_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_mc1"].transcripts["mRNA_mc1"]
        assert len(t.CDSs) == 1
        assert len(t.CDSs["CDS_mc1a"].CDS_segments) == 2

    def test_two_exons(self, multi_cds_ids_gff3_file):
        annot = Annotation(multi_cds_ids_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_mc1"].transcripts["mRNA_mc1"]
        assert len(t.exons) == 2


# ---- 8. Transcript without Parent (gene inferred) ----

class TestTranscriptNoParentGFF3:
    """mRNA with no Parent → gene auto-inferred"""

    def test_gene_auto_created(self, transcript_no_parent_gff3_file):
        annot = Annotation(transcript_no_parent_gff3_file, quiet=True)
        assert "mRNA_np1_gene" in annot.all_gene_ids

    def test_warning_raised(self, transcript_no_parent_gff3_file):
        annot = Annotation(transcript_no_parent_gff3_file, quiet=True)
        assert len(annot.warnings["transcript_with_no_parent"]) == 1


# ---- 9. Just CDS without Parent (gene + transcript inferred) ----

class TestCDSNoParentGFF3:
    """CDS subfeatures with no Parent → gene + transcript auto-inferred"""

    def test_gene_auto_created(self, cds_no_parent_gff3_file):
        annot = Annotation(cds_no_parent_gff3_file, quiet=True)
        assert "CDS_cnp1_gene" in annot.all_gene_ids

    def test_transcript_auto_created(self, cds_no_parent_gff3_file):
        annot = Annotation(cds_no_parent_gff3_file, quiet=True)
        assert "CDS_cnp1_transcript" in annot.all_transcript_ids

    def test_transcript_is_coding(self, cds_no_parent_gff3_file):
        annot = Annotation(cds_no_parent_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["CDS_cnp1_gene"].transcripts["CDS_cnp1_transcript"]
        assert t.coding is True

    def test_cds_has_two_segments(self, cds_no_parent_gff3_file):
        annot = Annotation(cds_no_parent_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["CDS_cnp1_gene"].transcripts["CDS_cnp1_transcript"]
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 2

    def test_exons_auto_generated(self, cds_no_parent_gff3_file):
        annot = Annotation(cds_no_parent_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["CDS_cnp1_gene"].transcripts["CDS_cnp1_transcript"]
        assert len(t.exons) == 2
        coords = [(e.start, e.end) for e in t.exons]
        assert coords == [(1200, 2000), (3000, 4800)]

    def test_warning_raised(self, cds_no_parent_gff3_file):
        annot = Annotation(cds_no_parent_gff3_file, quiet=True)
        assert len(annot.warnings["subfeature_with_no_parent"]) == 1


# ---- 10. Multiple isoforms ----

class TestMultipleIsoformsGFF3:
    """Gene with 3 mRNA isoforms sharing exon regions but different CDS"""

    def test_single_gene(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        assert len(annot.all_gene_ids) == 1
        assert "gene_iso1" in annot.all_gene_ids

    def test_three_transcripts(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_iso1"]
        assert len(g.transcripts) == 3

    def test_all_isoforms_coding(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_iso1"]
        for t in g.transcripts.values():
            assert t.coding is True

    def test_each_isoform_has_one_cds(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_iso1"]
        for t in g.transcripts.values():
            assert len(t.CDSs) == 1

    def test_isoform_a_has_three_exons(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"]
        assert len(t.exons) == 3

    def test_isoform_b_has_two_exons(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"]
        assert len(t.exons) == 2

    def test_isoform_c_has_two_exons(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        t = annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"]
        assert len(t.exons) == 2

    def test_transcript_ids_correct(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)
        g = annot.chrs["chr1"]["gene_iso1"]
        assert set(g.transcripts.keys()) == {"mRNA_iso1a", "mRNA_iso1b", "mRNA_iso1c"}


# ---- 11. Subfeatures reference gene as Parent (no transcript line) ----

class TestSubfeatureParentIsGene:
    """Exons and CDSs reference a gene ID as Parent with no mRNA/transcript
    feature in the GFF3. AEGIS should auto-create a transcript."""

    def test_gene_loaded(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        assert "g1" in annot.all_gene_ids

    def test_transcript_auto_created(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        g = annot.chrs["chr01"]["g1"]
        assert len(g.transcripts) == 1
        assert "g1_t1" in g.transcripts

    def test_transcript_is_coding(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        t = annot.chrs["chr01"]["g1"].transcripts["g1_t1"]
        assert t.coding is True

    def test_cds_present(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        t = annot.chrs["chr01"]["g1"].transcripts["g1_t1"]
        assert len(t.CDSs) == 1
        cds = list(t.CDSs.values())[0]
        assert len(cds.CDS_segments) == 3

    def test_exon_coordinates(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        t = annot.chrs["chr01"]["g1"].transcripts["g1_t1"]
        coords = [(e.start, e.end) for e in t.exons]
        assert coords == [(10000, 10400), (11000, 11300), (12100, 12500)]

    def test_transcript_boundaries(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        t = annot.chrs["chr01"]["g1"].transcripts["g1_t1"]
        assert t.start == 10000
        assert t.end == 12500

    def test_warning_raised(self, subfeature_parent_is_gene_gff3_file):
        annot = Annotation(subfeature_parent_is_gene_gff3_file, quiet=True)
        assert len(annot.warnings["subfeature_to_gene"]) > 0

# ---- 12. Non standard Parent attribute format ----

class TestNonStandardParentAttribute:
    """Tests proper reading of GFF3 files with no gene entries, geneID attributes and no ID on exons"""

    def test_four_genes_loaded(self, geneID_attribute_as_parent_gff3_file):
        annot = Annotation(geneID_attribute_as_parent_gff3_file, quiet=True)
        
        assert len(annot.all_gene_ids) == 4
        assert "g1" in annot.all_gene_ids
        assert "g2" in annot.all_gene_ids
        assert "g3" in annot.all_gene_ids
        assert "g4" in annot.all_gene_ids

    def test_transcripts_are_associated(self, geneID_attribute_as_parent_gff3_file):
        annot = Annotation(geneID_attribute_as_parent_gff3_file, quiet=True)
        
        g2 = annot.chrs["chr01"]["g2"]
        assert len(g2.transcripts) == 2
        assert "t2.1" in g2.transcripts
        assert "t2.2" in g2.transcripts
        
        g3 = annot.chrs["chr01"]["g3"]
        assert len(g3.transcripts) == 1
        assert "t3.1" in g3.transcripts

    def test_exons_are_parsed(self, geneID_attribute_as_parent_gff3_file):
        """Verify that exons are correctly assigned even without IDs in the GFF3"""
        annot = Annotation(geneID_attribute_as_parent_gff3_file, quiet=True)
        
        t = annot.chrs["chr01"]["g3"].transcripts["t3.1"]
        assert len(t.exons) == 4

        coords = sorted([(e.start, e.end) for e in t.exons])
        assert coords[0] == (3000, 3100)
        assert coords[-1] == (3900, 4000)

    def test_gene_boundaries(self, geneID_attribute_as_parent_gff3_file):
        """Gene boundaries should encompass all transcripts (min start and max end)"""
        annot = Annotation(geneID_attribute_as_parent_gff3_file, quiet=True)
        
        g = annot.chrs["chr01"]["g2"]
        assert g.start == 1000
        assert g.end == 2500

        g1 = annot.chrs["chr01"]["g1"]
        assert g1.start == 100
        assert g1.end == 500


# ============================================================
# Shared Exon Parents with multiple fw and rv transcripts also requiring collapse
# ============================================================

class TestAnnotationSharedExonParents:
    def test_shared_parents_are_detected_with_exons_to_collapse(self, shared_parents_gff3_file):
        annot = Annotation(shared_parents_gff3_file, quiet=True)

        assert len(annot.all_gene_ids) == 1

        assert len(annot.chrs["chr1"]["gene_ov1"].transcripts) == 5

        assert len(annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons) == 2

        assert len(annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov2"].exons) == 3

        assert len(annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA1"].exons) == 3

        assert len(annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA2"].exons) == 4
        
        assert len(annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA3"].exons) == 4

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov2"].exons[0].start == 500

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov2"].exons[0].end == 800

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov2"].exons[0].parents == ["mRNA_ov2"]

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons[0].start == 1000

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons[0].end == 1300

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov2"].exons[0].id == "ov1_e001"

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons[0].id == "ov1_e002"

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons[0].parents == ["mRNA_ov1", "mRNA_ov2"]

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA1"].exons[0].id == "ov1_e008"

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA1"].exons[0].parents == ["ncRNA1", "ncRNA2", "ncRNA3"]

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA1"].exons[-1].id == "ov1_e004"
        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA1"].exons[-1].parents == ["ncRNA1", "ncRNA2"]

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA3"].exons[-1].id == "ov1_e005"

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA3"].exons[-1].start == 4500

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["ncRNA3"].exons[-1].end == 6000

    def test_shared_parents_reformatting(self, shared_parents_gff3_file):
        annot = Annotation(shared_parents_gff3_file, quiet=True)

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons[0].parents == ["mRNA_ov1", "mRNA_ov2"]

        annot.single_parent_for_exons_utrs()

        assert annot.chrs["chr1"]["gene_ov1"].transcripts["mRNA_ov1"].exons[0].parents == ["mRNA_ov1"]


# ============================================================
# Clash of IDs and removal of genes with no transcripts, pseudogene or not, also some checks on shared parents of renamed transcripts
# ============================================================

class TestAnnotationClashOfIDs:
    annot: Annotation
    def test_clash_of_ids_and_transcriptless_gene_removal(self, clash_of_ids_gff3_file):
        annot = Annotation(clash_of_ids_gff3_file, quiet=True)

        assert len(annot.all_gene_ids) == 7

        assert "g1_1" not in annot.chrs["chr8"]
        assert "g1_1_1" in annot.chrs["chr8"]
        assert "g1" in annot.chrs["chr1"]
        assert "g1" not in annot.chrs["chr2"]
        assert "g1_1" in annot.chrs["chr2"]

        assert "g1.t1" in annot.chrs["chr1"]["g1"].transcripts
        assert "g1.t1_1" in annot.chrs["chr2"]["g1_1"].transcripts

        assert annot.chrs["chr2"]["g1_1"].transcripts["g1.t1_1"].exons[0].id == "exon_1"
        assert annot.chrs["chr2"]["g1_1"].transcripts["g1.t1_1"].exons[0].parents == ["g1.t1_1"]
        assert annot.chrs["chr2"]["g1_1"].transcripts["g1.t1_1"].exons[0].start == 500
        assert annot.chrs["chr2"]["g1_1"].transcripts["g1.t1_1"].exons[0].end == 2000

        assert "g1.t1" not in annot.chrs["chr2"]["g1_1"].transcripts

        assert "g1.t1" not in annot.chrs["chr3"]["g2"].transcripts
        assert "g1.t1_1" not in annot.chrs["chr3"]["g2"].transcripts
        assert "g1.t1_2" in annot.chrs["chr3"]["g2"].transcripts

        assert len(annot.chrs["chr3"]["g3"].transcripts) == 2
        assert "g3.t1" in annot.chrs["chr3"]["g3"].transcripts
        assert "g3.t1_1" in annot.chrs["chr3"]["g3"].transcripts

        assert len(annot.chrs["chr3"]["g3"].transcripts["g3.t1"].exons) == 2

        assert len(annot.chrs["chr3"]["g3"].transcripts["g3.t1_1"].exons) == 1

        assert annot.chrs["chr3"]["g3"].transcripts["g3.t1"].exons[0].parents == ["g3.t1"]
        assert annot.chrs["chr3"]["g3"].transcripts["g3.t1"].exons[1].parents == ["g3.t1", "g3.t1_1"]

        assert annot.chrs["chr3"]["g3"].transcripts["g3.t1_1"].exons[0].parents == ["g3.t1", "g3.t1_1"]
        

        annot.detect_genes_with_no_transcripts(remove=True, quiet=True)

        assert len(annot.all_gene_ids) == 5

        annot.detect_genes_with_no_transcripts(remove=True, remove_pseudogene=True, quiet=True)

        assert len(annot.all_gene_ids) == 4



# ============================================================
# MicroRNA GFF3 files in human and Arabidopsis format
# ============================================================

class TestAnnotationMiRNAs:
    annot: Annotation
    def test_miRNA_human_format(self, miRNA_human_format_gff3_file):
        annot = Annotation(miRNA_human_format_gff3_file, quiet=True)

        assert len(annot.orphaned_features) == 0

        assert "rna-NR_039937.1" in annot.all_transcript_ids
        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].exons[0].id == "exon-NR_039937.1-1"

        assert len(annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs) == 2

        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs[0].id == "rna-MIR4777"

        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs[0].start == 231362723
        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs[0].end == 231362744

        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs[1].id == "rna-MIR4777-2"

        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs[1].start == 231362764
        assert annot.chrs["NC_000002.12"]["gene-MIR4777"].transcripts["rna-NR_039937.1"].miRNAs[1].end == 231362785

        annot = Annotation(miRNA_human_format_gff3_file, skip_orphaned_features=False, quiet=True)

        assert len(annot.orphaned_features) == 2

        assert annot.orphaned_features[0].id == "exon-MIR4777-1"
        assert annot.orphaned_features[0].start == 231362723
        assert annot.orphaned_features[0].end == 231362744
        assert annot.orphaned_features[1].id == "exon-MIR4777-2-1"
        assert annot.orphaned_features[1].start == 231362764
        assert annot.orphaned_features[1].end == 231362785

    def test_miRNA_arabidopsis_format(self, miRNA_arabidopsis_format_gff3_file):
        annot = Annotation(miRNA_arabidopsis_format_gff3_file, quiet=True)

        assert len(annot.orphaned_features) == 0

        assert "AT4G04095.1" in annot.all_transcript_ids
        assert annot.chrs["Chr4"]["AT4G04095"].transcripts["AT4G04095.1"].exons[0].start == 1023912
        assert annot.chrs["Chr4"]["AT4G04095"].transcripts["AT4G04095.1"].exons[0].end == 1024149

        
        assert annot.chrs["Chr4"]["AT4G04095"].transcripts["AT4G04095.1"].miRNAs[0].id == "ath-miR5635d"
        assert annot.chrs["Chr4"]["AT4G04095"].transcripts["AT4G04095.1"].miRNAs[0].parents == ["AT4G04095.1"]
        assert annot.chrs["Chr4"]["AT4G04095"].transcripts["AT4G04095.1"].miRNAs[0].start == 1024077
        assert annot.chrs["Chr4"]["AT4G04095"].transcripts["AT4G04095.1"].miRNAs[0].end == 1024097
        

# ============================================================
# Testing different merging styles
# ============================================================

class TestAnnotationMerge:
    annot: Annotation
    def test_merge_gff3(self, merge_gff3_file_1, merge_gff3_file_2):
        annot = Annotation(merge_gff3_file_1, quiet=True)
        assert len(annot.all_gene_ids) == 5
        annot2 = Annotation(merge_gff3_file_2, quiet=True)
        assert len(annot2.all_gene_ids) == 6

        annot.merge(annot2, quiet=True)
        assert len(annot.all_gene_ids) == 11

        annot = Annotation(merge_gff3_file_1, quiet=True)
        annot.merge(annot2, rename_clashing_ids=False, quiet=True)
        assert len(annot.all_gene_ids) == 6

        annot = Annotation(merge_gff3_file_1, quiet=True)
        annot.merge(annot2, max_cds_overlap=0, max_exon_overlap=0, max_gene_overlap=0, quiet=True)
        assert len(annot.all_gene_ids) == 7

        annot = Annotation(merge_gff3_file_1, quiet=True)
        annot.merge(annot2, max_cds_overlap=0, quiet=True)
        assert len(annot.all_gene_ids) == 9

        annot = Annotation(merge_gff3_file_1, quiet=True)
        annot.merge(annot2, max_exon_overlap=50, quiet=True)
        assert len(annot.all_gene_ids) == 9

        annot = Annotation(merge_gff3_file_1, quiet=True)
        annot.merge(annot2, max_gene_overlap=50, quiet=True)
        assert len(annot.all_gene_ids) == 8

# ============================================================
# Testing self overlapping genes detection
# ============================================================

class TestAnnotationOverlaps:
    annot: Annotation
    def test_self_overlaps(self, self_overlapping_genes_gff3_file):
        annot = Annotation(self_overlapping_genes_gff3_file, quiet=True)
        annot.overlaps.detect(quiet=True)

        assert annot.overlapped_annotations == set()

        assert set(annot.overlaps.self_genes) == {"g2", "g3", "g4", "g5"}
        assert annot.overlaps.other_genes == set()
        assert annot.chrs["chr1"]["g1"].overlaps == {"self" : [], "other" : []}

        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].id == "g3"
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].CDSs_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].exons_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].orientation == True
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].gene_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].gene_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].min_gene_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].exon_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].exon_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].min_exon_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].min_CDS_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].score == 11

        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].id == "g2"
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].CDSs_in_both == True
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].exons_in_both == True
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].orientation == True
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].gene_query_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].gene_target_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].min_gene_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].exon_query_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].exon_target_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].min_exon_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].min_CDS_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].score == 11

        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].id == "g5"
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].gene_query_percent == 72.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].gene_target_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].min_gene_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].exon_query_percent == 86.9
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].exon_target_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].min_exon_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].CDS_query_percent == 57.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].min_CDS_percent == 57.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].score == 8

        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].id == "g4"
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].gene_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].gene_target_percent == 72.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].min_gene_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].exon_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].exon_target_percent == 86.9
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].min_exon_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].CDS_target_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].min_CDS_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].score == 8

        annot.overlaps.clear(keep_self=True)

        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].id == "g3"
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].CDSs_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].exons_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].orientation == True
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].gene_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].gene_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].min_gene_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].exon_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].exon_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].min_exon_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].min_CDS_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["self"][0].score == 11

        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].id == "g2"
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].CDSs_in_both == True
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].exons_in_both == True
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].orientation == True
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].gene_query_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].gene_target_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].min_gene_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].exon_query_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].exon_target_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].min_exon_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].min_CDS_percent == 100.0
        assert annot.chrs["chr2"]["g3"].overlaps["self"][0].score == 11

        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].id == "g5"
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].gene_query_percent == 72.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].gene_target_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].min_gene_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].exon_query_percent == 86.9
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].exon_target_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].min_exon_percent == 66.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].CDS_query_percent == 57.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].min_CDS_percent == 57.7
        assert annot.chrs["chr3"]["g4"].overlaps["self"][0].score == 8

        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].id == "g4"
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].gene_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].gene_target_percent == 72.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].min_gene_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].exon_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].exon_target_percent == 86.9
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].min_exon_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].CDS_target_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].min_CDS_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["self"][0].score == 8

        annot.overlaps.clear(keep_self=False)

        assert annot.chrs["chr2"]["g2"].overlaps["self"] == []
        

    def test_other_overlaps(self, other_overlapping_genes_gff3_file_1, other_overlapping_genes_gff3_file_2):
        annot = Annotation(other_overlapping_genes_gff3_file_1, quiet=True)
        annot2 = Annotation(other_overlapping_genes_gff3_file_2, quiet=True)
        annot.overlaps.detect(other=annot2, quiet=True)

        assert annot.overlaps.self_genes == set()
        assert annot2.overlaps.self_genes == set()

        assert annot.overlapped_annotations == {"other_overlapping_genes_2"}
        assert annot2.overlapped_annotations == {"other_overlapping_genes_1"}

        assert set(annot.overlaps.other_genes) == {"g3", "g4"}
        assert set(annot2.overlaps.other_genes) == {"g2", "g5"}

        assert annot2.chrs["chr1"]["g1"].overlaps == {"self" : [], "other" : []}

        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].id == "g3"
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].CDSs_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].exons_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].orientation == True
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].gene_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].gene_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].min_gene_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].exon_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].exon_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].min_exon_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].min_CDS_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].score == 11

        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].id == "g2"
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].CDSs_in_both == True
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].exons_in_both == True
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].orientation == True
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].gene_query_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].gene_target_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].min_gene_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].exon_query_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].exon_target_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].min_exon_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].CDS_query_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].CDS_target_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].min_CDS_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].score == 11

        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].id == "g4"
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].gene_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].gene_target_percent == 72.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].min_gene_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].exon_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].exon_target_percent == 86.9
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].min_exon_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].CDS_target_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].min_CDS_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].score == 8

        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].id == "g5"
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].gene_query_percent == 72.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].gene_target_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].min_gene_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].exon_query_percent == 86.9
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].exon_target_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].min_exon_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].CDS_query_percent == 57.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].CDS_target_percent == 100.0
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].min_CDS_percent == 57.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].score == 8

        annot.overlaps.clear(keep_other=True)

        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].id == "g3"
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].CDSs_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].exons_in_both == True
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].orientation == True
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].gene_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].gene_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].min_gene_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].exon_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].exon_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].min_exon_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].CDS_target_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].min_CDS_percent == 100.0
        assert annot.chrs["chr2"]["g2"].overlaps["other"][0].score == 11

        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].id == "g2"
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].CDSs_in_both == True
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].exons_in_both == True
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].orientation == True
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].gene_query_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].gene_target_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].min_gene_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].exon_query_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].exon_target_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].min_exon_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].CDS_query_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].CDS_target_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].min_CDS_percent == 100.0
        assert annot2.chrs["chr2"]["g3"].overlaps["other"][0].score == 11

        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].id == "g4"
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].gene_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].gene_target_percent == 72.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].min_gene_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].exon_query_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].exon_target_percent == 86.9
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].min_exon_percent == 66.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].CDS_query_percent == 100.0
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].CDS_target_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].min_CDS_percent == 57.7
        assert annot.chrs["chr3"]["g5"].overlaps["other"][0].score == 8

        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].id == "g5"
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].gene_query_percent == 72.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].gene_target_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].min_gene_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].exon_query_percent == 86.9
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].exon_target_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].min_exon_percent == 66.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].CDS_query_percent == 57.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].CDS_target_percent == 100.0
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].min_CDS_percent == 57.7
        assert annot2.chrs["chr3"]["g4"].overlaps["other"][0].score == 8

        annot.overlaps.clear(keep_other=False)

        assert annot.chrs["chr2"]["g2"].overlaps["other"] == []

        
# ============================================================
# Testing transcript combining
# ============================================================

class TestTranscriptCombining:
    def test_transcripts_to_combine(self, transcripts_to_combine_gff3_file, transcripts_to_combine_fasta_file):
        g = Genome("transcripts_to_combine", transcripts_to_combine_fasta_file, quiet=True)
        annot = Annotation(transcripts_to_combine_gff3_file, genome=g, quiet=True)
        assert len(annot.chrs["chr1"]["g1"].transcripts) == 1
        assert len(annot.chrs["chr2"]["g2"].transcripts) == 3
        assert len(annot.chrs["chr7"]["g3"].transcripts) == 4
        assert len(annot.chrs["chr17"]["g4"].transcripts) == 2

        assert annot.chrs["chr17"]["g4"].start == 200
        assert annot.chrs["chr17"]["g4"].end == 1800
        assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].start == 200
        assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].end == 1600

        annot.combine_transcripts(quiet=True, redetect_CDS=False)
        assert len(annot.chrs["chr1"]["g1"].transcripts) == 1
        assert len(annot.chrs["chr2"]["g2"].transcripts) == 1
        assert len(annot.chrs["chr7"]["g3"].transcripts) == 1
        assert len(annot.chrs["chr17"]["g4"].transcripts) == 1

        assert "g1_t001" in annot.chrs["chr1"]["g1"].transcripts
        assert "g2_t001" in annot.chrs["chr2"]["g2"].transcripts
        assert "g3_t001" in annot.chrs["chr7"]["g3"].transcripts
        assert "g4_t001" in annot.chrs["chr17"]["g4"].transcripts

        assert len(annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons) == 2
        assert len(annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons) == 3
        assert len(annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons) == 3
        assert len(annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons) == 2

        assert len(annot.chrs["chr1"]["g1"].transcripts["g1_t001"].CDSs) == 1
        assert len(annot.chrs["chr2"]["g2"].transcripts["g2_t001"].CDSs) == 0
        assert len(annot.chrs["chr7"]["g3"].transcripts["g3_t001"].CDSs) == 0
        assert len(annot.chrs["chr17"]["g4"].transcripts["g4_t001"].CDSs) == 1

        assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[0].start == 1000
        assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[0].end == 4000
        assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[1].start == 5000
        assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[1].end == 8000

        assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[0].start == 700
        assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[0].end == 2000
        assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[1].start == 2600
        assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[1].end == 2800
        assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[2].start == 3000
        assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[2].end == 3200

        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[0].start == 100
        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[0].end == 600
        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[1].start == 1000
        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[1].end == 1800
        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[2].start == 2000
        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[2].end == 2400

        assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[0].start == 200
        assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[0].end == 1000
        assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[1].start == 1200
        assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[1].end == 1800

        assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[1].parents == ["g3_t001"]

    # def test_transcripts_to_combine_redetect_CDS(self, transcripts_to_combine_gff3_file, transcripts_to_combine_fasta_file):
    #     g = Genome("transcripts_to_combine", transcripts_to_combine_fasta_file, quiet=True)
    #     annot = Annotation(transcripts_to_combine_gff3_file, genome=g, quiet=True)
    #     assert len(annot.chrs["chr1"]["g1"].transcripts) == 1
    #     assert len(annot.chrs["chr2"]["g2"].transcripts) == 3
    #     assert len(annot.chrs["chr7"]["g3"].transcripts) == 4
    #     assert len(annot.chrs["chr17"]["g4"].transcripts) == 2

    #     assert annot.chrs["chr17"]["g4"].start == 200
    #     assert annot.chrs["chr17"]["g4"].end == 1800
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].start == 200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].end == 1600

    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_mRNA1"].CDSs["g1_CDS1"].CDS_segments[0].start == 1200
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_mRNA1"].CDSs["g1_CDS1"].CDS_segments[0].end == 4000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_mRNA1"].CDSs["g1_CDS1"].CDS_segments[1].start == 5000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_mRNA1"].CDSs["g1_CDS1"].CDS_segments[1].end == 6000

    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].CDSs["g4_CDS1"].CDS_segments[0].start == 200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].CDSs["g4_CDS1"].CDS_segments[0].end == 1000
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].CDSs["g4_CDS1"].CDS_segments[1].start == 1200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_mRNA1"].CDSs["g4_CDS1"].CDS_segments[1].end == 1600

    #     annot.combine_transcripts(quiet=True, redetect_CDS=True, respect_non_coding=True, respect_non_combined=True)
    #     assert len(annot.chrs["chr1"]["g1"].transcripts) == 1
    #     assert len(annot.chrs["chr2"]["g2"].transcripts) == 1
    #     assert len(annot.chrs["chr7"]["g3"].transcripts) == 1
    #     assert len(annot.chrs["chr17"]["g4"].transcripts) == 1

    #     assert "g1_t001" in annot.chrs["chr1"]["g1"].transcripts
    #     assert "g2_t001" in annot.chrs["chr2"]["g2"].transcripts
    #     assert "g3_t001" in annot.chrs["chr7"]["g3"].transcripts
    #     assert "g4_t001" in annot.chrs["chr17"]["g4"].transcripts

    #     assert len(annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons) == 2
    #     assert len(annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons) == 3
    #     assert len(annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons) == 3
    #     assert len(annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons) == 2

    #     assert len(annot.chrs["chr1"]["g1"].transcripts["g1_t001"].CDSs) == 1
    #     assert len(annot.chrs["chr2"]["g2"].transcripts["g2_t001"].CDSs) == 0
    #     assert len(annot.chrs["chr7"]["g3"].transcripts["g3_t001"].CDSs) == 0
    #     assert len(annot.chrs["chr17"]["g4"].transcripts["g4_t001"].CDSs) == 1

    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].CDSs["g1_CDS1"].CDS_segments[0].start == 1200
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].CDSs["g1_CDS1"].CDS_segments[0].end == 4000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].CDSs["g1_CDS1"].CDS_segments[1].start == 5000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].CDSs["g1_CDS1"].CDS_segments[1].end == 6000

    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].CDSs["g4_CDS1"].CDS_segments[0].start != 200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].CDSs["g4_CDS1"].CDS_segments[0].end != 1000
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].CDSs["g4_CDS1"].CDS_segments[1].start != 1200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].CDSs["g4_CDS1"].CDS_segments[1].end != 1600

    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[0].start == 1000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[0].end == 4000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[1].start == 5000
    #     assert annot.chrs["chr1"]["g1"].transcripts["g1_t001"].exons[1].end == 8000

    #     assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[0].start == 700
    #     assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[0].end == 2000
    #     assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[1].start == 2600
    #     assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[1].end == 2800
    #     assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[2].start == 3000
    #     assert annot.chrs["chr2"]["g2"].transcripts["g2_t001"].exons[2].end == 3200

    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[0].start == 100
    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[0].end == 600
    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[1].start == 1000
    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[1].end == 1800
    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[2].start == 2000
    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[2].end == 2400

    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[0].start == 200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[0].end == 1000
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[1].start == 1200
    #     assert annot.chrs["chr17"]["g4"].transcripts["g4_t001"].exons[1].end == 1800

    #     assert annot.chrs["chr7"]["g3"].transcripts["g3_t001"].exons[1].parents == ["g3_t001"]

# ============================================================
# Testing CDS reformatting
# ============================================================

class TestCDSReformatting:
    def test_cds_reformatting(self, multiple_isoforms_gff3_file):
        annot = Annotation(multiple_isoforms_gff3_file, quiet=True)

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["CDS_iso1a"].id == "CDS_iso1a"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["CDS_iso1b"].id == "CDS_iso1b"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["CDS_iso1c"].id == "CDS_iso1c"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["CDS_iso1a"].CDS_segments[0].id == "CDS_iso1a"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["CDS_iso1a"].CDS_segments[1].id == "CDS_iso1a"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["CDS_iso1a"].CDS_segments[2].id == "CDS_iso1a"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["CDS_iso1b"].CDS_segments[0].id == "CDS_iso1b"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["CDS_iso1b"].CDS_segments[1].id == "CDS_iso1b"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["CDS_iso1c"].CDS_segments[0].id == "CDS_iso1c"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["CDS_iso1c"].CDS_segments[1].id == "CDS_iso1c"

        annot.CDS_to_CDS_segment_ids()

        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs) == 1
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs) == 1
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs) == 1

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].id == "mRNA_iso1a_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].id == "mRNA_iso1b_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].id == "mRNA_iso1c_CDS1"

        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments) == 3
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].CDS_segments) == 2
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].CDS_segments) == 2

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments[0].id == "mRNA_iso1a_CDS1_1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments[1].id == "mRNA_iso1a_CDS1_2"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments[2].id == "mRNA_iso1a_CDS1_3"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].CDS_segments[0].id == "mRNA_iso1b_CDS1_1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].CDS_segments[1].id == "mRNA_iso1b_CDS1_2"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].CDS_segments[0].id == "mRNA_iso1c_CDS1_1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].CDS_segments[1].id == "mRNA_iso1c_CDS1_2"

        annot.CDS_segment_to_CDS_ids()

        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs) == 1
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs) == 1
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs) == 1

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].id == "mRNA_iso1a_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].id == "mRNA_iso1b_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].id == "mRNA_iso1c_CDS1"

        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments) == 3
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].CDS_segments) == 2
        assert len(annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].CDS_segments) == 2

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments[0].id == "mRNA_iso1a_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments[1].id == "mRNA_iso1a_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1a"].CDSs["mRNA_iso1a_CDS1"].CDS_segments[2].id == "mRNA_iso1a_CDS1"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].CDS_segments[0].id == "mRNA_iso1b_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1b"].CDSs["mRNA_iso1b_CDS1"].CDS_segments[1].id == "mRNA_iso1b_CDS1"

        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].CDS_segments[0].id == "mRNA_iso1c_CDS1"
        assert annot.chrs["chr1"]["gene_iso1"].transcripts["mRNA_iso1c"].CDSs["mRNA_iso1c_CDS1"].CDS_segments[1].id == "mRNA_iso1c_CDS1"

# ============================================================
# Testing synteny
# ============================================================

class TestSynteny:
    def test_synteny(self, synteny_before_liftover_gff3_file, synteny_after_liftover_gff3_file):

        annot_before = Annotation(synteny_before_liftover_gff3_file, quiet=True, define_synteny=True)

        assert annot_before.chrs["chr1"]["g1"].synteny.previous == None
        assert annot_before.chrs["chr1"]["g2"].synteny.previous == "g1"
        assert annot_before.chrs["chr1"]["g3"].synteny.previous == "g2"
        assert annot_before.chrs["chr1"]["g4"].synteny.previous == "g3"
        assert annot_before.chrs["chr1"]["g5"].synteny.previous == "g4"
        assert annot_before.chrs["chr1"]["g6"].synteny.previous == "g5"
        assert annot_before.chrs["chr1"]["g7"].synteny.previous == "g6"
        assert annot_before.chrs["chr1"]["g8"].synteny.previous == "g7"
        assert annot_before.chrs["chr1"]["g9"].synteny.previous == "g8"

        assert annot_before.chrs["chr1"]["g1"].synteny.next == "g2"
        assert annot_before.chrs["chr1"]["g2"].synteny.next == "g3"
        assert annot_before.chrs["chr1"]["g3"].synteny.next == "g4"
        assert annot_before.chrs["chr1"]["g4"].synteny.next == "g5"
        assert annot_before.chrs["chr1"]["g5"].synteny.next == "g6"
        assert annot_before.chrs["chr1"]["g6"].synteny.next == "g7"
        assert annot_before.chrs["chr1"]["g7"].synteny.next == "g8"
        assert annot_before.chrs["chr1"]["g8"].synteny.next == "g9"
        assert annot_before.chrs["chr1"]["g9"].synteny.next == None

        assert annot_before.chrs["chr1"]["g1"].synteny.order == 0
        assert annot_before.chrs["chr1"]["g2"].synteny.order == 1
        assert annot_before.chrs["chr1"]["g3"].synteny.order == 2
        assert annot_before.chrs["chr1"]["g4"].synteny.order == 3
        assert annot_before.chrs["chr1"]["g5"].synteny.order == 4
        assert annot_before.chrs["chr1"]["g6"].synteny.order == 5
        assert annot_before.chrs["chr1"]["g7"].synteny.order == 6
        assert annot_before.chrs["chr1"]["g8"].synteny.order == 7
        assert annot_before.chrs["chr1"]["g9"].synteny.order == 8

        annot_after = Annotation(synteny_after_liftover_gff3_file, quiet=True, original_annotation=annot_before, define_synteny=True)
        
        assert annot_after.unmapped == ["g6"]

        assert annot_after.chrs["chr1"]["g1"].synteny.old_previous == None
        assert annot_after.chrs["chr1"]["g2"].synteny.old_previous == "g1"
        assert annot_after.chrs["chr1"]["g3"].synteny.old_previous == "g2"
        assert annot_after.chrs["chr1"]["g4"].synteny.old_previous == "g3"
        assert annot_after.chrs["chr1"]["g5"].synteny.old_previous == "g4"
        assert annot_after.chrs["chr1"]["g7"].synteny.old_previous == "g6"
        assert annot_after.chrs["chr1"]["g8"].synteny.old_previous == "g7"
        assert annot_after.chrs["chr1"]["g9"].synteny.old_previous == "g8"

        assert annot_after.chrs["chr1"]["g1"].synteny.old_next == "g2"
        assert annot_after.chrs["chr1"]["g2"].synteny.old_next == "g3"
        assert annot_after.chrs["chr1"]["g3"].synteny.old_next == "g4"
        assert annot_after.chrs["chr1"]["g4"].synteny.old_next == "g5"
        assert annot_after.chrs["chr1"]["g5"].synteny.old_next == "g6"
        assert annot_after.chrs["chr1"]["g7"].synteny.old_next == "g8"
        assert annot_after.chrs["chr1"]["g8"].synteny.old_next == "g9"
        assert annot_after.chrs["chr1"]["g9"].synteny.old_next == None

        assert annot_after.chrs["chr1"]["g1"].synteny.old_order == 0
        assert annot_after.chrs["chr1"]["g2"].synteny.old_order == 1
        assert annot_after.chrs["chr1"]["g3"].synteny.old_order == 2
        assert annot_after.chrs["chr1"]["g4"].synteny.old_order == 3
        assert annot_after.chrs["chr1"]["g5"].synteny.old_order == 4
        assert annot_after.chrs["chr1"]["g7"].synteny.old_order == 6
        assert annot_after.chrs["chr1"]["g8"].synteny.old_order == 7
        assert annot_after.chrs["chr1"]["g9"].synteny.old_order == 8

        assert annot_after.chrs["chr1"]["g1"].synteny.previous == None
        assert annot_after.chrs["chr1"]["g2"].synteny.previous == "g4"
        assert annot_after.chrs["chr1"]["g3"].synteny.previous == "g5"
        assert annot_after.chrs["chr1"]["g4"].synteny.previous == "g3"
        assert annot_after.chrs["chr1"]["g5"].synteny.previous == "g1"
        assert annot_after.chrs["chr1"]["g7"].synteny.previous == "g2"
        assert annot_after.chrs["chr1"]["g8"].synteny.previous == "g7"
        assert annot_after.chrs["chr1"]["g9"].synteny.previous == "g8"

        assert annot_after.chrs["chr1"]["g1"].synteny.next == "g5"
        assert annot_after.chrs["chr1"]["g2"].synteny.next == "g7"
        assert annot_after.chrs["chr1"]["g3"].synteny.next == "g4"
        assert annot_after.chrs["chr1"]["g4"].synteny.next == "g2"
        assert annot_after.chrs["chr1"]["g5"].synteny.next == "g3"
        assert annot_after.chrs["chr1"]["g7"].synteny.next == "g8"
        assert annot_after.chrs["chr1"]["g8"].synteny.next == "g9"
        assert annot_after.chrs["chr1"]["g9"].synteny.next == None

        assert annot_after.chrs["chr1"]["g1"].synteny.order == 0
        assert annot_after.chrs["chr1"]["g2"].synteny.order == 4
        assert annot_after.chrs["chr1"]["g3"].synteny.order == 2
        assert annot_after.chrs["chr1"]["g4"].synteny.order == 3
        assert annot_after.chrs["chr1"]["g5"].synteny.order == 1
        assert annot_after.chrs["chr1"]["g7"].synteny.order == 5
        assert annot_after.chrs["chr1"]["g8"].synteny.order == 6
        assert annot_after.chrs["chr1"]["g9"].synteny.order == 7

        assert annot_after.chrs["chr1"]["g1"].synteny.liftover_conserved == False
        assert annot_after.chrs["chr1"]["g2"].synteny.liftover_conserved == False
        assert annot_after.chrs["chr1"]["g3"].synteny.liftover_conserved == False
        assert annot_after.chrs["chr1"]["g4"].synteny.liftover_conserved == False
        assert annot_after.chrs["chr1"]["g5"].synteny.liftover_conserved == False
        assert annot_after.chrs["chr1"]["g7"].synteny.liftover_conserved == False
        assert annot_after.chrs["chr1"]["g8"].synteny.liftover_conserved == True
        assert annot_after.chrs["chr1"]["g9"].synteny.liftover_conserved == True

# ============================================================
# Testing subset
# ============================================================

class TestSubset:
    def test_subset(self, arabidopsis_araport11_gff3_file):

        annot = Annotation(arabidopsis_araport11_gff3_file, quiet=True)

        assert len(annot.all_gene_ids) == 3000

        annot.subset(gene_cap=300, quiet=True)

        assert len(annot.all_gene_ids) == 300

    def test_subset_no_gene_cap(self, arabidopsis_araport11_gff3_file):
        # gene_cap=None should retain all genes
        annot1 = Annotation(arabidopsis_araport11_gff3_file, quiet=True)
        assert len(annot1.all_gene_ids) == 3000
        annot1.subset(gene_cap=None, quiet=True)
        assert len(annot1.all_gene_ids) == 3000

        # no_gene_cap=True should retain all genes
        annot2 = Annotation(arabidopsis_araport11_gff3_file, quiet=True)
        annot2.subset(no_gene_cap=True, quiet=True)
        assert len(annot2.all_gene_ids) == 3000

        # gene_cap=0 should also disable the gene cap
        annot3 = Annotation(arabidopsis_araport11_gff3_file, quiet=True)
        annot3.subset(gene_cap=0, quiet=True)
        assert len(annot3.all_gene_ids) == 3000

    def test_subset_flexible_inputs(self, arabidopsis_araport11_gff3_file):
        # Passing list or tuple for chosen_features should work seamlessly
        annot_list = Annotation(arabidopsis_araport11_gff3_file, quiet=True)
        annot_list.subset(chosen_features=["Chr4"], no_gene_cap=True, quiet=True)
        assert len(annot_list.all_gene_ids) == 3000

        annot_tuple = Annotation(arabidopsis_araport11_gff3_file, quiet=True)
        annot_tuple.subset(chosen_features=("Chr4",), gene_cap=150, quiet=True)
        assert len(annot_tuple.all_gene_ids) == 150

    def test_subset_chr_cap_and_seed(self, test_data_dir):
        gff_path = str(test_data_dir / "input/annotation/for_merge_2.gff3")

        # chr_cap=2 with seed
        annot_a = Annotation(gff_path, quiet=True)
        chosen_a = annot_a.subset(chr_cap=2, no_gene_cap=True, no_min_genes=True, seed=42, quiet=True)
        assert len(chosen_a) == 2
        assert len(annot_a.chrs) == 2

        # Same seed should pick exact same chromosomes
        annot_b = Annotation(gff_path, quiet=True)
        chosen_b = annot_b.subset(chr_cap=2, no_gene_cap=True, no_min_genes=True, seed=42, quiet=True)
        assert chosen_a == chosen_b


# ============================================================
# Testing rework CDS
# ============================================================

class TestReworkCDS:
    def test_rework_cds(self, arabidopsis_tair10_fasta_file, arabidopsis_araport11_no_CDS_gff3_file, arabidopsis_araport11_with_CDS_gff3_file, tmp_path):

        output_dir = tmp_path

        genome = Genome("TAIR10", arabidopsis_tair10_fasta_file)
        annot = Annotation(annot_file_path=arabidopsis_araport11_no_CDS_gff3_file, name="araport11_no_CDS", rework_all_CDSs=True, genome=genome, quiet=True)

        annot.export.gff(output_dir=tmp_path, subfolder=False, quiet=True)

        # Compare contents
        with open(arabidopsis_araport11_with_CDS_gff3_file, "r") as f:
            expected_content = f.read()
            
        with open(output_dir / "araport11_no_CDS_on_TAIR10_clean.gff3", "r") as f:
            generated_content = f.read()
            
        assert generated_content == expected_content, f"Output mismatch for {arabidopsis_araport11_with_CDS_gff3_file}"