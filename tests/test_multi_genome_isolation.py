"""Tests for per-Annotation genome isolation (no class-global _ACTIVE_GENOME corruption)."""

import pickle
import pytest

from aegis.annotation import Annotation
from aegis.feature import Feature
from aegis.genome import Genome


def _write_fasta(path, records: dict[str, str]) -> str:
    lines = []
    for name, seq in records.items():
        lines.append(f">{name}\n{seq}\n")
    path.write_text("".join(lines), encoding="utf-8")
    return str(path)


def _write_gff(path, chrom: str, gene_id: str, tx_id: str, cds_id: str, exon1=(1, 12), exon2=(21, 32), strand="+") -> str:
    g_start = min(exon1[0], exon2[0])
    g_end = max(exon1[1], exon2[1])
    content = (
        "##gff-version 3\n"
        f"{chrom}\ttest\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\tID={gene_id}\n"
        f"{chrom}\ttest\tmRNA\t{g_start}\t{g_end}\t.\t{strand}\t.\tID={tx_id};Parent={gene_id}\n"
        f"{chrom}\ttest\texon\t{exon1[0]}\t{exon1[1]}\t.\t{strand}\t.\tID={tx_id}_e1;Parent={tx_id}\n"
        f"{chrom}\ttest\texon\t{exon2[0]}\t{exon2[1]}\t.\t{strand}\t.\tID={tx_id}_e2;Parent={tx_id}\n"
        f"{chrom}\ttest\tCDS\t{exon1[0]}\t{exon1[1]}\t.\t{strand}\t0\tID={cds_id};Parent={tx_id}\n"
        f"{chrom}\ttest\tCDS\t{exon2[0]}\t{exon2[1]}\t.\t{strand}\t0\tID={cds_id};Parent={tx_id}\n"
    )
    path.write_text(content, encoding="utf-8")
    return str(path)


class TestMultiGenomeIsolation:
    """Verify multiple Annotation instances bound to distinct Genomes never cross-contaminate."""

    def test_annotation_init_does_not_mutate_class_globals(self, tmp_path):
        prev_g = Feature._ACTIVE_GENOME
        prev_hg = Feature._ACTIVE_HARD_GENOME
        try:
            Feature._ACTIVE_GENOME = None
            Feature._ACTIVE_HARD_GENOME = None

            fa = _write_fasta(tmp_path / "g1.fa", {"chr1": "ATGAAATTTGGGGTAAAAAGCCCAAATAA" + "A" * 20})
            hfa = _write_fasta(tmp_path / "h1.fa", {"chr1": "NNN" * 20})
            gff = _write_gff(tmp_path / "a1.gff3", "chr1", "g1", "t1", "c1")

            g = Genome("g1", fa, quiet=True)
            hg = Genome("h1", hfa, quiet=True)
            annot = Annotation(gff, name="a1", genome=g, hard_masked_genome=hg, quiet=True)

            assert Feature._ACTIVE_GENOME is None
            assert Feature._ACTIVE_HARD_GENOME is None
            assert annot.genome is g
            assert annot.hard_masked_genome is hg
        finally:
            Feature._ACTIVE_GENOME = prev_g
            Feature._ACTIVE_HARD_GENOME = prev_hg

    def test_two_annotations_same_chrom_distinct_sequences(self, tmp_path):
        # Genome A: exon1 (1..12) = ATG AAA TTT GGG (MKFG), intron (13..20) = GTAAAAAG (GT-AG), exon2 (21..32) = CCC AAA CCC TAA (PKP*)
        seq_a = "ATGAAATTTGGG" + "GTAAAAAG" + "CCCAAACCCTAA" + "A" * 20
        # Genome B: exon1 (1..12) = ATG GGG CCC TTT (MGPF), intron (13..20) = GCAAAAAG (GC-AG), exon2 (21..32) = TTT GGG TTT TGA (FGF*)
        seq_b = "ATGGGGCCCTTT" + "GCAAAAAG" + "TTTGGGTTTTGA" + "C" * 20

        fa_a = _write_fasta(tmp_path / "ga.fa", {"chr1": seq_a})
        fa_b = _write_fasta(tmp_path / "gb.fa", {"chr1": seq_b})
        gff_a = _write_gff(tmp_path / "aa.gff3", "chr1", "gA", "tA", "cA")
        gff_b = _write_gff(tmp_path / "ab.gff3", "chr1", "gB", "tB", "cB")

        ga = Genome("ga", fa_a, quiet=True)
        gb = Genome("gb", fa_b, quiet=True)

        annot_a = Annotation(gff_a, name="aa", genome=ga, quiet=True)
        annot_b = Annotation(gff_b, name="ab", genome=gb, quiet=True)

        # Direct .seq / .seqs on annot_a after loading annot_b
        gene_a = annot_a.chrs["chr1"]["gA"]
        tx_a = gene_a.transcripts["tA"]
        cds_a = next(iter(tx_a.CDSs.values()))

        gene_b = annot_b.chrs["chr1"]["gB"]
        tx_b = gene_b.transcripts["tB"]
        cds_b = next(iter(tx_b.CDSs.values()))

        assert gene_a.seq == seq_a[:32]
        assert gene_b.seq == seq_b[:32]
        assert tx_a.seq == "ATGAAATTTGGGCCCAAACCCTAA"
        assert tx_b.seq == "ATGGGGCCCTTTTTTGGGTTTTGA"
        assert cds_a.seq == "ATGAAATTTGGGCCCAAACCCTAA"
        assert cds_b.seq == "ATGGGGCCCTTTTTTGGGTTTTGA"
        assert tx_a.seqs[0] == "ATGAAATTTGGGCCCAAACCCTAA"
        assert cds_a.seqs[0] == "ATGAAATTTGGGCCCAAACCCTAA"

        # generate_proteins() on annot_a must translate ga, not gb
        annot_a.generate_proteins(quiet=True)
        annot_b.generate_proteins(quiet=True)
        assert cds_a.protein is not None and cds_a.protein.seq == "MKFGPKP*"
        assert cds_b.protein is not None and cds_b.protein.seq == "MGPFFGF*"

        # stats.update() (introns + splice sites + GC content)
        annot_a.stats.update(export=False, quiet=True)
        annot_b.stats.update(export=False, quiet=True)
        assert tx_a.introns[0].boundary == "GT-AG"
        assert tx_b.introns[0].boundary == "GC-AG"
        assert annot_a.stats.data["intron-exon boundary: GT-AG"] == 1
        assert annot_a.stats.data["intron-exon boundary: GC-AG"] == 0
        assert annot_b.stats.data["intron-exon boundary: GC-AG"] == 1
        assert annot_b.stats.data["intron-exon boundary: GT-AG"] == 0

        # Promoters generated on annot_a must slice ga
        annot_a.generate_promoters(promoter_size=10)
        annot_b.generate_promoters(promoter_size=10)
        assert tx_a.promoter._genome is ga
        assert tx_b.promoter._genome is gb

        # Export methods on annot_a must slice ga
        out_dir = tmp_path / "exports"
        out_dir.mkdir()
        annot_a.clear_proteins()
        annot_a.export.proteins(output_dir=str(out_dir), filename="a_prot.fa", strip_stop=False, verbose=False)
        annot_a.export.CDSs(output_dir=str(out_dir), filename="a_cds.fa", verbose=False)
        annot_a.export.transcripts(output_dir=str(out_dir), filename="a_tx.fa", verbose=False)
        annot_a.export.genes(output_dir=str(out_dir), filename="a_gene.fa", verbose=False)

        assert "MKFGPKP*" in (out_dir / "a_prot.fa").read_text()
        assert "ATGAAATTTGGGCCCAAACCCTAA" in (out_dir / "a_cds.fa").read_text()
        assert "ATGAAATTTGGGCCCAAACCCTAA" in (out_dir / "a_tx.fa").read_text()
        assert seq_a[:32] in (out_dir / "a_gene.fa").read_text()

    def test_two_annotations_disjoint_scaffolds_and_rework_cdss(self, tmp_path):
        seq_a = "CCCATGAAATTTGGGTAA" + "A" * 20
        seq_b = "GGGATGGGGCCCTGACCC" + "C" * 20

        fa_a = _write_fasta(tmp_path / "ga.fa", {"scfA": seq_a})
        fa_b = _write_fasta(tmp_path / "gb.fa", {"scfB": seq_b})
        gff_a = tmp_path / "aa.gff3"
        gff_a.write_text(
            "##gff-version 3\n"
            "scfA\ttest\tgene\t1\t18\t.\t+\t.\tID=gA\n"
            "scfA\ttest\tmRNA\t1\t18\t.\t+\t.\tID=tA;Parent=gA\n"
            "scfA\ttest\texon\t1\t18\t.\t+\t.\tID=eA;Parent=tA\n"
        )
        gff_b = tmp_path / "ab.gff3"
        gff_b.write_text(
            "##gff-version 3\n"
            "scfB\ttest\tgene\t1\t18\t.\t+\t.\tID=gB\n"
            "scfB\ttest\tmRNA\t1\t18\t.\t+\t.\tID=tB;Parent=gB\n"
            "scfB\ttest\texon\t1\t18\t.\t+\t.\tID=eB;Parent=tB\n"
        )

        ga = Genome("ga", fa_a, quiet=True)
        gb = Genome("gb", fa_b, quiet=True)
        annot_a = Annotation(str(gff_a), name="aa", genome=ga, quiet=True)
        annot_b = Annotation(str(gff_b), name="ab", genome=gb, quiet=True)

        # Reworking CDSs on annot_a after annot_b was loaded must still find scfA's ORF
        annot_a.rework_CDSs(override=True, quiet=True)
        tx_a = annot_a.chrs["scfA"]["gA"].transcripts["tA"]
        assert tx_a.coding is True
        cds_a = next(iter(tx_a.CDSs.values()))
        assert cds_a.protein is not None
        assert cds_a.protein.seq == "MKFG*"
        # Newly created CDS segments and UTRs must also be bound to ga
        assert cds_a.seq == "ATGAAATTTGGGTAA"
        assert all(u._genome is ga for u in cds_a.UTRs)

        # Reworking CDSs on annot_b finds scfB's ORF
        annot_b.rework_CDSs(override=True, quiet=True)
        cds_b = next(iter(annot_b.chrs["scfB"]["gB"].transcripts["tB"].CDSs.values()))
        assert cds_b.protein is not None
        assert cds_b.protein.seq == "MGP*"

    def test_genome_none_annotation_never_slices_previous_genome(self, tmp_path):
        seq = "ATGAAATTTGGG" + "GTAAAAAG" + "CCCAAACCCTAA" + "A" * 20
        fa = _write_fasta(tmp_path / "g1.fa", {"chr1": seq})
        hfa = _write_fasta(tmp_path / "h1.fa", {"chr1": "N" * len(seq)})
        gff = _write_gff(tmp_path / "a1.gff3", "chr1", "g1", "t1", "c1")

        g1 = Genome("g1", fa, quiet=True)
        h1 = Genome("h1", hfa, quiet=True)
        annot_with_genome = Annotation(gff, name="with_g", genome=g1, hard_masked_genome=h1, quiet=True)
        annot_without_genome = Annotation(gff, name="no_g", genome=None, hard_masked_genome=None, quiet=True)

        # annot_without_genome MUST raise ValueError, never silently slice g1/h1
        g_none = annot_without_genome.chrs["chr1"]["g1"]
        t_none = g_none.transcripts["t1"]
        c_none = next(iter(t_none.CDSs.values()))

        with pytest.raises(ValueError, match="No genome loaded"):
            _ = g_none.seq
        with pytest.raises(ValueError, match="No genome loaded"):
            _ = t_none.seq
        with pytest.raises(ValueError, match="No genome loaded"):
            _ = c_none.seq
        with pytest.raises(ValueError, match="No genome loaded"):
            annot_without_genome.generate_proteins(quiet=True)
        with pytest.raises(ValueError, match="No hard masked genome loaded"):
            _ = g_none.hard_seq

        # Loading annot_without_genome must NOT affect annot_with_genome
        g_ok = annot_with_genome.chrs["chr1"]["g1"]
        assert g_ok.seq == seq[:32]
        assert g_ok.hard_seq == "N" * 32

    def test_hard_masked_genome_isolation(self, tmp_path):
        seq = "ATGAAATTTGGG" + "GTAAAAAG" + "CCCAAACCCTAA"
        fa = _write_fasta(tmp_path / "g.fa", {"chr1": seq})
        hfa_a = _write_fasta(tmp_path / "ha.fa", {"chr1": "N" * len(seq)})
        hfa_b = _write_fasta(tmp_path / "hb.fa", {"chr1": "A" * len(seq)})
        gff = _write_gff(tmp_path / "a.gff3", "chr1", "g1", "t1", "c1")

        g = Genome("g", fa, quiet=True)
        ha = Genome("ha", hfa_a, quiet=True)
        hb = Genome("hb", hfa_b, quiet=True)

        annot_a = Annotation(gff, name="a", genome=g, hard_masked_genome=ha, quiet=True)
        annot_b = Annotation(gff, name="b", genome=g, hard_masked_genome=hb, quiet=True)

        annot_a.stats.calculate_transcript_masking()
        annot_b.stats.calculate_transcript_masking()

        assert annot_a.chrs["chr1"]["g1"].quality.masked_fraction == 1.0
        assert annot_b.chrs["chr1"]["g1"].quality.masked_fraction == 0.0

    def test_dynamic_genome_reassignment(self, tmp_path):
        seq_a = "ATGAAATTTGGG" + "GTAAAAAG" + "CCCAAACCCTAA"
        seq_b = "ATGGGGCCCTTT" + "GTAAAAAG" + "TTTGGGTTTTGA"
        fa_a = _write_fasta(tmp_path / "ga.fa", {"chr1": seq_a})
        fa_b = _write_fasta(tmp_path / "gb.fa", {"chr1": seq_b})
        gff = _write_gff(tmp_path / "a.gff3", "chr1", "g1", "t1", "c1")

        ga = Genome("ga", fa_a, quiet=True)
        gb = Genome("gb", fa_b, quiet=True)

        annot = Annotation(gff, name="a", genome=ga, quiet=True)
        assert annot.chrs["chr1"]["g1"].seq == seq_a

        # Reassigning annot.genome immediately updates all features in O(1)
        annot.genome = gb
        assert annot.chrs["chr1"]["g1"].seq == seq_b
        assert annot.chrs["chr1"]["g1"].transcripts["t1"].exons[0].seq == seq_b[:12]

        annot.genome = None
        with pytest.raises(ValueError, match="No genome loaded"):
            _ = annot.chrs["chr1"]["g1"].seq

    def test_copy_merge_pickle_and_identical(self, tmp_path):
        seq_a = "ATGAAATTTGGG" + "GTAAAAAG" + "CCCAAACCCTAA"
        seq_b = "ATGGGGCCCTTT" + "GTAAAAAG" + "TTTGGGTTTTGA"
        fa_a = _write_fasta(tmp_path / "ga.fa", {"chr1": seq_a})
        fa_b = _write_fasta(tmp_path / "gb.fa", {"chr1": seq_b})
        gff_a = _write_gff(tmp_path / "aa.gff3", "chr1", "gA", "tA", "cA")
        gff_b = _write_gff(tmp_path / "ab.gff3", "chr1", "gB", "tB", "cB")

        ga = Genome("ga", fa_a, quiet=True)
        gb = Genome("gb", fa_b, quiet=True)
        annot_a = Annotation(gff_a, name="aa", genome=ga, quiet=True)
        annot_b = Annotation(gff_b, name="ab", genome=gb, quiet=True)

        # 1. gene.copy() must NOT deep-copy the Genome
        gene_a = annot_a.chrs["chr1"]["gA"]
        gene_copy = gene_a.copy()
        assert gene_copy._genome_context is gene_a._genome_context
        assert gene_copy.seq == seq_a

        # 2. annot.copy() has an independent _GenomeContext
        annot_a_copy = annot_a.copy()
        assert annot_a_copy._genome_context is not annot_a._genome_context
        annot_a_copy.genome = None
        with pytest.raises(ValueError):
            _ = annot_a_copy.chrs["chr1"]["gA"].seq
        assert annot_a.chrs["chr1"]["gA"].seq == seq_a

        # 3. merge() rebinds genes copied from other to self's genome context
        annot_merged = annot_a.copy()
        annot_merged.merge(annot_b, quiet=True)
        assert annot_merged.chrs["chr1"]["gB"]._genome_context is annot_merged._genome_context
        assert annot_merged.chrs["chr1"]["gB"].seq == annot_merged.genome.scaffolds["chr1"].seq

        # 4. pickle round-trip preserves genome & rebinds features
        restored = pickle.loads(pickle.dumps(annot_a))
        assert restored.genome is not None
        assert restored.chrs["chr1"]["gA"].seq == seq_a

        # 5. Feature.identical() ignores _genome_context
        annot_a_dup = Annotation(gff_a, name="aa2", genome=None, quiet=True)
        seg_with_g = next(iter(annot_a.chrs["chr1"]["gA"].transcripts["tA"].CDSs.values())).CDS_segments[0]
        seg_no_g = next(iter(annot_a_dup.chrs["chr1"]["gA"].transcripts["tA"].CDSs.values())).CDS_segments[0]
        assert seg_with_g._genome_context is not seg_no_g._genome_context
        assert seg_with_g.identical(seg_no_g)

    def test_sort_genes_multiprocess_keeps_genome_binding(self, tmp_path):
        seq = "ATGAAATTTGGG" + "GTAAAAAG" + "CCCAAACCCTAA"
        fa = _write_fasta(tmp_path / "g.fa", {"chr1": seq, "chr2": seq})
        gff = tmp_path / "a.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t32\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t1\t32\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\ttest\texon\t1\t32\t.\t+\t.\tID=e1;Parent=t1\n"
            "chr2\ttest\tgene\t1\t32\t.\t+\t.\tID=g2\n"
            "chr2\ttest\tmRNA\t1\t32\t.\t+\t.\tID=t2;Parent=g2\n"
            "chr2\ttest\texon\t1\t32\t.\t+\t.\tID=e2;Parent=t2\n"
        )
        g = Genome("g", fa, quiet=True)
        annot = Annotation(str(gff), name="a", genome=g, quiet=True)
        annot.sort_genes(processes=2, quiet=True)
        assert annot.chrs["chr1"]["g1"].seq == seq
        assert annot.chrs["chr2"]["g2"].seq == seq
