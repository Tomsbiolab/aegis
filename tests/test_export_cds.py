import pytest
from pathlib import Path
from aegis.annotation import Annotation
from aegis.genome import Genome


def test_export_cds_protein_oriented_and_raw(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    # 1. Default export: protein_oriented=True
    annot.export.CDSs(output_dir=str(tmp_path), verbose=True, quiet=True)
    default_file = tmp_path / f"{annot.id}_CDSs_c_id_main_coordinates.fasta"
    assert default_file.exists()
    default_text = default_file.read_text()
    assert "_raw" not in default_file.name

    # 2. Raw export: protein_oriented=False
    annot.export.CDSs(output_dir=str(tmp_path), verbose=True, protein_oriented=False, quiet=True)
    raw_file = tmp_path / f"{annot.id}_CDSs_raw_c_id_main_coordinates.fasta"
    assert raw_file.exists()
    assert "_raw" in raw_file.name

    # 3. Unique CDSs default
    annot.export.unique_CDSs(output_dir=str(tmp_path), quiet=True)
    unique_default = tmp_path / f"{annot.id}_unique_CDSs.fasta"
    assert unique_default.exists()

    # 4. Unique CDSs raw
    annot.export.unique_CDSs(output_dir=str(tmp_path), protein_oriented=False, quiet=True)
    unique_raw = tmp_path / f"{annot.id}_unique_CDSs_raw.fasta"
    assert unique_raw.exists()


def test_export_unique_proteins_per_gene(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    annot.export.proteins(output_dir=str(tmp_path), only_main=False, verbose=False, unique_proteins_per_gene=False, quiet=True)
    all_prot_file = tmp_path / f"{annot.id}_proteins_p_id_all.fasta"
    assert all_prot_file.exists()
    all_text = all_prot_file.read_text()

    annot.export.proteins(output_dir=str(tmp_path), only_main=False, verbose=False, unique_proteins_per_gene=True, quiet=True)
    uniq_prot_file = tmp_path / f"{annot.id}_proteins_p_id_unique_per_gene.fasta"
    assert uniq_prot_file.exists()
    uniq_text = uniq_prot_file.read_text()

    all_headers = [line for line in all_text.splitlines() if line.startswith(">")]
    uniq_headers = [line for line in uniq_text.splitlines() if line.startswith(">")]
    assert len(uniq_headers) < len(all_headers)


def test_export_proteins_detailed_headers(tmp_path):
    # Minus-strand CDS over two exons: mRNA ATGAAA|CCCTAA -> MKP*
    fasta_path = tmp_path / "genome.fasta"
    fasta_path.write_text(">chrX\nGGTTAGGGCCCCTTTCATGG\n")
    gff3_path = tmp_path / "annot.gff3"
    gff3_path.write_text(
        "##gff-version 3\n"
        "chrX\ttest\tgene\t3\t18\t.\t-\t.\tID=g1\n"
        "chrX\ttest\tmRNA\t3\t18\t.\t-\t.\tID=t1;Parent=g1\n"
        "chrX\ttest\texon\t3\t8\t.\t-\t.\tID=e1;Parent=t1\n"
        "chrX\ttest\texon\t13\t18\t.\t-\t.\tID=e2;Parent=t1\n"
        "chrX\ttest\tCDS\t3\t8\t.\t-\t0\tID=c1;Parent=t1;exception=ribosomal slippage\n"
        "chrX\ttest\tCDS\t13\t18\t.\t-\t0\tID=c1;Parent=t1;exception=ribosomal slippage\n"
    )
    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    annot.export.proteins(output_dir=str(tmp_path), verbose=True, quiet=True)
    header, seq = (tmp_path / f"{annot.id}_proteins_p_id_main_coordinates.fasta").read_text().splitlines()
    assert header == ">c1.prot|exception=ribosomal slippage|readthrough=end|table=1|-|chrX|location=complement(join(3..8,13..18))"
    assert seq == "MKP"

    annot.export.proteins(output_dir=str(tmp_path), verbose=False, quiet=True)
    assert (tmp_path / f"{annot.id}_proteins_p_id_main.fasta").read_text().splitlines()[0] == ">c1.prot"


def test_cds_transl_table_overrides_chromosome_table(tmp_path):
    # ATG TGA TGG AGA: TGA is W and AGA a stop in the vertebrate mitochondrial code (2), on a nuclear chromosome
    fasta_path = tmp_path / "genome.fasta"
    fasta_path.write_text(">chrX\nATGTGATGGAGAATGTGATGGAGA\n")
    gff3_path = tmp_path / "annot.gff3"
    gff3_path.write_text(
        "##gff-version 3\n"
        "chrX\ttest\tgene\t1\t24\t.\t+\t.\tID=g1\n"
        "chrX\ttest\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
        "chrX\ttest\texon\t1\t12\t.\t+\t.\tID=e1;Parent=t1\n"
        "chrX\ttest\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=t1;transl_table=2\n"
        "chrX\ttest\tmRNA\t13\t24\t.\t+\t.\tID=t2;Parent=g1\n"
        "chrX\ttest\texon\t13\t24\t.\t+\t.\tID=e2;Parent=t2\n"
        "chrX\ttest\tCDS\t13\t24\t.\t+\t0\tID=c2;Parent=t2\n"
    )
    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)
    c1 = annot.chrs["chrX"]["g1"].transcripts["t1"].CDSs["c1"]
    c2 = annot.chrs["chrX"]["g1"].transcripts["t2"].CDSs["c2"]
    assert (annot.translation_table("chrX", c1), annot.translation_table("chrX", c2)) == (2, 1)

    annot.generate_proteins(quiet=True)
    assert (c1.protein.seq, c1.protein.tags["table"]) == ("MWW*", "2")
    assert (c2.protein.seq, c2.protein.tags["table"]) == ("M*WR", "1")

    # Raw CDS export strips the stop codon of each CDS's own table
    annot.export.CDSs(output_dir=str(tmp_path), only_main=False, verbose=False, protein_oriented=False, strip_stop=True, quiet=True)
    out = (tmp_path / f"{annot.id}_CDSs_raw_c_id_all.fasta").read_text().splitlines()
    assert dict(zip(out[::2], out[1::2])) == {">c1": "ATGTGATGG", ">c2": "ATGTGATGGAGA"}


def test_set_genetic_codes_clears_proteins_and_exports_follow(tmp_path):
    # ATG TGA TGG TAA: TGA is a stop in table 1 but W in the vertebrate mitochondrial code (2)
    fasta_path = tmp_path / "genome.fasta"
    fasta_path.write_text(">chrX\nCCCATGTGATGGTAACCC\n")
    gff3_path = tmp_path / "annot.gff3"
    gff3_path.write_text(
        "##gff-version 3\n"
        "chrX\ttest\tgene\t4\t15\t.\t+\t.\tID=g1\n"
        "chrX\ttest\tmRNA\t4\t15\t.\t+\t.\tID=t1;Parent=g1\n"
        "chrX\ttest\texon\t4\t15\t.\t+\t.\tID=e1;Parent=t1\n"
        "chrX\ttest\tCDS\t4\t15\t.\t+\t0\tID=c1;Parent=t1\n"
    )
    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)
    out_file = tmp_path / f"{annot.id}_proteins_p_id_main.fasta"

    assert annot.translation_table("chrX") == 1
    annot.export.proteins(output_dir=str(tmp_path), verbose=False, quiet=True)
    assert out_file.read_text().splitlines()[1] == "M*W"

    annot.set_genetic_codes(taxonomy="vertebrate", mitochondria_chroms=["chrX"])
    assert annot.contains_protein_sequences is False
    assert annot.translation_table("chrX") == 2
    annot.export.proteins(output_dir=str(tmp_path), verbose=False, quiet=True)
    assert out_file.read_text().splitlines()[1] == "MWW"

    # Re-applying the same configuration keeps the proteins
    annot.set_genetic_codes(taxonomy="vertebrate", mitochondria_chroms="chrX")
    assert annot.contains_protein_sequences is True

    with pytest.raises(ValueError, match="Specified chloroplast chromosome 'missing' was not found"):
        annot.set_genetic_codes(chloroplast_chroms=["missing"])


def test_export_unique_CDSs_per_gene(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    # 1. All CDSs
    annot.export.CDSs(output_dir=str(tmp_path), only_main=False, verbose=False, unique_CDSs_per_gene=False, quiet=True)
    all_cds_file = tmp_path / f"{annot.id}_CDSs_c_id_all.fasta"
    assert all_cds_file.exists()
    all_text = all_cds_file.read_text()

    # 2. Unique CDSs per gene (protein-oriented)
    annot.export.CDSs(output_dir=str(tmp_path), only_main=False, verbose=False, unique_CDSs_per_gene=True, protein_oriented=True, quiet=True)
    uniq_cds_file = tmp_path / f"{annot.id}_CDSs_c_id_unique_per_gene.fasta"
    assert uniq_cds_file.exists()
    uniq_text = uniq_cds_file.read_text()

    all_headers = [line for line in all_text.splitlines() if line.startswith(">")]
    uniq_headers = [line for line in uniq_text.splitlines() if line.startswith(">")]
    assert len(uniq_headers) <= len(all_headers)

    # 3. Unique CDSs per gene (raw)
    annot.export.CDSs(output_dir=str(tmp_path), only_main=False, verbose=False, unique_CDSs_per_gene=True, protein_oriented=False, quiet=True)
    uniq_raw_file = tmp_path / f"{annot.id}_CDSs_raw_c_id_unique_per_gene.fasta"
    assert uniq_raw_file.exists()


def test_export_cds_and_proteins_with_taxonomy(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    # Every export works with each taxonomy preset set on the annotation
    for taxonomy in ("plant", "vertebrate", "yeast", "invertebrate"):
        annot.set_genetic_codes(taxonomy=taxonomy)
        annot.export.proteins(output_dir=str(tmp_path), quiet=True)
        annot.export.unique_proteins(output_dir=str(tmp_path), quiet=True)
        annot.export.CDSs(output_dir=str(tmp_path), quiet=True)
        annot.export.unique_CDSs(output_dir=str(tmp_path), quiet=True)


def test_export_protein_strip_stop_default(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    # Default strip_stop is True: no trailing '*'
    dir_stripped = tmp_path / "stripped"
    annot.export.proteins(output_dir=str(dir_stripped), verbose=False, quiet=True)
    prot_file = dir_stripped / f"{annot.id}_proteins_p_id_main.fasta"
    seqs = [line.strip() for line in prot_file.read_text().splitlines() if line and not line.startswith(">")]
    assert not any(s.endswith("*") for s in seqs)

    # Explicit strip_stop=False: contains trailing '*'
    dir_kept = tmp_path / "kept"
    annot.export.proteins(output_dir=str(dir_kept), verbose=False, strip_stop=False, quiet=True)
    prot_file_kept = dir_kept / f"{annot.id}_proteins_p_id_main.fasta"
    seqs_kept = [line.strip() for line in prot_file_kept.read_text().splitlines() if line and not line.startswith(">")]
    assert any(s.endswith("*") for s in seqs_kept)


def test_export_cds_strip_stop_organelle(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    # Export CDS with strip_stop=True vs strip_stop=False
    dir_no_strip = tmp_path / "cds_nostrip"
    dir_strip = tmp_path / "cds_strip"

    annot.export.CDSs(output_dir=str(dir_no_strip), strip_stop=False, quiet=True)
    annot.export.CDSs(output_dir=str(dir_strip), strip_stop=True, quiet=True)

    cds_nostrip_file = dir_no_strip / f"{annot.id}_CDSs_c_id_main_coordinates.fasta"
    cds_strip_file = dir_strip / f"{annot.id}_CDSs_c_id_main_coordinates.fasta"

    seqs_nostrip = [line.strip() for line in cds_nostrip_file.read_text().splitlines() if line and not line.startswith(">")]
    seqs_strip = [line.strip() for line in cds_strip_file.read_text().splitlines() if line and not line.startswith(">")]

    assert any(len(ns) == len(s) + 3 for ns, s in zip(seqs_nostrip, seqs_strip))




@pytest.mark.parametrize("export", ["CDSs", "proteins"])
def test_simple_headers_use_each_cds_own_transcript_id(tmp_path, export):
    """Without detailed headers only the chosen id is written, and transcript ids come from each CDS's own transcript."""
    seq = "ATGAAATTTGGGCCCTAA" + "A" * 40 + "ATGCCCTAA" + "A" * 20
    (tmp_path / "g.fasta").write_text(f">chr1\n{seq}\n")
    (tmp_path / "g.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t1\t67\t.\t+\t.\tID=g1\n"
        "chr1\tt\tmRNA\t1\t18\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tt\texon\t1\t18\t.\t+\t.\tParent=t1\n"
        "chr1\tt\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tt\tmRNA\t59\t67\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\tt\texon\t59\t67\t.\t+\t.\tParent=t2\n"
        "chr1\tt\tCDS\t59\t67\t.\t+\t0\tID=c2;Parent=t2\n"
    )
    genome = Genome(name="g", genome_file_path=str(tmp_path / "g.fasta"), quiet=True)
    annot = Annotation(name="g", annot_file_path=str(tmp_path / "g.gff3"), genome=genome, quiet=True)

    expected = {"CDS": [">c1", ">c2"], "transcript": [">t1", ">t2"], "gene": [">g1"]}
    for used_id, headers in expected.items():
        out = tmp_path / f"{used_id}.fasta"
        getattr(annot.export, export)(only_main=used_id == "gene", verbose=False, used_id=used_id, filepath=str(out))
        assert [line for line in out.read_text().splitlines() if line.startswith(">")] == headers


@pytest.mark.parametrize("seq, strip_stop, expected", [
    ("MK*", True, "MK"),
    ("MK*", False, "MK*"),
    ("M*K*", True, "MXK"),
    ("M*K*", False, "MXK*"),
    ("MK**", False, "MKX*"),  # the stop right before the terminal one is internal
    ("M*K", True, "MXK"),     # 3' partial: no terminal stop to remove
    ("M*K", False, "MXK"),    # nor to add back
    ("MK", False, "MK"),
])
def test_protein_export_seq_internal_stops_as_x(seq, strip_stop, expected):
    from aegis.annotation_components.export import _protein_export_seq
    assert _protein_export_seq(seq, strip_stop, "x") == expected
    assert _protein_export_seq(seq, strip_stop, "keep") == (seq[:-1] if strip_stop and seq.endswith("*") else seq)


def _internal_stops_annotation(tmp_path):
    # t1: complete with an internal stop (M*K*); t2: 3' partial with an internal stop (M*K); t3: complete (MK*)
    seq = "ATGTAAAAATAA" + "ATGTAAAAA" + "ATGAAATAA"
    fasta = tmp_path / "internal_stops.fasta"
    fasta.write_text(f">chr1\n{seq}\n")
    lines = ["##gff-version 3", f"chr1\ttest\tgene\t1\t{len(seq)}\t.\t+\t.\tID=g1"]
    for n, (start, end) in enumerate(((1, 12), (13, 21), (22, 30)), 1):
        lines += [
            f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID=t{n};Parent=g1",
            f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tID=e{n};Parent=t{n}",
            f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tID=c{n};Parent=t{n}",
        ]
    gff = tmp_path / "internal_stops.gff3"
    gff.write_text("\n".join(lines) + "\n")
    genome = Genome(name="internal_stops", genome_file_path=str(fasta), quiet=True)
    return Annotation(name="internal_stops", annot_file_path=str(gff), genome=genome, quiet=True)


def _read_fasta(path):
    records, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[0]
            records[name] = ""
        elif line:
            records[name] += line.strip()
    return records


@pytest.mark.parametrize("strip_stop, internal_stops, expected", [
    (True, "keep", ["M*K", "M*K", "MK"]),
    (False, "keep", ["M*K*", "M*K", "MK*"]),
    (True, "x", ["MXK", "MXK", "MK"]),
    (False, "x", ["MXK*", "MXK", "MK*"]),
    (True, "skip", ["MK"]),
    (False, "skip", ["MK*"]),
])
def test_export_proteins_internal_stops(tmp_path, strip_stop, internal_stops, expected):
    annot = _internal_stops_annotation(tmp_path)
    annot.export.proteins(output_dir=str(tmp_path), filename="all.fasta", only_main=False, only_cds_main=False, used_id="CDS", verbose=False, strip_stop=strip_stop, internal_stops=internal_stops, quiet=True)
    annot.export.unique_proteins(output_dir=str(tmp_path), filename="unique.fasta", strip_stop=strip_stop, internal_stops=internal_stops, quiet=True)
    assert list(_read_fasta(tmp_path / "all.fasta").values()) == expected
    assert sorted(_read_fasta(tmp_path / "unique.fasta").values()) == sorted(expected)


def test_export_proteins_invalid_internal_stops(tmp_path):
    annot = _internal_stops_annotation(tmp_path)
    with pytest.raises(ValueError):
        annot.export.proteins(output_dir=str(tmp_path), internal_stops="drop", quiet=True)
    with pytest.raises(ValueError):
        annot.export.unique_proteins(output_dir=str(tmp_path), internal_stops="drop", quiet=True)
