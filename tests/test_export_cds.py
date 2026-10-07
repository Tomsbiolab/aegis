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


