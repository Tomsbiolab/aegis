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


def test_export_cds_with_table_none(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"

    genome = Genome(name="test_genome", genome_file_path=str(fasta_path), quiet=True)
    annot = Annotation(name="test_annot", annot_file_path=str(gff3_path), genome=genome, quiet=True)

    annot.export.CDSs(output_dir=str(tmp_path), table=None, quiet=True)
    annot.export.unique_CDSs(output_dir=str(tmp_path), table=None, quiet=True)

