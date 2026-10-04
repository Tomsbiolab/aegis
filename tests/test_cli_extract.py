import pytest
from typer.testing import CliRunner

from aegis.cli.extract import app

runner = CliRunner()

# Parametrized test cases mapping generated outputs to expected reference files
@pytest.mark.parametrize(
    "options, expected_filename",
    [
        # Gene
        (["-f", "gene"], "extract_test_on_extract_test_genes.fasta"),
        
        # CDS
        (["-f", "CDS", "-m", "all", "--feature-id", "CDS"], "extract_test_on_extract_test_CDSs_c_id_all.fasta"),
        (["-f", "CDS", "-m", "main", "--feature-id", "CDS"], "extract_test_on_extract_test_CDSs_c_id_main.fasta"),
        (["-f", "CDS", "-m", "unique_per_gene", "--feature-id", "CDS"], "extract_test_on_extract_test_CDSs_c_id_unique_per_gene.fasta"),
        (["-f", "CDS", "-m", "unique", "--feature-id", "CDS"], "extract_test_on_extract_test_unique_CDSs.fasta"),
        
        # Protein
        (["-f", "protein", "-m", "all", "--feature-id", "feature"], "extract_test_on_extract_test_proteins_p_id_all.fasta"),
        (["-f", "protein", "-m", "main", "--feature-id", "feature"], "extract_test_on_extract_test_proteins_p_id_main.fasta"),
        (["-f", "protein", "-m", "unique_per_gene", "--feature-id", "feature"], "extract_test_on_extract_test_proteins_p_id_unique_per_gene.fasta"),
        (["-f", "protein", "-m", "unique", "--feature-id", "feature"], "extract_test_on_extract_test_unique_proteins.fasta"),
        
        # Transcript
        (["-f", "transcript", "-m", "all", "--feature-id", "transcript"], "extract_test_on_extract_test_transcripts_t_id_all.fasta"),
        (["-f", "transcript", "-m", "main", "--feature-id", "transcript"], "extract_test_on_extract_test_transcripts_t_id_main.fasta"),
        (["-f", "transcript", "-m", "unique_per_gene", "--feature-id", "transcript"], "extract_test_on_extract_test_transcripts_t_id_unique_per_gene.fasta"),
        (["-f", "transcript", "-m", "unique", "--feature-id", "transcript"], "extract_test_on_extract_test_unique_transcripts.fasta"),
    ]
)
def test_aegis_extract_cli(test_data_dir, tmp_path, options, expected_filename):
    """
    Test the aegis extract CLI with different options and compare with reference outputs.
    """
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    
    output_dir = tmp_path / "aegis_output" / "features"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    args = [
        str(gff3_path),
        str(fasta_path),
        "-a", "extract_test",
        "-g", "extract_test",
        "-d", str(output_dir),
        "-q"
    ] + options
    
    result = runner.invoke(app, args)
    
    assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}. Output: {result.stdout}"
    
    # Check generated files against reference expected output
    expected_file_path = test_data_dir / "features_output" / expected_filename
    generated_file_path = output_dir / expected_filename
    
    assert generated_file_path.exists(), f"Expected output file not generated: {generated_file_path}"
    
    # Compare contents
    with open(expected_file_path, "r") as f:
        expected_content = f.read()
        
    with open(generated_file_path, "r") as f:
        generated_content = f.read()
        
    assert generated_content == expected_content, f"Output mismatch for {expected_filename}"


def test_extract_cli_translation_options(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "aegis_output" / "features"

    args = [
        str(gff3_path),
        str(fasta_path),
        "-f", "protein",
        "-d", str(output_dir),
        "-gc", "1",
        "--auto-organelle-codes",
        "--mito-code", "2",
        "--plastid-code", "11",
        "--adjust-internal-shifts", "intra_exon",
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"


def test_extract_cli_unknown_mito_contig(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "aegis_output" / "features"

    args = [
        str(gff3_path),
        str(fasta_path),
        "-f", "protein",
        "-d", str(output_dir),
        "--mitochondria-chroms", "nonexistent_mito",
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code != 0
    assert "Specified mitochondrial chromosome 'nonexistent_mito' was not found" in str(result.exception or result.stdout)


def test_extract_cli_invalid_adjust_shifts(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "aegis_output" / "features"

    args = [
        str(gff3_path),
        str(fasta_path),
        "-f", "protein",
        "-d", str(output_dir),
        "--adjust-internal-shifts", "bad_shift_mode",
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code != 0
    assert "Invalid adjust_internal_shifts" in result.output


def test_extract_cli_keep_stop(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "aegis_output" / "features"

    args = [
        str(gff3_path),
        str(fasta_path),
        "-a", "extract_test",
        "-g", "extract_test",
        "-f", "protein",
        "-m", "main",
        "--keep-stop",
        "-d", str(output_dir),
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"
    output_file = output_dir / "extract_test_on_extract_test_proteins_p_id_main.fasta"
    assert output_file.exists()
    content = output_file.read_text()
    # When --keep-stop is passed, protein sequences should have trailing '*'
    seq_lines = [line.strip() for line in content.splitlines() if line and not line.startswith(">")]
    assert any(s.endswith("*") for s in seq_lines)


def test_extract_cli_strip_stop_cds(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    dir_default = tmp_path / "default"
    dir_stripped = tmp_path / "stripped"

    # Default CDS export (retains stop codon)
    args_default = [
        str(gff3_path),
        str(fasta_path),
        "-a", "extract_test",
        "-g", "extract_test",
        "-f", "CDS",
        "-m", "main",
        "--feature-id", "CDS",
        "-d", str(dir_default),
        "-q",
    ]
    res_def = runner.invoke(app, args_default)
    assert res_def.exit_code == 0

    # Stripped CDS export
    args_stripped = [
        str(gff3_path),
        str(fasta_path),
        "-a", "extract_test",
        "-g", "extract_test",
        "-f", "CDS",
        "-m", "main",
        "--feature-id", "CDS",
        "--strip-stop-cds",
        "-d", str(dir_stripped),
        "-q",
    ]
    res_strip = runner.invoke(app, args_stripped)
    assert res_strip.exit_code == 0

    file_def = dir_default / "extract_test_on_extract_test_CDSs_c_id_main.fasta"
    file_strip = dir_stripped / "extract_test_on_extract_test_CDSs_c_id_main.fasta"
    assert file_def.exists() and file_strip.exists()

    seqs_def = [line.strip() for line in file_def.read_text().splitlines() if line and not line.startswith(">")]
    seqs_strip = [line.strip() for line in file_strip.read_text().splitlines() if line and not line.startswith(">")]
    assert len(seqs_def) == len(seqs_strip)
    # At least one CDS with a stop codon is 3 nt shorter
    assert any(len(d) == len(s) + 3 for d, s in zip(seqs_def, seqs_strip))


def test_extract_cli_skip_coordinate_polishing(test_data_dir, tmp_path):
    gff3_path = test_data_dir / "input/annotation/extract_test.gff3"
    fasta_path = test_data_dir / "input/fasta/extract_test.fasta"
    output_dir = tmp_path / "aegis_output" / "features"

    args = [
        str(gff3_path),
        str(fasta_path),
        "-f", "gene",
        "--skip-coordinate-polishing",
        "-d", str(output_dir),
        "-q",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

