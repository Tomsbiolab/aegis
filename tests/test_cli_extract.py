import pytest
from pathlib import Path
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
