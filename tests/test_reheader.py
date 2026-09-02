import pytest
from typer.testing import CliRunner

from aegis.genome import Genome
from aegis.annotation import Annotation
from aegis.cli.tidy_genome import app as tidy_genome_app
from aegis.cli.extract import app as extract_app
from aegis.cli.subset import app as subset_app
from aegis.cli.split import app as split_app
from aegis.cli.summary import app as summary_app

runner = CliRunner()

GWH_FASTA_CONTENT = """>GWHBAVD00000001\tOriSeqID=chr1\tLen=120
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>GWHBAVD00000002\tOriSeqID=chr2\tLen=120
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>GWHBAVD00000003\tOriSeqID=ctg100000\tLen=60
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
"""

GWH_GFF_CONTENT = """##gff-version 3
chr1\tTEST\tgene\t1\t60\t.\t+\t.\tID=gene1;Name=gene1
chr1\tTEST\tmRNA\t1\t60\t.\t+\t.\tID=tx1;Parent=gene1
chr1\tTEST\texon\t1\t60\t.\t+\t.\tID=exon1;Parent=tx1
chr1\tTEST\tCDS\t1\t60\t.\t+\t0\tID=cds1;Parent=tx1
chr2\tTEST\tgene\t1\t60\t.\t+\t.\tID=gene2;Name=gene2
chr2\tTEST\tmRNA\t1\t60\t.\t+\t.\tID=tx2;Parent=gene2
chr2\tTEST\texon\t1\t60\t.\t+\t.\tID=exon2;Parent=tx2
chr2\tTEST\tCDS\t1\t60\t.\t+\t0\tID=cds2;Parent=tx2
ctg100000\tTEST\tgene\t1\t30\t.\t+\t.\tID=gene3;Name=gene3
ctg100000\tTEST\tmRNA\t1\t30\t.\t+\t.\tID=tx3;Parent=gene3
ctg100000\tTEST\texon\t1\t30\t.\t+\t.\tID=exon3;Parent=tx3
ctg100000\tTEST\tCDS\t1\t30\t.\t+\t0\tID=cds3;Parent=tx3
"""

@pytest.fixture
def gwh_sample_files(tmp_path):
    fasta_path = tmp_path / "gwh_sample.fasta"
    gff_path = tmp_path / "gwh_sample.gff3"
    fasta_path.write_text(GWH_FASTA_CONTENT, encoding="utf-8")
    gff_path.write_text(GWH_GFF_CONTENT, encoding="utf-8")
    return fasta_path, gff_path


def test_genome_gwh_reheader(gwh_sample_files):
    fasta_path, _ = gwh_sample_files
    
    # 1. Without reheadering (default)
    g_raw = Genome(name="raw", genome_file_path=str(fasta_path), quiet=True)
    assert "GWHBAVD00000001" in g_raw.scaffolds
    assert "chr1" not in g_raw.scaffolds
    assert g_raw.scaffolds["GWHBAVD00000001"].chromosome is False

    # 2. With gwh=True
    g_gwh = Genome(name="gwh", genome_file_path=str(fasta_path), gwh=True, quiet=True)
    assert "chr1" in g_gwh.scaffolds
    assert "chr2" in g_gwh.scaffolds
    assert "ctg100000" in g_gwh.scaffolds
    assert "GWHBAVD00000001" not in g_gwh.scaffolds
    # Alias lookup
    assert g_gwh.get_scaffold("GWHBAVD00000001") is not None
    assert g_gwh.get_scaffold("GWHBAVD00000001").name == "chr1"
    assert "GWHBAVD00000001" in g_gwh
    # Chromosome detection
    assert g_gwh.scaffolds["chr1"].chromosome is True
    assert g_gwh.scaffolds["chr2"].chromosome is True
    assert g_gwh.scaffolds["ctg100000"].chromosome is False
    assert g_gwh.scaffolds["chr1"].original_name == "GWHBAVD00000001"

    # 3. With header_id_tag="OriSeqID"
    g_tag = Genome(name="tag", genome_file_path=str(fasta_path), header_id_tag="OriSeqID", quiet=True)
    assert "chr1" in g_tag.scaffolds

    # 4. With header_id_regex
    g_reg = Genome(name="reg", genome_file_path=str(fasta_path), header_id_regex=r"OriSeqID=([^\s]+)", quiet=True)
    assert "chr1" in g_reg.scaffolds


def test_genome_export_descriptions(gwh_sample_files, tmp_path):
    fasta_path, _ = gwh_sample_files
    g = Genome(name="gwh", genome_file_path=str(fasta_path), gwh=True, quiet=True)
    
    out_default = tmp_path / "out_default.fasta"
    g.export(filepath=str(out_default), quiet=True, keep_description=False)
    lines_default = [line.strip() for line in out_default.read_text().splitlines() if line.startswith(">")]
    assert lines_default == [">chr1", ">chr2", ">ctg100000"]

    out_keep = tmp_path / "out_keep.fasta"
    g.export(filepath=str(out_keep), quiet=True, keep_description=True)
    lines_keep = [line.strip() for line in out_keep.read_text().splitlines() if line.startswith(">")]
    # Leading token must be the new sequence ID
    assert lines_keep[0].startswith(">chr1")
    assert "OriSeqID=chr1" in lines_keep[0]


def test_annotation_compatibility_diagnostic(gwh_sample_files, tmp_path):
    fasta_path, gff_path = gwh_sample_files
    
    # Raw genome (0 overlap with GFF, but OriSeqID in header) raises ValueError with OriSeqID guidance
    g_raw = Genome(name="raw", genome_file_path=str(fasta_path), quiet=True)
    with pytest.raises(ValueError, match="match 'OriSeqID' found in FASTA header descriptions"):
        Annotation(annot_file_path=str(gff_path), genome=g_raw, quiet=True)

    # Raw genome (0 overlap with GFF, no OriSeqID in header) raises ValueError
    no_overlap_fasta = tmp_path / "no_overlap.fasta"
    no_overlap_fasta.write_text(">scaffold_999\nATGCATGCATGC\n", encoding="utf-8")
    g_no_overlap = Genome(name="no_overlap", genome_file_path=str(no_overlap_fasta), quiet=True)
    with pytest.raises(ValueError, match=r"None of the 3 chromosome IDs in annotation .* match the 1 scaffold IDs in genome"):
        Annotation(annot_file_path=str(gff_path), genome=g_no_overlap, quiet=True)

    # GWH genome (perfect overlap) succeeds
    g_gwh = Genome(name="gwh", genome_file_path=str(fasta_path), gwh=True, quiet=True)
    annot_ok = Annotation(annot_file_path=str(gff_path), genome=g_gwh, quiet=True)
    assert annot_ok is not None


def test_cli_tidy_genome_gwh(gwh_sample_files, tmp_path):
    fasta_path, gff_path = gwh_sample_files
    out_dir = tmp_path / "tidy_out"

    args = [
        str(fasta_path),
        str(gff_path),
        "--gwh",
        "--remove-scaffolds",
        "-d", str(out_dir),
        "-og", "tidy.fasta",
        "-oa", "tidy.gff3"
    ]
    result = runner.invoke(tidy_genome_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    out_fa = out_dir / "tidy.fasta"
    out_gff = out_dir / "tidy.gff3"
    assert out_fa.exists()
    assert out_gff.exists()

    fa_headers = [line.strip() for line in out_fa.read_text().splitlines() if line.startswith(">")]
    # Scaffold ctg100000 should be removed, chr1 and chr2 kept
    assert fa_headers == [">chr1", ">chr2"]

    gff_content = out_gff.read_text()
    assert "chr1" in gff_content
    assert "chr2" in gff_content
    assert "ctg100000" not in gff_content


def test_cli_extract_gwh(gwh_sample_files, tmp_path):
    fasta_path, gff_path = gwh_sample_files
    out_dir = tmp_path / "extract_out"

    args = [
        str(gff_path),
        str(fasta_path),
        "--gwh",
        "-f", "gene",
        "-d", str(out_dir),
        "-q"
    ]
    result = runner.invoke(extract_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    extracted_files = list(out_dir.glob("*.fasta"))
    assert len(extracted_files) >= 1
    content = extracted_files[0].read_text()
    assert ">gene1" in content
    assert ">gene2" in content


def test_cli_subset_gwh(gwh_sample_files, tmp_path):
    fasta_path, gff_path = gwh_sample_files
    out_dir = tmp_path / "subset_out"

    # 1. Test subset without --gwh gives informative error
    args_fail = [
        str(gff_path),
        str(fasta_path),
        "-c", "chr1",
        "-d", str(out_dir),
        "-q"
    ]
    result_fail = runner.invoke(subset_app, args_fail)
    assert result_fail.exit_code != 0
    assert "OriSeqID" in str(result_fail.exception)

    # 2. Test subset with --gwh succeeds
    args_ok = [
        str(gff_path),
        str(fasta_path),
        "--gwh",
        "-c", "chr1",
        "--no-gene-cap",
        "-d", str(out_dir),
        "-q"
    ]
    result_ok = runner.invoke(subset_app, args_ok)
    assert result_ok.exit_code == 0, f"Error: {result_ok.stdout}"

    out_fa = out_dir / "gwh_sample_subset.fasta"
    out_gff = out_dir / "gwh_sample_subset.gff3"
    assert out_fa.exists()
    assert out_gff.exists()
    assert ">chr1" in out_fa.read_text()
    assert "chr2" not in out_fa.read_text()


def test_cli_split_gwh(gwh_sample_files, tmp_path):
    fasta_path, gff_path = gwh_sample_files
    out_dir = tmp_path / "split_out"

    args = [
        str(gff_path),
        "-g", str(fasta_path),
        "--gwh",
        "-s", "chr1,chr2",
        "-d", str(out_dir),
        "-q"
    ]
    result = runner.invoke(split_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"

    fa_c1 = out_dir / "gwh_sample_splitchr1.fasta"
    gff_c1 = out_dir / "gwh_sample_splitchr1.gff3"
    assert fa_c1.exists()
    assert gff_c1.exists()
    assert ">chr1" in fa_c1.read_text()
    assert "chr1" in gff_c1.read_text()


def test_cli_summary_gwh(gwh_sample_files, tmp_path):
    fasta_path, gff_path = gwh_sample_files
    out_dir = tmp_path / "summary_out"

    args = [
        str(gff_path),
        str(fasta_path),
        "--gwh",
        "-d", str(out_dir),
        "-q"
    ]
    result = runner.invoke(summary_app, args)
    assert result.exit_code == 0, f"Error: {result.stdout}"
    assert (out_dir / "gwh_sample_on_gwh_sample_full_stats.csv").exists()
