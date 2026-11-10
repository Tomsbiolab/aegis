from typer.testing import CliRunner

from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.cli.extract import app

runner = CliRunner()

# extract_test: gene2 has mRNA2/3/4 sharing a TSS (identical promoters);
# mRNA1 starts at base 1, so it has an empty promoter and is never written.


def _records(path):
    records, header = [], None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
        elif header is not None:
            records.append((header, line))
            header = None
    return records


def _load(test_data_dir):
    genome = Genome(name="test_genome", genome_file_path=str(test_data_dir / "input/fasta/extract_test.fasta"), quiet=True)
    return Annotation(name="test_annot", annot_file_path=str(test_data_dir / "input/annotation/extract_test.gff3"), genome=genome, quiet=True)


def _extract(test_data_dir, out_dir, *extra):
    res = runner.invoke(app, [
        str(test_data_dir / "input/annotation/extract_test.gff3"),
        str(test_data_dir / "input/fasta/extract_test.fasta"),
        "-f", "promoter", "-d", str(out_dir), "-q", *extra,
    ])
    assert res.exit_code == 0, res.stdout
    files = list(out_dir.iterdir())
    assert len(files) == 1, files
    return files[0]


def test_cli_promoter_filename_uses_requested_size_and_type(test_data_dir, tmp_path):
    out = _extract(test_data_dir, tmp_path / "a", "-p", "upstream_ATG", "-ps", "1500", "-m", "all")
    assert out.name.endswith("_promoters_prom_id_all_1500_upstream_ATG.fasta")


def test_cli_promoter_filename_defaults(test_data_dir, tmp_path):
    out = _extract(test_data_dir, tmp_path / "a", "-m", "all")
    assert out.name.endswith("_promoters_prom_id_all_2000_standard.fasta")


def test_cli_promoter_modes_collapse_identical_sequences(test_data_dir, tmp_path):
    all_recs = _records(_extract(test_data_dir, tmp_path / "all", "-m", "all"))
    assert sorted(h for h, _ in all_recs) == ["mRNA2_promoter", "mRNA3_promoter", "mRNA4_promoter"]
    assert len({s for _, s in all_recs}) == 1

    per_gene = _extract(test_data_dir, tmp_path / "upg", "-m", "unique_per_gene")
    assert "_unique_per_gene_" in per_gene.name
    assert _records(per_gene) == all_recs[:1]

    unique = _extract(test_data_dir, tmp_path / "u", "-m", "unique")
    assert "_unique_promoters_2000_standard" in unique.name
    assert _records(unique) == all_recs[:1]


def test_export_promoters_regenerates_on_new_size_or_type(test_data_dir, tmp_path):
    annot = _load(test_data_dir)

    annot.export.promoters(only_main=False, verbose=False, output_dir=str(tmp_path), promoter_size=100, quiet=True)
    first = _records(tmp_path / f"{annot.id}_promoters_prom_id_all_100_standard.fasta")
    assert all(len(s) == 100 for _, s in first)

    annot.export.promoters(only_main=False, verbose=False, output_dir=str(tmp_path), promoter_size=50, quiet=True)
    second = _records(tmp_path / f"{annot.id}_promoters_prom_id_all_50_standard.fasta")
    assert all(len(s) == 50 for _, s in second)


def test_export_promoters_reuses_pre_generated_promoters(test_data_dir, tmp_path):
    annot = _load(test_data_dir)
    annot.generate_promoters(promoter_size=80, promoter_type="upstream_ATG")

    annot.export.promoters(only_main=False, verbose=False, output_dir=str(tmp_path), quiet=True)
    recs = _records(tmp_path / f"{annot.id}_promoters_prom_id_all_80_upstream_ATG.fasta")
    assert recs and all(len(s) == 80 for _, s in recs)
