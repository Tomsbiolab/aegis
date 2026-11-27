from aegis.annotation import Annotation
from aegis.genome import Genome


def _parse_ids_and_parents(gff_path):
    ids = set()
    parent_refs = []
    with open(gff_path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            attrs = dict(item.split("=", 1) for item in cols[8].split(";") if "=" in item)
            if "ID" in attrs:
                ids.add(attrs["ID"])
            for key in ("Parent", "Derives_from"):
                if key in attrs:
                    for parent_id in attrs[key].split(","):
                        parent_refs.append((cols[2], attrs.get("ID"), key, parent_id))
    return ids, parent_refs


def test_export_gff_main_only_has_no_dangling_parents_and_preserves_in_memory_parents(
    shared_parents_gff3_file, tmp_path
):
    annot = Annotation(shared_parents_gff3_file, quiet=True)
    before_parents = {
        (g.id, t.id, e.id): list(e.parents)
        for genes in annot.chrs.values()
        for g in genes.values()
        for t in g.transcripts.values()
        for e in t.exons
    }
    assert any(len(parents) > 1 for parents in before_parents.values())

    out_path = tmp_path / "main_only.gff3"
    annot.export.gff(filepath=str(out_path), main_only=True, UTRs=True, quiet=True)

    ids, parent_refs = _parse_ids_and_parents(out_path)
    dangling = [ref for ref in parent_refs if ref[3] not in ids]
    assert dangling == []

    after_parents = {
        (g.id, t.id, e.id): list(e.parents)
        for genes in annot.chrs.values()
        for g in genes.values()
        for t in g.transcripts.values()
        for e in t.exons
    }
    assert after_parents == before_parents


def test_for_dapseq_handles_already_dapfit_genome_and_exclude_non_coding(tmp_path):
    fasta_path = tmp_path / "dap.fa"
    fasta_path.write_text(">chr01\n" + "ACGT" * 250 + "\n", encoding="utf-8")

    gff_path = tmp_path / "annot.gff3"
    gff_path.write_text(
        "##gff-version 3\n"
        "chr01\tsrc\tgene\t1\t300\t.\t+\t.\tID=g_coding\n"
        "chr01\tsrc\tmRNA\t1\t300\t.\t+\t.\tID=t_coding;Parent=g_coding\n"
        "chr01\tsrc\texon\t1\t300\t.\t+\t.\tID=e1;Parent=t_coding\n"
        "chr01\tsrc\tCDS\t1\t300\t.\t+\t0\tID=c1;Parent=t_coding\n"
        "chr01\tsrc\tgene\t400\t600\t.\t+\t.\tID=g_nc\n"
        "chr01\tsrc\tlnc_RNA\t400\t600\t.\t+\t.\tID=t_nc;Parent=g_nc\n"
        "chr01\tsrc\texon\t400\t600\t.\t+\t.\tID=e2;Parent=t_nc\n",
        encoding="utf-8",
    )

    genome = Genome("dap", str(fasta_path), quiet=True)
    assert genome.dapfit is True
    annot = Annotation(str(gff_path), genome=genome, quiet=True)

    out_gff = tmp_path / "out_for_dap.gff3"
    annot.export.for_dapseq(
        genome=genome,
        genome_output_dir=str(tmp_path),
        gff_filepath=str(out_gff),
        exclude_non_coding=True,
    )

    rows = [
        line.split("\t")
        for line in out_gff.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    feature_types = [r[2] for r in rows]
    ids, parent_refs = _parse_ids_and_parents(out_gff)

    assert "g_coding" in ids
    assert "t_coding" in ids
    assert "e1" in ids
    assert "c1" in ids
    assert "g_nc" not in ids
    assert "t_nc" not in ids
    assert "e2" not in ids
    assert {"gene", "mRNA", "exon", "CDS"} <= set(feature_types)
    assert [ref for ref in parent_refs if ref[3] not in ids] == []

