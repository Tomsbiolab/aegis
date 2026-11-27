from aegis.annotation import Annotation


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
