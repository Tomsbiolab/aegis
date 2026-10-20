"""
Example protein tests: each example GFF (NCBI RefSeq GRCh38 loci trimmed into input/fasta/)
is translated and every protein is compared with its expected sequence and tags. Expected
sequences are what the GFF and genome give (NCBI translated_cds), which differ from the RefSeq
proteins where the mRNA differs from the genome (FOXO6, MUC1, the MOB2 NP_443731.2 leader).
"""

import functools
from pathlib import Path

import pytest

from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.misc_features import Protein
from aegis.utils.genefunctions import reverse_complement


TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"

SLIPPAGE_MINUS = {"frameshifts": "-1", "exception": "ribosomal slippage"}
SLIPPAGE_PLUS = {"frameshifts": "+1", "exception": "ribosomal slippage"}
GENOMIC_GAP = {"frameshifts": "+1", "exception": "annotated by transcript or proteomic data"}

# Example -> protein id -> tags besides readthrough and table
EXPECTED_TAGS = {
    "ribosomal_slippage_negative": {
        "NP_001165908.1": SLIPPAGE_MINUS,
        "NP_001171890.1": {"alt_start": "CTG", **SLIPPAGE_MINUS},
        "NP_055883.2": SLIPPAGE_MINUS,
        "NP_001165909.1": {},
        "NP_001035242.1": {},
        "NP_001171891.1": {"alt_start": "CTG"},
    },
    "ribosomal_slippage_positive": {
        "NP_004143.1": SLIPPAGE_PLUS,
        "NP_001287949.1": SLIPPAGE_PLUS,
        "NP_002528.1": SLIPPAGE_PLUS,
        "NP_001288231.1": SLIPPAGE_PLUS,
    },
    "genomic_gap_adjustment": {
        "NP_001278210.2": GENOMIC_GAP,
        "NP_001358649.2": GENOMIC_GAP,
    },
    "translation_exception": {
        "NP_001034936.1": {"transl_except": "(pos:217..219,aa:Sec)"},
        "NP_065184.2": {"transl_except": "(pos:379..381,aa:Sec),(pos:1384..1386,aa:Sec)"},
        "YP_003024026.1": {"transl_except": "(pos:955..956,aa:TERM)", "polyA_bases": "1", "table": "2"},
    },
    "human_partial_cds": {
        "XP_047295430.1": {"partial": "5'"},
        "NP_001392892.1": {},
        "XP_024305386.1": {"partial": "3'"},
        "NP_751944.1": {},
        "NP_443731.2": {"partial": "5'", "exception": "annotated by transcript or proteomic data"},
        "NP_001165694.1": {},
        "XP_011535754.3": {"partial": "3'"},
    },
}

PROTEINS = [(example, protein_id) for example, proteins in EXPECTED_TAGS.items() for protein_id in proteins]


@functools.cache
def generated_proteins(example: str) -> dict[str, Protein]:
    genome = Genome(name=example, genome_file_path=str(TEST_DATA_DIR / f"input/fasta/{example}.fasta"), quiet=True)
    annot = Annotation(name=example, annot_file_path=str(TEST_DATA_DIR / f"input/annotation/{example}.gff3"), genome=genome, quiet=True)
    annot.generate_proteins(quiet=True)
    return {c.id.removeprefix("cds-"): c.protein for genes in annot.chrs.values() for g in genes.values() for t in g.transcripts.values() for c in t.CDSs.values()}


@functools.cache
def expected_proteins(example: str) -> dict[str, str]:
    lines = (TEST_DATA_DIR / f"features_output/{example}_proteins.fasta").read_text().splitlines()
    return {header[1:].split()[0]: seq for header, seq in zip(lines[::2], lines[1::2])}


@pytest.mark.parametrize("example", EXPECTED_TAGS)
def test_every_protein_is_expected(example):
    assert set(generated_proteins(example)) == set(expected_proteins(example)) == set(EXPECTED_TAGS[example])


@pytest.mark.parametrize("example, protein_id", PROTEINS)
def test_protein_sequence(example, protein_id):
    assert generated_proteins(example)[protein_id].seq.rstrip("*") == expected_proteins(example)[protein_id]


@pytest.mark.parametrize("example, protein_id", PROTEINS)
def test_protein_tags(example, protein_id):
    expected = {"readthrough": "end", "table": "1", **EXPECTED_TAGS[example][protein_id]}
    assert generated_proteins(example)[protein_id].tags == expected


@pytest.mark.parametrize("example, protein_id", PROTEINS)
def test_nuc_seq_is_the_codons_of_the_protein(example, protein_id):
    protein = generated_proteins(example)[protein_id]
    assert len(protein.nuc_seq) == 3 * len(protein.seq)


def test_polyadenylation_completes_the_stop_in_nuc_seq():
    # ND1's genome ends in TA; its stop codon is TAA in the mRNA
    protein = generated_proteins("translation_exception")["YP_003024026.1"]
    assert protein.nuc_seq.endswith("ACCTAA") and protein.seq.endswith("T*")


@pytest.mark.parametrize("strip_stop, length, end", [(False, 957, "CAAACCTAA"), (True, 954, "CAAACC")])
def test_exported_CDS_with_a_stop_completed_by_polyadenylation(tmp_path, strip_stop, length, end):
    # ND1 ends ...CAA ACC TA(A): exported with its full stop codon, or without it, keeping the ACC (Thr) codon
    example = "translation_exception"
    genome = Genome(name=example, genome_file_path=str(TEST_DATA_DIR / f"input/fasta/{example}.fasta"), quiet=True)
    annot = Annotation(name=example, annot_file_path=str(TEST_DATA_DIR / f"input/annotation/{example}.gff3"), genome=genome, quiet=True)
    for export, filename in ((annot.export.CDSs, "cds.fasta"), (annot.export.unique_CDSs, "unique_cds.fasta")):
        export(output_dir=str(tmp_path), filename=filename, strip_stop=strip_stop, quiet=True)
        records = dict(r.split("\n", 1) for r in (tmp_path / filename).read_text().split(">")[1:])
        seq = next(v for k, v in records.items() if "YP_003024026.1" in k).replace("\n", "")
        assert len(seq) == length and seq.endswith(end)


def _rebuild(genome_seq: str, header: str, stop_stripped: bool) -> str:
    """The sequence a detailed CDS header describes: its location's genomic bases, plus the polyA_bases A's that end a kept stop."""
    fields = header.split("|")
    strand = fields[-3]
    location = fields[-1].removeprefix("location=").removeprefix("complement(").removeprefix("join(").rstrip(")")
    segments = [tuple(int(x) for x in part.split("..")) for part in location.split(",")]
    if strand == "-":
        segments = segments[::-1]
    bases = "".join(genome_seq[s - 1:e] if strand != "-" else reverse_complement(genome_seq[s - 1:e]) for s, e in segments)
    polyA = next((int(f.split("=")[1]) for f in fields if f.startswith("polyA_bases=")), 0)
    return bases + ("" if stop_stripped else "A" * polyA)


@pytest.mark.parametrize("example", EXPECTED_TAGS)
@pytest.mark.parametrize("protein_oriented", [True, False])
@pytest.mark.parametrize("strip_stop", [False, True])
def test_CDS_headers_locate_their_sequence(tmp_path, example, protein_oriented, strip_stop):
    genome = Genome(name=example, genome_file_path=str(TEST_DATA_DIR / f"input/fasta/{example}.fasta"), quiet=True)
    annot = Annotation(name=example, annot_file_path=str(TEST_DATA_DIR / f"input/annotation/{example}.gff3"), genome=genome, quiet=True)
    annot.export.CDSs(only_main=False, only_cds_main=False, output_dir=str(tmp_path), filename="cds.fasta", protein_oriented=protein_oriented, strip_stop=strip_stop, quiet=True)
    records = [r.split("\n", 1) for r in (tmp_path / "cds.fasta").read_text().split(">")[1:]]
    assert records
    for header, seq in records:
        chrom = header.split("|")[-2]
        assert _rebuild(genome.scaffolds[chrom].seq, header, strip_stop).upper() == seq.replace("\n", "").upper(), header
