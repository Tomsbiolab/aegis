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
        "YP_003024026.1": {"transl_except": "(pos:955..956,aa:TERM)", "table": "2"},
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
