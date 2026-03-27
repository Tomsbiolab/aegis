"""
Shared Pytest fixtures for the Aegis test suite.
"""

import pytest
from pathlib import Path

from aegis.utils.gtf_gff import parse_gff_parts
from aegis.feature import Feature
from aegis.gene import Gene
from aegis.transcript import Transcript
from aegis.subfeatures import CDS, Exon, UTR, Intron
from aegis.conf import default_features, default_introns


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA_DIR = REPO_ROOT / "tests" / "test_data"


@pytest.fixture
def test_data_dir():
    """Return the path to the real test_data/ directory."""
    return TEST_DATA_DIR


# ---------------------------------------------------------------------------
# GFF3 fixtures — loaded from test_data/
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_gff3_file():
    """Return path to the minimal GFF3 test file."""
    return str(TEST_DATA_DIR / "input/annotation/minimal.gff3")


@pytest.fixture
def multi_gene_gff3_file():
    """Return path to the multi-gene GFF3 test file."""
    return str(TEST_DATA_DIR / "input/annotation/multi_gene.gff3")


@pytest.fixture
def rich_gff3_file():
    """Return path to the rich GFF3 test file (3 genes, 2 chromosomes, alt transcripts)."""
    return str(TEST_DATA_DIR / "input/annotation/rich.gff3")


# ---------------------------------------------------------------------------
# Edge cases GFF3 files
# ---------------------------------------------------------------------------

@pytest.fixture
def exon_only_gff3_file():
    """Gene with only exons (no CDS/UTR) — non-coding mRNA."""
    return str(TEST_DATA_DIR / "input/annotation/exon_only.gff3")


@pytest.fixture
def cds_only_gff3_file():
    """Gene with only CDS segments (no explicit exons or UTRs)."""
    return str(TEST_DATA_DIR / "input/annotation/cds_only.gff3")


@pytest.fixture
def no_subfeatures_gff3_file():
    """Gene with mRNA but zero subfeatures."""
    return str(TEST_DATA_DIR / "input/annotation/no_subfeatures.gff3")


@pytest.fixture
def noncoding_transcripts_gff3_file():
    """Several non-coding transcript types (lnc_RNA, tRNA, rRNA, snoRNA)."""
    return str(TEST_DATA_DIR / "input/annotation/noncoding_transcripts.gff3")


@pytest.fixture
def pseudogene_gff3_file():
    """Pseudogene with exons pointing directly at the gene (no transcript)."""
    return str(TEST_DATA_DIR / "input/annotation/pseudogene.gff3")


@pytest.fixture
def exons_to_collapse_gff3_file():
    """Gene with overlapping exons that need collapsing."""
    return str(TEST_DATA_DIR / "input/annotation/exons_to_collapse.gff3")


@pytest.fixture
def multi_cds_ids_gff3_file():
    """Transcript with CDS segments having different IDs (polycistronic)."""
    return str(TEST_DATA_DIR / "input/annotation/multi_cds_ids.gff3")


@pytest.fixture
def transcript_no_parent_gff3_file():
    """Transcript without a Parent attribute (gene inferred)."""
    return str(TEST_DATA_DIR / "input/annotation/transcript_no_parent.gff3")


@pytest.fixture
def cds_no_parent_gff3_file():
    """CDS subfeatures with no Parent attribute (gene+transcript inferred)."""
    return str(TEST_DATA_DIR / "input/annotation/cds_no_parent.gff3")


@pytest.fixture
def multiple_isoforms_gff3_file():
    """Gene with 3 mRNA isoforms sharing exon regions but different CDS."""
    return str(TEST_DATA_DIR / "input/annotation/multiple_isoforms.gff3")


@pytest.fixture
def subfeature_parent_is_gene_gff3_file():
    """Exons and CDSs reference a gene ID as Parent (no mRNA/transcript line)."""
    return str(TEST_DATA_DIR / "input/annotation/subfeature_parent_is_gene.gff3")


@pytest.fixture
def geneID_attribute_as_parent_gff3_file():
    """GFF3 without gene features and with geneID rather than Parent."""
    return str(TEST_DATA_DIR / "input/annotation/geneID_attribute_as_parent.gff3")

@pytest.fixture
def shared_parents_gff3_file():
    """GFF3 with exons that should share parents but are neither named uniformly nor sharing any parents."""
    return str(TEST_DATA_DIR / "input/annotation/shared_parents_plus_negative_strand.gff3")


@pytest.fixture
def clash_of_ids_gff3_file():
    """ID clashes at the gene and transcript level of different styles."""
    return str(TEST_DATA_DIR / "input/annotation/clash_of_ids.gff3")

@pytest.fixture
def miRNA_human_format_gff3_file():
    """miRNA GFF3 file in human format."""
    return str(TEST_DATA_DIR / "input/annotation/miRNA_human_format.gff3")

@pytest.fixture
def miRNA_arabidopsis_format_gff3_file():
    """miRNA GFF3 file in Arabidopsis format."""
    return str(TEST_DATA_DIR / "input/annotation/miRNA_arabidopsis_format.gff3")

@pytest.fixture
def merge_gff3_file_1():
    """GFF3 file with features that will merge differently depending on the parameters."""
    return str(TEST_DATA_DIR / "input/annotation/for_merge_1.gff3")

@pytest.fixture
def merge_gff3_file_2():
    """GFF3 file with features that will merge differently depending on the parameters."""
    return str(TEST_DATA_DIR / "input/annotation/for_merge_2.gff3")

@pytest.fixture
def self_overlapping_genes_gff3_file():
    """GFF3 file with self overlapping genes."""
    return str(TEST_DATA_DIR / "input/annotation/self_overlapping_genes.gff3")

@pytest.fixture
def other_overlapping_genes_gff3_file_1():
    """GFF3 file 1 with other overlapping genes."""
    return str(TEST_DATA_DIR / "input/annotation/other_overlapping_genes_1.gff3")

@pytest.fixture
def other_overlapping_genes_gff3_file_2():
    """GFF3 file 2 with other overlapping genes."""
    return str(TEST_DATA_DIR / "input/annotation/other_overlapping_genes_2.gff3")

# ---------------------------------------------------------------------------
# FASTA fixture — loaded from test_data/
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_fasta_file():
    """Return path to the minimal FASTA test file."""
    return str(TEST_DATA_DIR / "input/fasta/minimal.fasta")


# ---------------------------------------------------------------------------
# Ready-made class instances and creation helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def make_feature():
    def _make(feature_id="feat001", ch="chr1", source="aegis", feature="gene", strand="+", start=100, end=500, score=".", parents=None, attributes=None) -> Feature:
        if parents is None: parents = []
        if attributes is None: attributes = {"Name": "TestFeature", "Alias": ["TF1", "TF2"], "Symbol": "TFS"}
        return Feature(feature_id, ch, source, feature, strand, start, end, score, parents[:], attributes.copy())
    return _make

@pytest.fixture
def make_gene():
    def _make(pseudogene=False, transposable=False, feature_id="gene001", ch="chr1", source="aegis", feature="gene", strand="+", start=100, end=2000, score=".", parents=None, attributes=None) -> Gene:
        if parents is None: parents = []
        if attributes is None: attributes = {"Name": "TestGene"}
        return Gene(pseudogene, transposable, feature_id, ch, source, feature, strand, start, end, score, parents[:], attributes.copy())
    return _make

@pytest.fixture
def make_transcript():
    def _make(feature_id="mRNA1", ch="chr1", source="aegis", feature="mRNA", strand="+", start=1000, end=5000, score=".", parents=None, attributes=None) -> Transcript:
        if parents is None: parents = ["gene1"]
        if attributes is None: attributes = {}
        return Transcript(feature_id, ch, source, feature, strand, start, end, score, parents[:], attributes.copy())
    return _make

@pytest.fixture
def make_CDS_segment():
    """Create a Feature segment to use as a CDS segment."""
    def _make(feature_id="seg1", ch="chr1", source="aegis", feature="mRNA", strand="+", start=100, end=300, score=".", parents=None, attributes=None):
        if parents is None: parents = ["mRNA1"]
        if attributes is None: attributes = {}
        return Feature(feature_id, ch, source, feature, strand, start, end, score, parents[:], attributes.copy())
    return _make

@pytest.fixture
def make_CDS(make_CDS_segment):
    def _make(segments=None, feature_id="cds1", source="aegis", ch="chr1", strand="+", score=".", parents=None, attributes=None, feature="CDS"):
        if parents is None: parents = ["mRNA1"]
        if attributes is None: attributes = {}
        if segments is None:
            segments = [
                make_CDS_segment("seg1", start=1200, end=2000), 
                make_CDS_segment("seg2", start=3000, end=4500),
            ]
        first = segments[0]
        last = segments[-1]
        return CDS(CDS_segments=segments, feature_id=feature_id, ch=ch, source=source, feature=feature, strand=strand, start=first.start, end=last.end, score=score, parents=parents[:], attributes=attributes.copy())
    return _make

@pytest.fixture
def make_exon():
    def _make(feature_id="e1", start=1000, end=2000, strand="+", parents=None, score=".", source="aegis", ch="chr1", feature="exon", attributes=None):
        if parents is None: parents = ["mRNA1"]
        if attributes is None: attributes = {}
        return Exon(feature_id=feature_id, ch=ch, source=source, feature=feature, strand=strand, start=start, end=end, score=score, parents=parents[:], attributes=attributes.copy())
    return _make

@pytest.fixture
def create_test_feature():
    def _make(line: str) -> Feature|Gene|Transcript|Intron|UTR|Exon:
        parts = line.strip().split("\t")
        if len(parts) != 9:
            raise ValueError(f"GFF line must have 9 columns, got {len(parts)}")
        
        entry = parse_gff_parts(parts)
        
        args = {
            "feature_id": entry.id,
            "ch": entry.ch,
            "source": entry.source,
            "feature": entry.feature,
            "strand": entry.strand,
            "start": entry.start,
            "end": entry.end,
            "score": entry.score,
            "parents": entry.parents,
            "attributes": entry.attributes
        }

        if entry.feature in default_features["gene"]:
            return Gene(pseudogene=entry.pseudogene, transposable=entry.transposable, **args)
        elif entry.feature in default_features["transcript"]:
            return Transcript(**args)
        elif entry.feature in default_features["exon"]:
            return Exon(**args)
        elif entry.feature in default_features["UTR"]:
            return UTR(**args)
        elif entry.feature in default_introns:
            return Intron(**args)
        else:
            return Feature(**args)
    return _make