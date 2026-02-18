import pytest
import pathlib
from aegis.annotation import Annotation
from aegis.genome import Genome

TEST_FOLDER = pathlib.Path(__file__).parent
DATA_FOLDER = TEST_FOLDER / "test_data"

def test_with_gff_local_file():
    input_path = DATA_FOLDER / "arabidopsis_araport11.gff3"

    ara = Annotation(str(input_path), "ara")

    num_genes = len(ara.all_gene_ids)

    assert num_genes == 3000

def test_with_fasta_local_file():
    input_path = DATA_FOLDER / "arabidopsis_tair10.fasta"

    g = Genome("ara_genome", str(input_path), quiet=True)

    num_chr = len(g.scaffolds)

    assert num_chr == 1
    assert g.scaffolds["Chr4"].size == 18585056