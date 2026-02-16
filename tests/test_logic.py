import pytest
import pathlib
from aegis.annotation import Annotation
from aegis.tools import pickle_load

# This finds the directory where this script is located
HERE = pathlib.Path(__file__).parent
DATA_FOLDER = HERE / "test_data"

def test_with_local_file():
    # Point to the small version of your big file
    input_path = DATA_FOLDER / "arabidopsis_araport11.gff3"
    input_path_gold = DATA_FOLDER / "arabidopsis_araport11.pkl"

    ara = Annotation(str(input_path), "ara")

    num_genes = len(ara.all_gene_ids)

    ara_gold = pickle_load(str(input_path_gold))

    num_genes_gold = len(ara_gold.all_gene_ids)

    assert num_genes == num_genes_gold