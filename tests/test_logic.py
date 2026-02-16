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

    state = vars(ara)

    ara_gold = pickle_load(str(input_path_gold))

    state_gold = vars(ara_gold)

    assert state == state_gold