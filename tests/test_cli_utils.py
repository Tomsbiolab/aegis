import pytest
from aegis.cli.utils import split_callback


def test_split_callback_string_comma():
    assert split_callback("gene,CDS,protein") == ["gene", "CDS", "protein"]
    assert split_callback(" gene , CDS , protein ") == ["gene", "CDS", "protein"]


def test_split_callback_single_string():
    assert split_callback("gene") == ["gene"]
    assert split_callback(" gene ") == ["gene"]


def test_split_callback_empty_and_none():
    assert split_callback("") == []
    assert split_callback(None) == []
    assert split_callback([]) == []


def test_split_callback_list_inputs():
    assert split_callback(["gene", "CDS"]) == ["gene", "CDS"]
    assert split_callback(["gene,CDS", "protein"]) == ["gene", "CDS", "protein"]
    assert split_callback([" gene , CDS ", " protein "]) == ["gene", "CDS", "protein"]
