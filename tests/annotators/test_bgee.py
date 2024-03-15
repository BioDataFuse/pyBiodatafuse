#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Bgee annotator."""

import json
import os
from unittest.mock import patch

import pandas as pd
import pytest

from pyBiodatafuse.annotators.bgee import get_gene_expression, get_version_bgee
from pyBiodatafuse.constants import BGEE


def test_sparql_get_version_bgee():
    """Test that the SPARQL endpoint returns the expected get_version_bgee results."""
    data_file_folder = os.path.join(os.path.dirname(__file__), "data")
    bgee_version_file_path = os.path.join(data_file_folder, "bgee_version_data.json")
    expected = pd.read_json(bgee_version_file_path)

    obtained_version = get_version_bgee()

    expected_version = {
        "bgee_version": expected["results"]["bindings"][0]["date_modified"]["value"]
    }

    assert obtained_version == expected_version


@patch("pyBiodatafuse.annotators.bgee.SPARQLWrapper.queryAndConvert")
def test_get_version_bgee(mock_sparql_request):
    """Test the get_version_bgee."""
    # TODO: need to find a fix for getting metadata from stable Bgee SPARQL endpoint
    data_file_folder = os.path.join(os.path.dirname(__file__), "data")
    bgee_version_file_path = os.path.join(data_file_folder, "bgee_version_data.json")

    mock_sparql_request.return_value = pd.read_json(bgee_version_file_path)

    obtained_version = get_version_bgee()

    expected_version = {
        "bgee_version": mock_sparql_request.return_value["results"]["bindings"][0]["date_modified"][
            "value"
        ]
    }

    assert obtained_version == expected_version


@patch("pyBiodatafuse.annotators.bgee.SPARQLWrapper.queryAndConvert")
def test_get_gene_expression(mock_sparql_request, bridgedb_dataframe):
    """Test the get_gene_expression function."""
    data_file_folder = os.path.join(os.path.dirname(__file__), "data")
    data_file_path = os.path.join(data_file_folder, "bgee_mock_data.json")

    mock_data_list = []
    with open(data_file_path) as f:
        mock_data = json.load(f)

    for json_response in mock_data:
        mock_data_list.append(pd.DataFrame.from_dict(json_response))

    bgee_version_file_path = os.path.join(data_file_folder, "bgee_version_data.json")
    mock_version_data = pd.read_json(bgee_version_file_path)

    mock_data_list.append(mock_version_data)

    mock_sparql_request.side_effect = mock_data_list

    obtained_data, metadata = get_gene_expression(bridgedb_dataframe)

    obtained_sorted = pd.DataFrame(obtained_data[BGEE][0])
    obtained_sorted = obtained_sorted.astype({"expression_level": float})

    expected_data = pd.read_csv("bgee_expected_data.tsv", sep="\t")

    pd.testing.assert_series_equal(obtained_sorted, expected_data)


@pytest.fixture(scope="module")
def bridgedb_dataframe():
    """Reusable sample Pandas DataFrame to be used as input for the tests."""
    return pd.DataFrame(
        {
            "identifier": ["AGRN"],
            "identifier.source": ["HGNC"],
            "target": ["ENSG00000188157"],
            "target.source": ["Ensembl"],
        }
    )
