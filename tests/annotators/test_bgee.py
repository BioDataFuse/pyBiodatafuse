#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Bgee annotator."""

import json
import os
from unittest.mock import patch

import pandas as pd
import pytest
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.annotators.bgee import get_gene_expression, get_version_bgee


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

def test_sparql_endpoint_bgee():
    """Test the availability of the Bgee SPARQL endpoint."""
    endpoint = "https://www.bgee.org/sparql/"
    sparql_query = "ASK WHERE {?s ?p ?o}"

    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        assert True
    except SPARQLWrapperException:
        assert False

@patch("pyBiodatafuse.annotators.bgee.SPARQLWrapper.queryAndConvert")
def test_get_gene_expression(mock_sparql_request, bridgedb_dataframe):
    """Test the get_gene_expression function."""
    data_file_folder = os.path.join(os.path.dirname(__file__), "data")
    data_file_path = os.path.join(data_file_folder, "bgee_mock_data.json")
    with open(data_file_path) as f:
        mock_data = json.load(f)

    mock_bgee_data = pd.DataFrame(mock_data)

    bgee_version_file_path = os.path.join(data_file_folder, "bgee_version_data.json")
    mock_version_data = pd.read_json(bgee_version_file_path)

    mocked_data = [mock_bgee_data, mock_version_data]

    mock_sparql_request.side_effect = mocked_data

    obtained_data, metadata = get_gene_expression(bridgedb_dataframe)

    expected_data = pd.read_json(os.path.join(data_file_folder, "bgee_expected_data.json"))
    expected_data = set(expected_data["anatomical_entity_id"].unique())

    obtained_sorted = pd.DataFrame(obtained_data["Bgee"][0])
    obtained_sorted = set(obtained_sorted["anatomical_entity_id"].unique())

    assert obtained_sorted == expected_data


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
