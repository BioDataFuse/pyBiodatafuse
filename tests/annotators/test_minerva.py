#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the MINERVA annotator."""

from unittest.mock import Mock, patch
import pandas as pd
import pytest

from pyBiodatafuse.annotators.minerva import get_version_minerva, get_gene_minerva_pathways, get_minerva_components, list_projects


@patch("pyBiodatafuse.annotators.minerva.requests.get")
def test_get_version_minerva(mock_requests_get):
    # Mocking the response from the MINERVA API
    mock_response = Mock()
    mock_requests_get.return_value.json.return_value = "16.4.1"

    # Call the function under test
    obtained_version = get_version_minerva()

    # Assert that the obtained version matches the expected version
    expected_version = "16.4.1"
    assert obtained_version == expected_version


@pytest.fixture
def bridgedb_dataframe():
    """Reusable sample Pandas DataFrame to be used as input for the tests."""
    return pd.DataFrame(
        {
            "identifier": ["ABCG2", "AK1", "AHR"],
            "identifier.source": ["HGNC", "HGNC", "HGNC"],
            "target": ["9429", "203", "196"],
            "target.source": ["NCBI Gene", "NCBI Gene", "NCBI Gene"],
        }
    )


@patch("pyBiodatafuse.annotators.minerva.requests.get")
def test_get_gene_minerva_pathways(bridgedb_dataframe):
    # Mocking the response from the MINERVA API

    # Call the function under test (test for COVID19 disease map)
    project_list_df = list_projects()
    map_components = get_minerva_components(project_list_df, map_name='COVID19 Disease Map', get_reactions=False )
    obtained_df = get_gene_minerva_pathways(bridgedb_dataframe, map_components)
    
    obtained_df=obtained_df.sort_values(by='identifier')
    
    
    # Define the expected DataFrame
    expected_df = pd.DataFrame({
        "identifier": ["ABCG2", "AK1","AHR"],
        "identifier.source": ["HGNC", "HGNC", "HGNC"],
        "target": ["9429", "203", "196"],
        "target.source": ["NCBI Gene", "NCBI Gene", "NCBI Gene"],
        "MINERVA": ["[{'pathwayId': 952, 'pathwayLabel': 'HMOX1 pathway', 'pathwayGeneCount': 113}]",  "[{'pathwayId': 942, 'pathwayLabel': 'Nsp14 and metabolism', 'pathwayGeneCount': 96}]","[{'pathwayId': 953, 'pathwayLabel': 'Kynurenine synthesis pathway', 'pathwayGeneCount': 45}]"],
    
    })
    import ast
    expected_df['MINERVA'] = expected_df['MINERVA'].apply(ast.literal_eval)
    expected_df=expected_df.sort_values(by='identifier')

    # Assert that the obtained DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(obtained_df, expected_df)


# Run the tests
if __name__ == "__main__":
    pytest.main()
