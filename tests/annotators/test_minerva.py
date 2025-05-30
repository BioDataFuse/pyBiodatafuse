#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the MINERVA annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import minerva
from pyBiodatafuse.annotators.minerva import get_gene_pathways, get_version_minerva
from pyBiodatafuse.constants import MINERVA_PATHWAY_COL

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestMinerva(unittest.TestCase):
    """Test the MINERVA class."""

    @patch("pyBiodatafuse.annotators.minerva.requests.get")
    def test_get_version_minerva(self, mock_requests_get):
        """Test the get_version_minerva."""
        # Mocking the response from the MINERVA API
        mock_requests_get.return_value.json.return_value = {
            "version": "16.4.1",
            "buildDate": "2024-02-20T10:25:34+0000",
        }

        map_url = "https://covid19map.elixir-luxembourg.org/minerva/"
        # Call the function under test
        obtained_version = get_version_minerva(map_url)

        # Assert that the obtained version matches the expected version
        expected_version = {"source_version": "16.4.1"}
        assert obtained_version == expected_version

    def test_get_gene_pathways(self):
        """Test the get_gene_pathways function."""
        with open(os.path.join(data_file_folder, "minerva_components.json")) as f:
            mock_data = json.load(f)

        # Mock the request call in the function
        minerva.get_version_bgee = Mock(return_value={"source_version": "16.4.1"})
        minerva.check_endpoint_minerva = Mock(return_value=True)
        minerva.get_minerva_components = Mock(
            return_value=("https://covid19map.elixir-luxembourg.org/minerva/", mock_data)
        )

        # Mocking the response from the MINERVA API
        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["ABCG2", "AK1", "AHR"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["ENSG00000118777", "ENSG00000106992", "ENSG00000106546"],
                "target.source": ["Ensembl", "Ensembl", "Ensembl"],
            }
        )

        obtained_df, metadata = get_gene_pathways(bridgedb_dataframe, "COVID19 Disease Map")

        # Define the expected DataFrame
        expected_df = pd.Series(
            [
                [
                    {
                        "pathway_id": "MINERVA:952",
                        "pathway_label": "HMOX1 pathway",
                        "pathway_gene_counts": 113,
                    }
                ],
                [
                    {
                        "pathway_id": "MINERVA:953",
                        "pathway_label": "Kynurenine synthesis pathway",
                        "pathway_gene_counts": 45,
                    }
                ],
                [
                    {
                        "pathway_id": "MINERVA:942",
                        "pathway_label": "Nsp14 and metabolism",
                        "pathway_gene_counts": 96,
                    }
                ],
            ]
        )
        expected_df.name = MINERVA_PATHWAY_COL

        pd.testing.assert_series_equal(obtained_df[MINERVA_PATHWAY_COL], expected_df)
        self.assertIsInstance(metadata, dict)
