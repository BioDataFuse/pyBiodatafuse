#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Pubchem annotator."""
import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators import pubchem
from pyBiodatafuse.annotators.pubchem import get_protein_compound_screened
from pyBiodatafuse.constants import PUBCHEM_COMPOUND_ASSAYS_COL

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestPubchem(unittest.TestCase):
    """Test the PubChem class."""

    # TODO after Pubchem update
    # @patch("pyBiodatafuse.annotators.pubchem.SPARQLWrapper.queryAndConvert")
    # def test_get_version_bgee(self, mock_sparql_request):
    #     """Test the get_version_bgee function."""
    #     version_data = "{"results": {"bindings": [{"date_modified": {"type": "literal", "value": "2023-11-01"}}]}}"
    #     mock_sparql_request.return_value = json.loads(version_data)

    #     obtained_version = get_version_bgee()

    #     expected_version = {"bgee_version": "2023-11-01"}

    #     self.assertEqual(obtained_version, expected_version)

    @patch("pyBiodatafuse.annotators.pubchem.SPARQLWrapper.queryAndConvert")
    def test_get_protein_compound_screened(self, mock_sparql_request):
        """Test the get_protein_compound_screened."""
        with open(os.path.join(data_file_folder, "pubchem_mock_data.json")) as f:
            mock_data = json.load(f)

        # The mock_data is already in the correct SPARQL response format
        mock_sparql_request.return_value = mock_data
        pubchem.check_endpoint_pubchem = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["EGFR", "HTR3A"],
                "identifier.source": ["HGNC", "HGNC"],
                "target": ["P00533", "P46098"],
                "target.source": ["Uniprot-TrEMBL", "Uniprot-TrEMBL"],
            }
        )

        obtained_data, metadata = get_protein_compound_screened(bridgedb_dataframe)

        expected_data = pd.Series([None, None])
        expected_data.name = PUBCHEM_COMPOUND_ASSAYS_COL

        pd.testing.assert_series_equal(obtained_data[PUBCHEM_COMPOUND_ASSAYS_COL], expected_data)
        self.assertIsInstance(metadata, dict)
