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

        mock_data_list = []
        for json_response in mock_data:
            mock_data_list.append(pd.DataFrame.from_dict(json_response))

        mock_sparql_request.side_effect = mock_data_list
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

        expected_data = pd.Series(
            [
                [
                    {
                        "assay_type": nan,
                        "outcome": nan,
                        "compound_cid": nan,
                        "compound_name": nan,
                        "smiles": nan,
                        "inchi": nan,
                        "pubchem_assay_id": nan,
                    }
                ],
                [
                    {
                        "assay_type": "Ki",
                        "outcome": "active",
                        "compound_cid": "CID9911844",
                        "compound_name": "DR-4485 free base",
                        "smiles": "C1CC2=C(C=CC3=C2C(C1)(C(=O)N3)CCCCN4CCC(=CC4)C5=CC=C(C=C5)Cl)Cl",
                        "inchi": "InChI=1S/C26H28Cl2N2O/c27-20-7-5-18(6-8-20)19-11-16-30(17-12-19)15-2-1-13-26-14-3-4-21-22(28)9-10-23(24(21)26)29-25(26)31/h5-11H,1-4,12-17H2,(H,29,31)",
                        "pubchem_assay_id": "AID6505",
                    }
                ],
            ]
        )
        expected_data.name = PUBCHEM_COMPOUND_ASSAYS_COL

        pd.testing.assert_series_equal(obtained_data[PUBCHEM_COMPOUND_ASSAYS_COL], expected_data)
        self.assertIsInstance(metadata, dict)
