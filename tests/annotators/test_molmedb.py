#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the MolMeDB annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators import molmedb
from pyBiodatafuse.annotators.molmedb import (
    get_compound_gene_inhibitor,
    get_gene_compound_inhibitor,
)
from pyBiodatafuse.constants import MOLMEDB_COMPOUND_PROTEIN_COL, MOLMEDB_PROTEIN_COMPOUND_COL

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestMolMeDb(unittest.TestCase):
    """Test the MolMeDB class."""

    # TODO after MolMeDB update
    # @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
    # def test_get_version_bgee(self, mock_sparql_request):
    #     """Test the get_version_bgee function."""
    #     version_data = '{"results": {"bindings": [{"date_modified": {"type": "literal", "value": "2023-11-01"}}]}}'
    #     mock_sparql_request.return_value = json.loads(version_data)

    #     obtained_version = get_version_bgee()

    #     expected_version = {"bgee_version": "2023-11-01"}

    #     self.assertEqual(obtained_version, expected_version)

    @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
    def test_get_gene_compound_inhibitor(self, mock_sparql_request):
        """Test the get_gene_compound_inhibitor."""
        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "identifier": ["SLC17A9", "SLC17A9", "SLC25A1", "SLC25A1", "KCNJ5", "KCNJ5"],
                "identifier.source": ["HGNC", "HGNC", "HGNC", "HGNC", "HGNC", "HGNC"],
                "target": ["Q5W197", "Q9BYT1", "P53007", "B4DP62", "P48544", "A0A5J6E2W8"],
                "target.source": [
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                ],
            }
        )

        with open(os.path.join(data_file_folder, "molmedb_mock_data.json")) as f:
            mock_data = json.load(f)

        mock_sparql_request.side_effect = [mock_data]

        molmedb.check_endpoint_molmedb = Mock(return_value=True)

        obtained_data, metadata = get_gene_compound_inhibitor(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "MolMeDB_compound_name": "Euphorbiaproliferin C",
                        "MolMeDB_inchikey": "MEMULCZBXUZFOZ-UHFFFAOYSA-N",
                        "MolMeDB_smiles": "CC(=O)OC1C2(C)OCC3(C(=O)C=CC(C(C)(C)OC(C)=O)C23)C(OC(C)=O)C2C(OC(=O)C(C)C)C(C)CC21OC(C)=O",
                        "MolMeDB_id": "MM470852",
                        "source_pmid": "27441737",
                    },
                    {
                        "MolMeDB_compound_name": "Euphornin",
                        "MolMeDB_inchikey": "BRVXVMOWTHQKHC-LVYIKVSWSA-N",
                        "MolMeDB_smiles": "CC(=O)OC1CC(OC(C)=O)C(C)(C)/C=C\\C(C)C(OC(C)=O)C2(O)CC(C)C(OC(=O)c3ccccc3)C2/C=C\\1C",
                        "MolMeDB_id": "MM470853",
                        "source_pmid": "30411614",
                    },
                ],
                [
                    {
                        "MolMeDB_compound_name": "Euphorbiaproliferin C",
                        "MolMeDB_inchikey": "MEMULCZBXUZFOZ-UHFFFAOYSA-N",
                        "MolMeDB_smiles": "CC(=O)OC1C2(C)OCC3(C(=O)C=CC(C(C)(C)OC(C)=O)C23)C(OC(C)=O)C2C(OC(=O)C(C)C)C(C)CC21OC(C)=O",
                        "MolMeDB_id": "MM470852",
                        "source_pmid": "27441737",
                    },
                    {
                        "MolMeDB_compound_name": "Euphornin",
                        "MolMeDB_inchikey": "BRVXVMOWTHQKHC-LVYIKVSWSA-N",
                        "MolMeDB_smiles": "CC(=O)OC1CC(OC(C)=O)C(C)(C)/C=C\\C(C)C(OC(C)=O)C2(O)CC(C)C(OC(=O)c3ccccc3)C2/C=C\\1C",
                        "MolMeDB_id": "MM470853",
                        "source_pmid": "30411614",
                    },
                ],
                [
                    {
                        "MolMeDB_compound_name": "MM17483",
                        "MolMeDB_inchikey": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                        "MolMeDB_smiles": "O=P(O)(O)C(Cl)(Cl)P(=O)(O)O",
                        "MolMeDB_id": "MM17483",
                        "source_pmid": "28720702",
                    }
                ],
                [
                    {
                        "MolMeDB_compound_name": "MM17483",
                        "MolMeDB_inchikey": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                        "MolMeDB_smiles": "O=P(O)(O)C(Cl)(Cl)P(=O)(O)O",
                        "MolMeDB_id": "MM17483",
                        "source_pmid": "28720702",
                    }
                ],
                [
                    {
                        "MolMeDB_compound_name": nan,
                        "MolMeDB_inchikey": nan,
                        "MolMeDB_smiles": nan,
                        "MolMeDB_id": nan,
                        "source_pmid": nan,
                    }
                ],
                [
                    {
                        "MolMeDB_compound_name": nan,
                        "MolMeDB_inchikey": nan,
                        "MolMeDB_smiles": nan,
                        "MolMeDB_id": nan,
                        "source_pmid": nan,
                    }
                ],
            ]
        )
        expected_data.name = MOLMEDB_PROTEIN_COMPOUND_COL

        pd.testing.assert_series_equal(obtained_data[MOLMEDB_PROTEIN_COMPOUND_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
    def test_get_compound_gene_inhibitor(self, mock_sparql_request):
        """Test the get_compound_gene_inhibitor."""
        with open(os.path.join(data_file_folder, "molmedb_compound_gene_mock.json")) as f:
            mock_data = json.load(f)

        mock_sparql_request.side_effect = [mock_data]

        molmedb.check_endpoint_molmedb = Mock(return_value=True)

        bridgedb_dataframe_compounds = pd.DataFrame(
            {
                "identifier": ["10041551", "10025195", "2153"],
                "identifier.source": ["PubChem Compound", "PubChem Compound", "PubChem Compound"],
                "target": [
                    "OVVBIIBBRZVPAL-UHFFFAOYSA-N",
                    "LEJRLSZVESQKJK-UHFFFAOYSA-N",
                    "ZFXYFBGIUFBOJW-UHFFFAOYSA-N",
                ],
                "target.source": ["InChIKey", "InChIKey", "InChIKey"],
            }
        )

        obtained_data, metadata = get_compound_gene_inhibitor(bridgedb_dataframe_compounds)

        expected_data = pd.Series(
            [
                [
                    {
                        "MolMeDB_uniprot_trembl_id": nan,
                        "MolMeDB_hgnc_symbol": nan,
                        "source_pmid": nan,
                    }
                ],
                [
                    {
                        "MolMeDB_uniprot_trembl_id": "Uniprot-TrEMBL:P23975",
                        "MolMeDB_hgnc_symbol": "SLC6A2",
                        "source_pmid": "20223878",
                    },
                    {
                        "MolMeDB_uniprot_trembl_id": "Uniprot-TrEMBL:P31645",
                        "MolMeDB_hgnc_symbol": "SLC6A4",
                        "source_pmid": "20223878",
                    },
                ],
                [
                    {
                        "MolMeDB_uniprot_trembl_id": nan,
                        "MolMeDB_hgnc_symbol": nan,
                        "source_pmid": nan,
                    }
                ],
            ]
        )
        expected_data.name = MOLMEDB_COMPOUND_PROTEIN_COL

        pd.testing.assert_series_equal(obtained_data[MOLMEDB_COMPOUND_PROTEIN_COL], expected_data)

        self.assertIsInstance(metadata, dict)
