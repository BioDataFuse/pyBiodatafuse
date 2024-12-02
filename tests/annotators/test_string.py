#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Stringdb annotator."""

import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import stringdb
from pyBiodatafuse.annotators.stringdb import get_ppi, get_version_stringdb
from pyBiodatafuse.constants import STRING_PPI_COL


class TestString(unittest.TestCase):
    """Test the String class."""

    @patch("pyBiodatafuse.annotators.stringdb.requests.get")
    def test_get_version_stringdb(self, mock_requests_get):
        """Test the get_version_stringdb."""
        mock_requests_get.return_value.json.return_value = [
            {"string_version": "12.0", "stable_address": "https://version-12-0.string-db.org"}
        ]

        obtained_version = get_version_stringdb()

        expected_version = {"source_version": "12.0"}

        assert obtained_version == expected_version

    def test_get_ppi(self):
        """Test the get_ppi function."""
        stringdb.check_endpoint_stringdb = Mock(return_value=True)
        stringdb.get_version_stringdb = Mock(return_value={"source_version": "12.0"})
        stringdb.get_string_ids = Mock(
            return_value=pd.DataFrame(
                [
                    {
                        "queryIndex": 0,
                        "queryItem": "ENSG00000119523",
                        "stringId": "9606.ENSP00000417764",
                        "ncbiTaxonId": 9606,
                        "taxonName": "Homo sapiens",
                        "preferredName": "ALG2",
                    },
                    {
                        "queryIndex": 1,
                        "queryItem": "ENSG00000138435",
                        "stringId": "9606.ENSP00000261007",
                        "ncbiTaxonId": 9606,
                        "taxonName": "Homo sapiens",
                        "preferredName": "CHRNA1",
                    },
                    {
                        "queryIndex": 2,
                        "queryItem": "ENSG00000172339",
                        "stringId": "9606.ENSP00000359224",
                        "ncbiTaxonId": 9606,
                        "taxonName": "Homo sapiens",
                        "preferredName": "ALG14",
                    },
                ]
            )
        )

        stringdb.get_ppi_data = Mock(
            return_value=pd.DataFrame(
                [
                    {
                        "stringId_A": "9606.ENSP00000261007",
                        "stringId_B": "9606.ENSP00000359224",
                        "preferredName_A": "CHRNA1",
                        "preferredName_B": "ALG14",
                        "ncbiTaxonId": "9606",
                        "score": 0.543,
                        "nscore": 0,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0,
                        "escore": 0,
                        "dscore": 0,
                        "tscore": 0.543,
                    },
                    {
                        "stringId_A": "9606.ENSP00000359224",
                        "stringId_B": "9606.ENSP00000417764",
                        "preferredName_A": "ALG14",
                        "preferredName_B": "ALG2",
                        "ncbiTaxonId": "9606",
                        "score": 0.633,
                        "nscore": 0,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0.067,
                        "escore": 0,
                        "dscore": 0.119,
                        "tscore": 0.589,
                    },
                ]
            )
        )

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["ALG14", "ALG2", "CHRNA1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["ENSG00000172339", "ENSG00000119523", "ENSG00000138435"],
                "target.source": ["Ensembl", "Ensembl", "Ensembl"],
            }
        )

        obtained_data, metadata = get_ppi(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "stringdb_link_to": "CHRNA1",
                        "Ensembl": "Ensembl:ENSP00000261007",
                        "score": 0.543,
                        "Uniprot-TrEMBL": "ALG14",
                    },
                    {
                        "stringdb_link_to": "ALG2",
                        "Ensembl": "Ensembl:ENSP00000417764",
                        "score": 0.633,
                        "Uniprot-TrEMBL": "ALG14",
                    },
                ],
                [
                    {
                        "stringdb_link_to": "ALG14",
                        "Ensembl": "Ensembl:ENSP00000359224",
                        "score": 0.633,
                        "Uniprot-TrEMBL": "ALG2",
                    }
                ],
                [
                    {
                        "stringdb_link_to": "ALG14",
                        "Ensembl": "Ensembl:ENSP00000359224",
                        "score": 0.543,
                        "Uniprot-TrEMBL": "CHRNA1",
                    }
                ],
            ]
        )
        expected_data.name = STRING_PPI_COL

        pd.testing.assert_series_equal(obtained_data[STRING_PPI_COL], expected_data)
        self.assertIsInstance(metadata, dict)
