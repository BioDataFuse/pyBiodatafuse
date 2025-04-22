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
                ]
            )
        )

        stringdb.get_ppi_data = Mock(
            return_value=pd.DataFrame(
                [
                    {
                        "stringId_A": "9606.ENSP00000262374",
                        "stringId_B": "9606.ENSP00000332247",
                        "preferredName_A": "ALG1",
                        "preferredName_B": "ATP6V0A2",
                        "ncbiTaxonId": "9606",
                        "score": 0.415,
                        "nscore": 0.053,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0.049,
                        "escore": 0.0,
                        "dscore": 0.0,
                        "tscore": 0.402,
                    },
                    {
                        "stringId_A": "9606.ENSP00000262374",
                        "stringId_B": "9606.ENSP00000360124",
                        "preferredName_A": "ALG1",
                        "preferredName_B": "PGM1",
                        "ncbiTaxonId": "9606",
                        "score": 0.677,
                        "nscore": 0.0,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0.196,
                        "escore": 0.085,
                        "dscore": 0.287,
                        "tscore": 0.458,
                    },
                    {
                        "stringId_A": "9606.ENSP00000262374",
                        "stringId_B": "9606.ENSP00000268261",
                        "preferredName_A": "ALG1",
                        "preferredName_B": "PMM2",
                        "ncbiTaxonId": "9606",
                        "score": 0.839,
                        "nscore": 0.053,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0.074,
                        "escore": 0.05,
                        "dscore": 0.286,
                        "tscore": 0.771,
                    },
                    {
                        "stringId_A": "9606.ENSP00000262374",
                        "stringId_B": "9606.ENSP00000380793",
                        "preferredName_A": "ALG1",
                        "preferredName_B": "ALG3",
                        "ncbiTaxonId": "9606",
                        "score": 0.872,
                        "nscore": 0.0,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0.095,
                        "escore": 0.094,
                        "dscore": 0.178,
                        "tscore": 0.833,
                    },
                    {
                        "stringId_A": "9606.ENSP00000262374",
                        "stringId_B": "9606.ENSP00000430236",
                        "preferredName_A": "ALG1",
                        "preferredName_B": "ALG11",
                        "ncbiTaxonId": "9606",
                        "score": 0.936,
                        "nscore": 0.0,
                        "fscore": 0,
                        "pscore": 0,
                        "ascore": 0.049,
                        "escore": 0.092,
                        "dscore": 0.0,
                        "tscore": 0.932,
                    },
                ]
            )
        )

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": [
                    "ENSG00000119523",
                ],
                "identifier.source": [
                    "Ensembl",
                ],
                "target": [
                    "ALG2",
                ],
                "target.source": [
                    "HGNC",
                ],
            }
        )

        obtained_data, metadata = get_ppi(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "stringdb_link_to": "ALG1",
                        "Ensembl": "ENSP00000262374",
                        "score": 0.982,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q9BT22",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "ATP4A",
                        "Ensembl": "ENSP00000262623",
                        "score": 0.731,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "P20648",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "ALG6",
                        "Ensembl": "ENSP00000263440",
                        "score": 0.838,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q9Y672",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "PMM2",
                        "Ensembl": "ENSP00000268261",
                        "score": 0.938,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "O15305",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "MUS81",
                        "Ensembl": "ENSP00000307853",
                        "score": 0.793,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q96NY9",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "ATP6V0A2",
                        "Ensembl": "ENSP00000332247",
                        "score": 0.928,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q9Y487",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "PGM1",
                        "Ensembl": "ENSP00000360124",
                        "score": 0.744,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "P36871",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "ALG3",
                        "Ensembl": "ENSP00000380793",
                        "score": 0.913,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q92685",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "PDCD6IP",
                        "Ensembl": "ENSP00000411825",
                        "score": 0.8,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q8WUM4",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                    {
                        "stringdb_link_to": "ALG11",
                        "Ensembl": "ENSP00000430236",
                        "score": 0.969,
                        "Ensembl_link": "ENSP00000417764",
                        "Uniprot-TrEMBL": "Q2TAA5",
                        "Uniprot-TrEMBL_link": "Q9H553",
                    },
                ]
            ]
        )
        expected_data.name = STRING_PPI_COL
        pd.testing.assert_series_equal(obtained_data[STRING_PPI_COL], expected_data)
        self.assertIsInstance(metadata, dict)
