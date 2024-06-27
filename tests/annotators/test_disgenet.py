#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Disgenet annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import disgenet
from pyBiodatafuse.annotators.disgenet import get_gene_disease, get_version_disgenet
from pyBiodatafuse.constants import DISGENET

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestDisgenet(unittest.TestCase):
    """Test the Disgenet class."""

    @patch("pyBiodatafuse.annotators.disgenet.SPARQLWrapper.queryAndConvert")
    def test_get_version_disgenet(self, mock_sparql_version):
        """Test that the SPARQL endpoint returns the expected get_version_bgee results."""
        mock_sparql_version.side_effect = [
            {
                "head": {"link": [], "vars": ["title"]},
                "results": {
                    "distinct": False,
                    "ordered": True,
                    "bindings": [
                        {"title": {"type": "literal", "value": "DisGeNET v7.0 RDF Distribution"}},
                    ],
                },
            }
        ]

        obtained_version = get_version_disgenet()

        expected_version = {"source_version": "DisGeNET v7.0 RDF Distribution"}

        assert obtained_version == expected_version

    @patch("pyBiodatafuse.annotators.bgee.SPARQLWrapper.queryAndConvert")
    def test_get_gene_disease(self, mock_sparql_request):
        """Test the get_gene_disease function."""
        with open(os.path.join(data_file_folder, "disgenet_mock_data.json")) as f:
            mock_data = json.load(f)

        mock_data_list = []
        for json_response in mock_data:
            mock_data_list.append(pd.DataFrame.from_dict(json_response))

        mock_sparql_request.side_effect = mock_data_list  # For the two queries

        disgenet.get_version_disgenet = Mock(
            return_value={"source_version": "DisGeNET v7.0 RDF Distribution"}
        )  # Mock the version call
        disgenet.check_endpoint_disgenet = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["SLC25A1", "ENSG00000005108"],
                "identifier.source": ["HGNC", "Ensembl"],
                "target": ["6576", "221981"],
                "target.source": ["NCBI Gene", "NCBI Gene"],
            }
        )

        obtained_data, metadata = get_gene_disease(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "disease_id": "umls:C0338831",
                        "disease_name": "Manic",
                        "score": 0.3,
                        "evidence_source": "CTD_human",
                    },
                    {
                        "disease_id": "umls:C0005587",
                        "disease_name": "Depression, Bipolar",
                        "score": 0.3,
                        "evidence_source": "CTD_human",
                    },
                    {
                        "disease_id": "umls:C0005586",
                        "disease_name": "Bipolar Disorder",
                        "score": 0.4,
                        "evidence_source": "CTD_human",
                    },
                    {
                        "disease_id": "umls:C0024713",
                        "disease_name": "Manic Disorder",
                        "score": 0.3,
                        "evidence_source": "CTD_human",
                    },
                ],
                [
                    {
                        "disease_id": "umls:C2746066",
                        "disease_name": "Combined D-2- and L-2-hydroxyglutaric aciduria",
                        "score": 0.72,
                        "evidence_source": "UNIPROT",
                    },
                    {
                        "disease_id": "umls:C2746066",
                        "disease_name": "Combined D-2- and L-2-hydroxyglutaric aciduria",
                        "score": 0.72,
                        "evidence_source": "ORPHANET",
                    },
                    {
                        "disease_id": "umls:C2746066",
                        "disease_name": "Combined D-2- and L-2-hydroxyglutaric aciduria",
                        "score": 0.72,
                        "evidence_source": "CTD_human",
                    },
                    {
                        "disease_id": "umls:C2746066",
                        "disease_name": "Combined D-2- and L-2-hydroxyglutaric aciduria",
                        "score": 0.72,
                        "evidence_source": "GENOMICS_ENGLAND",
                    },
                    {
                        "disease_id": "umls:C0751884",
                        "disease_name": "Congenital Myasthenic Syndromes, Presynaptic",
                        "score": 0.3,
                        "evidence_source": "ORPHANET",
                    },
                    {
                        "disease_id": "umls:C4748678",
                        "disease_name": "MYASTHENIC SYNDROME, CONGENITAL, 23, PRESYNAPTIC",
                        "score": 0.6,
                        "evidence_source": "UNIPROT",
                    },
                    {
                        "disease_id": "umls:C4748678",
                        "disease_name": "MYASTHENIC SYNDROME, CONGENITAL, 23, PRESYNAPTIC",
                        "score": 0.6,
                        "evidence_source": "GENOMICS_ENGLAND",
                    },
                ],
            ]
        )
        expected_data.name = DISGENET

        pd.testing.assert_series_equal(obtained_data[DISGENET], expected_data)
        self.assertIsInstance(metadata, dict)
