#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the DISGENET annotator."""

import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
from requests import Response, Session

from pyBiodatafuse.annotators import disgenet
from pyBiodatafuse.constants import DISGENET_DISEASE_COL

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestDisgenet(unittest.TestCase):
    """Test the DISGENET class."""

    @patch.object(Session, "get")
    def test_get_version_disgenet(self, requests_get):
        """Test that the API endpoint returns the expected get_version_disgenet results."""
        successful_response = Mock(Response)
        successful_response.status_code = 406
        successful_response.json.return_value = {
            "status": "OK",
            "payload": {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"},
            "httpStatus": 200,
        }
        requests_get.return_value = successful_response

        obtained_version = disgenet.get_version_disgenet(api_key="test")

        expected_version = {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"}

        assert obtained_version == expected_version

    @patch.object(Session, "get")
    def test_get_gene_disease(self, mock_post_gene_disease):
        """Test the get_gene_disease function."""
        disgenet.get_version_disgenet = Mock(
            return_value={"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"}
        )  # Mock the version call
        disgenet.check_endpoint_disgenet = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["ALG14"],
                "identifier.source": ["HGNC"],
                "target": ["199857"],
                "target.source": ["NCBI Gene"],
            }
        )

        mock_post_gene_disease.return_value.ok = True
        with open(os.path.join(data_file_folder, "disgenet_gene_disease_mock.json")) as f:
            mock_post_gene_disease.return_value.text = f.read()

        obtained_data, metadata = disgenet.get_gene_disease(
            api_key="test", bridgedb_df=bridgedb_dataframe
        )

        expected_data = pd.Series(
            [
                [
                    {
                        "disease_name": "Carbohydrate Deficient Glycoprotein Syndrome",
                        "HPO": None,
                        "NCI": "NCI:C84615",
                        "OMIM": None,
                        "MONDO": "MONDO:0015286",
                        "ORDO": "ORDO:137",
                        "EFO": None,
                        "DO": "DOID:5212",
                        "MESH": "MESH:D018981",
                        "UMLS": "UMLS:C0282577",
                        "disease_type": "disease",
                        "score": 0.6000000000000001,
                        "ei": 1.0,
                        "el": "Limited",
                    },
                    {
                        "disease_name": "CONGEN MYASTHENIA GRAVIS",
                        "HPO": None,
                        "NCI": "NCI:C84647",
                        "OMIM": None,
                        "MONDO": "MONDO:0018940",
                        "ORDO": "ORDO:590",
                        "EFO": None,
                        "DO": "DOID:3635",
                        "MESH": "MESH:D020294",
                        "UMLS": "UMLS:C0751882",
                        "disease_type": "disease",
                        "score": 0.6000000000000001,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "disease_name": "CMSWTA",
                        "HPO": None,
                        "NCI": None,
                        "OMIM": "MIM:616227, MIM:612866",
                        "MONDO": "MONDO:0014542",
                        "ORDO": None,
                        "EFO": None,
                        "DO": "DOID:0110658",
                        "MESH": None,
                        "UMLS": "UMLS:C4015596",
                        "disease_type": "disease",
                        "score": 0.6,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "disease_name": "MEPCA",
                        "HPO": None,
                        "NCI": None,
                        "OMIM": "MIM:619036, MIM:612866",
                        "MONDO": "MONDO:0033619",
                        "ORDO": None,
                        "EFO": None,
                        "DO": None,
                        "MESH": None,
                        "UMLS": "UMLS:C5436652",
                        "disease_type": "disease",
                        "score": 0.5,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "disease_name": "Congenital myasthenic syndromes with glycosylation defect",
                        "HPO": None,
                        "NCI": None,
                        "OMIM": None,
                        "MONDO": None,
                        "ORDO": "ORDO:353327",
                        "EFO": "EFO:0700079",
                        "DO": None,
                        "MESH": None,
                        "UMLS": "UMLS:C5680989",
                        "disease_type": "disease",
                        "score": 0.4,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "disease_name": "IDDEBF",
                        "HPO": None,
                        "NCI": None,
                        "OMIM": "MIM:619031, MIM:612866",
                        "MONDO": "MONDO:0033572",
                        "ORDO": None,
                        "EFO": None,
                        "DO": None,
                        "MESH": None,
                        "UMLS": "UMLS:C5436646",
                        "disease_type": "disease",
                        "score": 0.4,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "disease_name": "Congenital myopathies",
                        "HPO": None,
                        "NCI": None,
                        "OMIM": None,
                        "MONDO": "MONDO:0012138, MONDO:0011246, MONDO:0013177, MONDO:0013178, MONDO:0011486, MONDO:0011688, MONDO:0019952, MONDO:0018276",
                        "ORDO": "ORDO:97245",
                        "EFO": None,
                        "DO": "DOID:0112374, DOID:0110632, DOID:0110633, DOID:0110634, DOID:0110635, DOID:0110637, DOID:0110639, DOID:0110640",
                        "MESH": None,
                        "UMLS": "UMLS:C0270960",
                        "disease_type": "disease",
                        "score": 0.4,
                        "ei": None,
                        "el": None,
                    },
                ]
            ]
        )
        expected_data.name = DISGENET_DISEASE_COL

        pd.testing.assert_series_equal(obtained_data[DISGENET_DISEASE_COL], expected_data)
        self.assertIsInstance(metadata, dict)
