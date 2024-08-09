#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the DISGENET annotator."""

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
    """Test the DISGENET class."""

    @patch("pyBiodatafuse.annotators.disgenet.requests.post")
    def test_get_version_disgenet(self, mock_sparql_version):
        """Test that the API endpoint returns the expected get_version_bgee results."""
        mock_sparql_version.side_effect = {
            "status": "OK",
            "payload": {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"},
        }

        obtained_version = get_version_disgenet()

        expected_version = {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"}

        assert obtained_version == expected_version

    @patch("pyBiodatafuse.annotators.disgenet.requests.post")
    def test_get_gene_disease(self, mock_post_gene_disease):
        """Test the get_gene_disease function."""
        disgenet.get_version_disgenet = Mock(
            return_value={
                "status": "OK",
                "payload": {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"},
            }
        )  # Mock the version call
        disgenet.check_endpoint_disgenet = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["CHRNG"],
                "identifier.source": ["HGNC"],
                "target": ["1146"],
                "target.source": ["NCBI Gene"],
            }
        )

        with open(os.path.join(data_file_folder, "disgenet_mock_data.json")) as f:
            mock_post_gene_disease.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_gene_disease(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "HPO": "",
                        "NCI": "NCI_C101039",
                        "OMIM": "OMIM_265000, OMIM_100730, OMIM_163950",
                        "MONDO": "MONDO_0009926, MONDO_0020746, MONDO_0017415",
                        "ORDO": "ORDO_2990, ORDO_294060",
                        "ICD10": "",
                        "EFO": "",
                        "DO": "DO_0080110, DO_0081322",
                        "MESH": "MESH_C537377",
                        "UMLS": "UMLS_C0265261",
                        "ICD9CM": "",
                        "disease_name": "Multiple pterygium syndrome",
                        "disease_type": "disease",
                        "disease_umlscui": "C0265261",
                        "score": 1.0,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "HPO": "",
                        "NCI": "NCI_C101038",
                        "OMIM": "OMIM_253290, OMIM_100730, OMIM_100690, OMIM_100720",
                        "MONDO": "MONDO_0009668, MONDO_0017415",
                        "ORDO": "ORDO_294060, ORDO_33108",
                        "ICD10": "",
                        "EFO": "",
                        "DO": "DO_0080110",
                        "MESH": "MESH_C537377",
                        "UMLS": "UMLS_C1854678",
                        "ICD9CM": "",
                        "disease_name": "MULTIPLE PTERYGIUM SYNDROME, LETHAL TYPE",
                        "disease_type": "disease",
                        "disease_umlscui": "C1854678",
                        "score": 0.95,
                        "ei": 1.0,
                        "el": None,
                    },
                ]
            ]
        )
        expected_data.name = DISGENET

        pd.testing.assert_series_equal(obtained_data[DISGENET], expected_data)
        self.assertIsInstance(metadata, dict)
