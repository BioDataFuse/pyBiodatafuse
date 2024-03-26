#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Bgee annotator."""

import json
import os
import unittest
from unittest.mock import patch, Mock
import pandas as pd

from pyBiodatafuse.annotators import bgee
from pyBiodatafuse.annotators.bgee import get_gene_expression, get_version_bgee
from pyBiodatafuse.constants import BGEE, ANATOMICAL_ENTITIES_LIST

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestBgee(unittest.TestCase):
    """Test the Bgee class."""

    @patch("pyBiodatafuse.annotators.bgee.SPARQLWrapper.queryAndConvert")
    def test_get_version_bgee(self, mock_sparql_version):
        """Test that the SPARQL endpoint returns the expected get_version_bgee results."""
        mock_sparql_version.side_effect = [
            {
                "head": {"link": [], "vars": ["date_modified"]},
                "results": {
                    "distinct": False,
                    "ordered": True,
                    "bindings": [{"date_modified": {"type": "literal", "value": "2023-11-01"}}],
                },
            }
        ]

        obtained_version = get_version_bgee()

        expected_version = {"source_version": "2023-11-01"}

        assert obtained_version == expected_version

    @patch("pyBiodatafuse.annotators.bgee.SPARQLWrapper.queryAndConvert")
    def test_get_gene_expression(self, mock_sparql_request):
        """Test the get_gene_expression function."""
        with open(os.path.join(data_file_folder, "bgee_mock_data.json")) as f:
            mock_data = json.load(f)

        mock_data_list = []
        for i in range(len(ANATOMICAL_ENTITIES_LIST)):
            for json_response in mock_data:
                mock_data_list.append(pd.DataFrame.from_dict(json_response))

        mock_sparql_request.side_effect = mock_data_list
        bgee.get_version_bgee = Mock(
            return_value={"source_version": "2023-11-01"}
        )  # Mock the version call
        bgee.check_endpoint_bgee = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["AGRN", "ATXN7"],
                "identifier.source": ["HGNC", "HGNC"],
                "target": ["ENSG00000188157", "ENSG00000163635"],
                "target.source": ["Ensembl", "Ensembl"],
            }
        )

        obtained_data, metadata = get_gene_expression(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "anatomical_entity_id": "UBERON_0000178",
                        "anatomical_entity_name": "blood",
                        "developmental_stage_id": "UBERON_0000104",
                        "developmental_stage_name": "life cycle",
                        "expression_level": 67.3533,
                        "confidence_level_id": "CIO_0000029",
                        "confidence_level_name": "high confidence level",
                    }
                ],
                [
                    {
                        "anatomical_entity_id": "UBERON_0000178",
                        "anatomical_entity_name": "blood",
                        "developmental_stage_id": "UBERON_0000104",
                        "developmental_stage_name": "life cycle",
                        "expression_level": 91.3264,
                        "confidence_level_id": "CIO_0000029",
                        "confidence_level_name": "high confidence level",
                    }
                ],
            ]
        )
        expected_data.name = BGEE

        pd.testing.assert_series_equal(obtained_data[BGEE], expected_data)
        self.assertIsInstance(metadata, dict)
