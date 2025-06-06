#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Bgee annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import bgee
from pyBiodatafuse.annotators.bgee import get_gene_expression, get_version_bgee
from pyBiodatafuse.constants import ANATOMICAL_ENTITIES_LIST, BGEE_GENE_EXPRESSION_LEVELS_COL

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
        for _ in range(len(ANATOMICAL_ENTITIES_LIST)):
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
                        "anatomical_entity_id": "UBERON:0000178",
                        "anatomical_entity_name": "blood",
                        "expression_level": 67.35328,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002371",
                        "anatomical_entity_name": "bone marrow",
                        "expression_level": 71.3253,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000955",
                        "anatomical_entity_name": "brain",
                        "expression_level": 88.21244,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000310",
                        "anatomical_entity_name": "breast",
                        "expression_level": 85.39745,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0004535",
                        "anatomical_entity_name": "cardiovascular system",
                        "expression_level": 84.48897999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001007",
                        "anatomical_entity_name": "digestive system",
                        "expression_level": 85.71852,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000948",
                        "anatomical_entity_name": "heart",
                        "expression_level": 84.08238,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0005057",
                        "anatomical_entity_name": "immune organ",
                        "expression_level": 80.20202,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002113",
                        "anatomical_entity_name": "kidney",
                        "expression_level": 91.80176,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002107",
                        "anatomical_entity_name": "liver",
                        "expression_level": 84.28287,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002048",
                        "anatomical_entity_name": "lung",
                        "expression_level": 90.17769,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001016",
                        "anatomical_entity_name": "nervous system",
                        "expression_level": 87.82321,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001264",
                        "anatomical_entity_name": "pancreas",
                        "expression_level": 89.98209,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001987",
                        "anatomical_entity_name": "placenta",
                        "expression_level": 76.73503999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000990",
                        "anatomical_entity_name": "reproductive system",
                        "expression_level": 85.26,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001004",
                        "anatomical_entity_name": "respiratory system",
                        "expression_level": 87.16410999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001434",
                        "anatomical_entity_name": "skeletal system",
                        "expression_level": 70.55227,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                ],
                [
                    {
                        "anatomical_entity_id": "UBERON:0002371",
                        "anatomical_entity_name": "bone marrow",
                        "expression_level": 93.19957,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000955",
                        "anatomical_entity_name": "brain",
                        "expression_level": 77.72678,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000310",
                        "anatomical_entity_name": "breast",
                        "expression_level": 86.06812999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0004535",
                        "anatomical_entity_name": "cardiovascular system",
                        "expression_level": 82.83533,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001007",
                        "anatomical_entity_name": "digestive system",
                        "expression_level": 84.02469,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000948",
                        "anatomical_entity_name": "heart",
                        "expression_level": 81.59009,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0005057",
                        "anatomical_entity_name": "immune organ",
                        "expression_level": 82.4269,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002113",
                        "anatomical_entity_name": "kidney",
                        "expression_level": 82.23363999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002107",
                        "anatomical_entity_name": "liver",
                        "expression_level": 85.3249,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0002048",
                        "anatomical_entity_name": "lung",
                        "expression_level": 85.92834999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001016",
                        "anatomical_entity_name": "nervous system",
                        "expression_level": 78.75558999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001264",
                        "anatomical_entity_name": "pancreas",
                        "expression_level": 86.74894,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001987",
                        "anatomical_entity_name": "placenta",
                        "expression_level": 87.91509,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000990",
                        "anatomical_entity_name": "reproductive system",
                        "expression_level": 84.27241999999998,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001004",
                        "anatomical_entity_name": "respiratory system",
                        "expression_level": 86.36332,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0001434",
                        "anatomical_entity_name": "skeletal system",
                        "expression_level": 91.51794,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                    {
                        "anatomical_entity_id": "UBERON:0000178",
                        "anatomical_entity_name": "blood",
                        "expression_level": 91.3264,
                        "confidence_level_id": "CIO:0000029",
                        "confidence_level_name": "high confidence level",
                        "developmental_stage_id": "UBERON:0000104",
                        "developmental_stage_name": "life cycle",
                    },
                ],
            ]
        )
        expected_data.name = BGEE_GENE_EXPRESSION_LEVELS_COL

        pd.testing.assert_series_equal(
            obtained_data[BGEE_GENE_EXPRESSION_LEVELS_COL], expected_data
        )
        self.assertIsInstance(metadata, dict)
