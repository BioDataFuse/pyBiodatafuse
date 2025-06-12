#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the WikiPathways annotator."""

import unittest
from unittest.mock import Mock, patch

import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators import wikipathways
from pyBiodatafuse.annotators.wikipathways import get_gene_wikipathways, get_version_wikipathways
from pyBiodatafuse.constants import WIKIPATHWAYS_PATHWAY_COL


class TestWikipathway(unittest.TestCase):
    """Test the WikiPathway class."""

    @patch("pyBiodatafuse.annotators.wikipathways.SPARQLWrapper.queryAndConvert")
    def test_wikipathways_version(self, mock_sparql_version):
        """Test the get_version_wikipathways."""
        mock_sparql_version.side_effect = [
            {
                "head": {"link": [], "vars": ["title"]},
                "results": {
                    "distinct": False,
                    "ordered": True,
                    "bindings": [
                        {"title": {"type": "literal", "value": "WikiPathways RDF 20240310"}}
                    ],
                },
            }
        ]

        obtained_version = get_version_wikipathways()

        expected_version = {"source_version": "WikiPathways RDF 20240310"}

        assert obtained_version == expected_version

    @patch("pyBiodatafuse.annotators.wikipathways.SPARQLWrapper.queryAndConvert")
    def test_get_gene_wikipathways(self, mock_sparql_request):
        """Test the get_gene_wikipathways."""
        mock_sparql_request.side_effect = [
            {
                "head": {
                    "link": [],
                    "vars": ["gene_id", "pathway_id", "pathway_label", "pathway_gene_count"],
                },
                "results": {
                    "distinct": False,
                    "ordered": True,
                    "bindings": [
                        {
                            "gene_id": {"type": "literal", "value": "85365"},
                            "pathway_id": {"type": "literal", "value": "WP5153"},
                            "pathway_label": {
                                "type": "literal",
                                "xml:lang": "en",
                                "value": "N-glycan biosynthesis",
                            },
                            "pathway_gene_count": {
                                "type": "typed-literal",
                                "datatype": "http: //www.w3.org/2001/XMLSchema#integer",
                                "value": "57",
                            },
                        },
                        {
                            "gene_id": {"type": "literal", "value": "199857"},
                            "pathway_id": {"type": "literal", "value": "WP5153"},
                            "pathway_label": {
                                "type": "literal",
                                "xml:lang": "en",
                                "value": "N-glycan biosynthesis",
                            },
                            "pathway_gene_count": {
                                "type": "typed-literal",
                                "datatype": "http://www.w3.org/2001/XMLSchema#integer",
                                "value": "57",
                            },
                        },
                    ],
                },
            }
        ]

        wikipathways.get_version_wikipathways = Mock(
            return_value={"source_version": "WikiPathways RDF 20240310"}
        )  # Mock the version call
        wikipathways.check_endpoint_wikipathways = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["ALG14", "ALG2", "CHRNA1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["199857", "85365", "1134"],
                "target.source": ["NCBI Gene", "NCBI Gene", "NCBI Gene"],
            }
        )

        obtained_data, metadata = get_gene_wikipathways(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "pathway_id": "WP:WP5153",
                        "pathway_label": "N-glycan biosynthesis",
                        "pathway_gene_counts": 57.0,
                    }
                ],
                [
                    {
                        "pathway_id": "WP:WP5153",
                        "pathway_label": "N-glycan biosynthesis",
                        "pathway_gene_counts": 57.0,
                    }
                ],
                [{"pathway_id": nan, "pathway_label": nan, "pathway_gene_counts": nan}],
            ]
        )
        expected_data.name = WIKIPATHWAYS_PATHWAY_COL

        pd.testing.assert_series_equal(obtained_data[WIKIPATHWAYS_PATHWAY_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    # test for molecular pathways missing
