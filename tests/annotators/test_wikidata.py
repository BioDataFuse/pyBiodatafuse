#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the WikiData annotator."""

import unittest
from unittest.mock import Mock, patch
import datetime
import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators import wikidata
from pyBiodatafuse.annotators.wikidata import get_version_wikidata, get_gene_cellular_component
from pyBiodatafuse.constants import WIKIDATA


class TestWikidata(unittest.TestCase):
    """Test the WikiData class."""

    # TODO: Fix this when the main function changes
    def test_wikipathways_version(self):
        """Test the get_version_wikipathways."""
        now = str(datetime.datetime.now())

        obtained_version = get_version_wikidata()

        expected_version = {
            "metadata": {
                "data_version": {
                    "dataVersion": {
                        "year": now[0:4],
                        "month": now[5:7],
                    }
                },
            },
        }

        assert obtained_version == expected_version

    # @patch("pyBiodatafuse.annotators.wikidata.SPARQLWrapper.queryAndConvert")
    # def test_get_gene_cellular_component(self, mock_sparql_request):
    #     """Test the get_gene_cellular_component."""
    #     mock_sparql_request.side_effect = []

    #     wikidata.get_version_wikidata = Mock(
    #         return_value={
    #             "metadata": {"data_version": {"dataVersion": {"year": "2021", "month": "09"}}}
    #         }
    #     )  # Mock the version call
    #     wikidata.check_endpoint_wikidata = Mock(return_value=True)

    #     bridgedb_dataframe = pd.DataFrame(
    #         {
    #             "identifier": ["ALG14", "ALG2", "CHRNA1"],
    #             "identifier.source": ["HGNC", "HGNC", "HGNC"],
    #             "target": ["199857", "85365", "1134"],
    #             "target.source": ["NCBI Gene", "NCBI Gene", "NCBI Gene"],
    #         }
    #     )

    #     obtained_data, metadata = get_gene_cellular_component(bridgedb_dataframe)

    #     expected_data = pd.Series([])
    #     expected_data.name = WIKIDATA

    #     pd.testing.assert_series_equal(obtained_data[WIKIDATA], expected_data)
    #     self.assertIsInstance(metadata, dict)
