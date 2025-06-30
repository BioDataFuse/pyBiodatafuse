#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the AOP annotator."""

import json
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators.aopwiki import get_aops
from pyBiodatafuse.constants import AOPWIKI_COL, AOPWIKIRDF


class TestAOPAnnotator(unittest.TestCase):
    """Test the AOP annotator."""

    @patch("pyBiodatafuse.annotators.aopwiki.SPARQLWrapper.queryAndConvert")
    def test_get_aops(self, mock_sparql_request):
        """Test the get_aops function."""
        with open("tests/annotators/data/aop_mock_data.json", "r", encoding="utf-8") as file:
            mock_sparql_request.return_value = json.load(file)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["7350"],
                "identifier.source": ["Entrez Gene"],
                "target": ["ENSG00000109424"],
                "target.source": ["Ensembl"],
            }
        )

        obtained_data, metadata = get_aops(
            bridgedb_df=bridgedb_dataframe,
            input_type="gene",
            input_identifier="Ensembl",
        )

        # Validate the content of the DataFrame
        with open("tests/annotators/data/aop_mock_res.json", "r", encoding="utf-8") as file:
            expected_data = pd.Series(json.load(file))
        expected_data.name = AOPWIKI_COL

        pd.testing.assert_series_equal(obtained_data[AOPWIKI_COL], expected_data)
        self.assertIsInstance(metadata, dict)
