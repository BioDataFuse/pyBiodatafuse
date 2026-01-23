#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the AOP annotator."""

import json
import unittest
from unittest.mock import Mock, patch

import pandas as pd

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.annotators.aopwiki import get_aops


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
        )

        # Validate the structure of the returned DataFrame
        self.assertIn(Cons.AOPWIKI_GENE_COL, obtained_data.columns)
        self.assertEqual(len(obtained_data), len(bridgedb_dataframe))

        # Validate the metadata
        self.assertIsInstance(metadata, dict)
        self.assertIn("datasource", metadata)
        self.assertEqual(metadata["datasource"], Cons.AOPWIKIRDF)
        self.assertIn("query", metadata)
        self.assertIn("size", metadata["query"])
        self.assertEqual(metadata["query"]["size"], 1)

        # Validate the content of the DataFrame
        with open("tests/annotators/data/aop_mock_res.json", "r", encoding="utf-8") as file:
            expected_data = pd.Series(json.load(file)[Cons.AOPWIKI_GENE_COL])
        expected_data.name = Cons.AOPWIKI_GENE_COL
        expected_data.index = obtained_data.index

        pd.testing.assert_series_equal(obtained_data[Cons.AOPWIKI_GENE_COL], expected_data)


if __name__ == "__main__":
    unittest.main()
