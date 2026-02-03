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

    def test_get_aops(self):
        """Test the get_aops function."""
        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["4193"],
                "identifier.source": ["Entrez Gene"],
                "target": ["ENSG00000135679"],
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
        with open("tests/annotators/data/aop_mock_res_simple.json", "r", encoding="utf-8") as file:
            expected_data = pd.DataFrame(json.load(file))
        expected_data.name = Cons.AOPWIKI_GENE_COL
        expected_data.index = obtained_data.index
        self.assertEqual(type(obtained_data), pd.DataFrame)
        self.assertEqual(type(expected_data), pd.DataFrame)
        pd.testing.assert_frame_equal(obtained_data, expected_data)


if __name__ == "__main__":
    unittest.main()
