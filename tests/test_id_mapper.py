"""Unit tests for the id_mapper module."""

import unittest
import warnings

import pandas as pd

from pyBiodatafuse.id_mapper import match_input_datasource

# Suppress specific warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pydantic")


class TestMatchInputDatasource(unittest.TestCase):
    """Test cases for the match_input_datasource function."""

    def test_single_match(self):
        """Test that a single identifier matches the correct data source."""
        identifiers = pd.Series(["ENSG00000139618"])
        result = match_input_datasource(identifiers)
        self.assertEqual(result, "Ensembl", "Expected 'Ensembl' for ENSG identifier")

    def test_multiple_matches(self):
        """Test that multiple matches raise a ValueError."""
        identifiers = pd.Series(["12345"])
        with self.assertRaises(ValueError) as context:
            match_input_datasource(identifiers)
        self.assertIn(
            "Multiple data sources match the provided identifiers", str(context.exception)
        )

    def test_empty_series(self):
        """Test that an empty series raises a ValueError."""
        identifiers = pd.Series([])
        with self.assertRaises(ValueError) as context:
            match_input_datasource(identifiers)
        self.assertIn("The identifiers series is empty.", str(context.exception))


if __name__ == "__main__":
    unittest.main()
