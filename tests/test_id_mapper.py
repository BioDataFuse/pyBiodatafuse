import unittest
import pandas as pd
from pyBiodatafuse.id_mapper import match_input_datasource


class TestMatchInputDatasource(unittest.TestCase):
    def test_single_match(self):
        identifiers = pd.Series(["ENSG00000139618"])
        result = match_input_datasource(identifiers)
        self.assertEqual(result, "Ensembl", "Expected 'Ensembl' for ENSG identifier")

    def test_multiple_matches(self):
        identifiers = pd.Series(["12345"])
        with self.assertRaises(ValueError) as context:
            match_input_datasource(identifiers)
        self.assertIn("Multiple data sources match the provided identifiers", str(context.exception))

    def test_empty_series(self):
        identifiers = pd.Series([])
        with self.assertRaises(ValueError) as context:
            match_input_datasource(identifiers)
        self.assertIn("The identifiers series is empty.", str(context.exception))

# Missing bridgedb tests


if __name__ == "__main__":
    unittest.main()
