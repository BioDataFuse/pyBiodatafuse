"""
Unit tests for the BDFGraph class in the pyBiodatafuse package.

This module contains unit tests for the BDFGraph class, which is responsible for generating RDF graphs from dataframes.

Classes:
    TestBDFGraph: Unit test case for the BDFGraph class.

Constants:
    DATA: Mock data loaded from a JSON file for testing.
    METADATA: Mock metadata loaded from a JSON file for testing.
"""

import json
import os
import unittest

import pandas as pd

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.rdf import BDFGraph

DATA = pd.read_json(os.path.join("tests", "graph", "data", "combined_df_mock_data.json"))
with open(
    os.path.join("tests", "graph", "data", "combined_metadata_mock_data.json"),
    "r",
    encoding="utf-8",
) as file:
    METADATA = json.load(file)


class TestBDFGraph(unittest.TestCase):
    """
    Set unit test case for the BDFGraph class.

    Methods:
        setUp: Sets up the test environment.
        test_generate_rdf: Tests the generate_rdf method of the BDFGraph class.
    """

    def setUp(self):
        """
        Set up the test environment.

        This method initializes the BDFGraph instance with a base URI, version IRI, author name, and ORCID.
        """
        self.base_uri = "http://example.org/"
        self.version_iri = "http://example.org/version"
        self.author = "Author Name"
        self.orcid = "0000-0000-0000-0000"
        self.graph = BDFGraph(self.base_uri, self.version_iri, self.author, self.orcid)

    def test_generate_rdf(self):
        """
        Test the generate_rdf method of the BDFGraph class.

        This method verifies that the RDF graph is generated correctly from the provided dataframe and metadata.
        """
        df = pd.DataFrame(DATA)
        self.graph.generate_rdf(df, metadata=METADATA)
        self.assertTrue(len(self.graph) > 0)


if __name__ == "__main__":
    unittest.main()
