#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the KEGG annotator."""

import unittest
from unittest.mock import Mock, patch

import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators import kegg
from pyBiodatafuse.annotators.kegg import check_version_kegg, get_pathways
from pyBiodatafuse.constants import KEGG_COL


class TestKEGG(unittest.TestCase):
    """Test the KEGG annotation module."""

    def test_get_pathways(self):
        """Test the get_pathways function."""
        kegg.check_version_kegg = Mock(return_value={"source_version": "KEGG RDF 20240315"})
        kegg.check_endpoint_kegg = Mock(return_value=True)

        input_df = pd.DataFrame(
            {
                "identifier": ["ENSMUSG00000067274"],
                "identifier.source": ["Ensembl"],
                "target": ["11837"],
                "target.source": ["NCBI Gene"],
            }
        )

        result_df, metadata = get_pathways(input_df)

        expected_pathways = [
            {
                "pathway_id": "path:mmu03010",
                "pathway_label": "Ribosome - Mus musculus (house mouse)",
                "gene_count": 274,
                "compounds": [{"KEGG_identifier": None}],
            },
            {
                "pathway_id": "path:mmu05171",
                "pathway_label": "Coronavirus disease - COVID-19 - Mus musculus (house mouse)",
                "gene_count": 333,
                "compounds": [
                    {"KEGG_identifier": "C00027"},
                    {"KEGG_identifier": "C00046"},
                    {"KEGG_identifier": "C00165"},
                    {"KEGG_identifier": "C00290"},
                    {"KEGG_identifier": "C00533"},
                    {"KEGG_identifier": "C00704"},
                    {"KEGG_identifier": "C00873"},
                    {"KEGG_identifier": "C02135"},
                    {"KEGG_identifier": "C15850"},
                    {"KEGG_identifier": "C20640"},
                ],
            },
            {"gene_count": 0, "compounds": [{"KEGG_identifier": None}]},
        ]

        expected_series = pd.Series([[*expected_pathways]])
        expected_series.name = KEGG_COL

        pd.testing.assert_series_equal(result_df[KEGG_COL], expected_series)
        self.assertIsInstance(metadata, dict)


if __name__ == "__main__":
    unittest.main()
