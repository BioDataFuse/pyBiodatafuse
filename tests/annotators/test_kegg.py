#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the KEGG annotator."""

import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import kegg
from pyBiodatafuse.constants import KEGG_PATHWAY_COL


class TestKEGG(unittest.TestCase):
    """Test the KEGG annotation module."""

    def test_get_pathways(self):
        """Test the get_pathways function."""
        kegg.check_version_kegg = Mock(return_value={"source_version": "KEGG RDF 20240315"})
        kegg.check_endpoint_kegg = Mock(return_value=True)

        kegg.batch_request = Mock(
            return_value=open("tests/annotators/data/kegg_mock_data.txt").read()
        )

        input_df = pd.DataFrame(
            {
                "identifier": ["ENSMUSG00000067274"],
                "identifier.source": ["Ensembl"],
                "target": ["11837"],
                "target.source": ["NCBI Gene"],
            }
        )

        result_df, metadata = kegg.get_pathways(input_df)

        expected_pathways = pd.Series(
            [
                [
                    {
                        "pathway_id": "path:mmu03010",
                        "pathway_label": "Ribosome - Mus musculus (house mouse)",
                        "pathway_gene_counts": 278,
                        "pathway_compounds": [{"KEGG_id": None}],
                    },
                    {
                        "pathway_id": "path:mmu05171",
                        "pathway_label": "Coronavirus disease - COVID-19 - Mus musculus (house mouse)",
                        "pathway_gene_counts": 333,
                        "pathway_compounds": [
                            {"KEGG_id": "C00027"},
                            {"KEGG_id": "C00046"},
                            {"KEGG_id": "C00165"},
                            {"KEGG_id": "C00290"},
                            {"KEGG_id": "C00533"},
                            {"KEGG_id": "C00704"},
                            {"KEGG_id": "C00873"},
                            {"KEGG_id": "C02135"},
                            {"KEGG_id": "C15850"},
                            {"KEGG_id": "C20640"},
                        ],
                    },
                ]
            ]
        )

        expected_pathways.name = KEGG_PATHWAY_COL

        pd.testing.assert_series_equal(result_df[KEGG_PATHWAY_COL], expected_pathways)
        self.assertIsInstance(metadata, dict)
