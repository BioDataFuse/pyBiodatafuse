#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the WikiPathways annotator."""

import pandas as pd
import pytest
from numpy import nan

from pyBiodatafuse.annotators.wikipathways import get_gene_wikipathway, get_version_wikipathways


def test_wikipathways_version():
    """Test the get_version_wikipathways."""
    obtained_version = get_version_wikipathways()

    expected_version = {"wikipathways_version": "WikiPathways RDF 20231210"}

    assert obtained_version == expected_version


def test_get_gene_wikipathway(bridgedb_dataframe):
    """Test the get_gene_wikipathway."""
    obtained_data, metadata = get_gene_wikipathway(bridgedb_dataframe)

    expected_data = pd.Series(
        [
            [
                {
                    "pathwayId": "WP5153",
                    "pathwayLabel": "N-glycan biosynthesis",
                    "pathwayGeneCount": "57",
                }
            ],
            [
                {
                    "pathwayId": "WP5153",
                    "pathwayLabel": "N-glycan biosynthesis",
                    "pathwayGeneCount": "57",
                }
            ],
            [{"pathwayId": nan, "pathwayLabel": nan, "pathwayGeneCount": nan}],
        ]
    )
    expected_data.name = "WikiPathways"

    pd.testing.assert_series_equal(obtained_data["WikiPathways"], expected_data)


@pytest.fixture(scope="module")
def bridgedb_dataframe():
    """Reusable sample Pandas DataFrame to be used as input for the tests."""
    return pd.DataFrame(
        {
            "identifier": ["ALG14", "ALG2", "CHRNA1"],
            "identifier.source": ["HGNC", "HGNC", "HGNC"],
            "target": ["199857", "85365", "1134"],
            "target.source": ["NCBI Gene", "NCBI Gene", "NCBI Gene"],
        }
    )
