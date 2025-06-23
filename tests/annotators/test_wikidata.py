#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Wikidata annotator."""

import datetime
import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import wikidata
from pyBiodatafuse.annotators.wikidata import get_gene_cellular_component, get_version_wikidata
from pyBiodatafuse.constants import WIKIDATA_CC_COL

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestWikidata(unittest.TestCase):
    """Test the Wikidata class."""

    # TODO: Fix this when the main function changes
    def test_wikidata_version(self):
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

    @patch("pyBiodatafuse.annotators.wikidata.SPARQLWrapper.queryAndConvert")
    def test_get_gene_cellular_component(self, mock_sparql_request):
        """Test the get_gene_cellular_component."""
        with open(os.path.join(data_file_folder, "wikidata_mock_data.json")) as f:
            mock_data = json.load(f)

        mock_sparql_request.side_effect = [mock_data]
        wikidata.check_endpoint_wikidata = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["AGRN", "ATXN7"],
                "identifier.source": ["HGNC", "HGNC"],
                "target": ["375790", "6314"],
                "target.source": ["NCBI Gene", "NCBI Gene"],
            }
        )

        obtained_data, metadata = get_gene_cellular_component(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "wikidata_id": "Q29548",
                        "wikidata_label": "plasma membrane",
                        "go_id": "GO:0005886",
                    },
                    {
                        "wikidata_id": "Q32846",
                        "wikidata_label": "basement membrane",
                        "go_id": "GO:0005604",
                    },
                    {
                        "wikidata_id": "Q189073",
                        "wikidata_label": "cell junction",
                        "go_id": "GO:0030054",
                    },
                    {"wikidata_id": "Q187181", "wikidata_label": "synapse", "go_id": "GO:0045202"},
                    {
                        "wikidata_id": "Q193825",
                        "wikidata_label": "extracellular matrix",
                        "go_id": "GO:0031012",
                    },
                    {"wikidata_id": "Q220599", "wikidata_label": "cytosol", "go_id": "GO:0005829"},
                    {
                        "wikidata_id": "Q903634",
                        "wikidata_label": "extracellular exosome",
                        "go_id": "GO:0070062",
                    },
                    {
                        "wikidata_id": "Q54810453",
                        "wikidata_label": "collagen-containing extracellular matrix",
                        "go_id": "GO:0062023",
                    },
                    {
                        "wikidata_id": "Q14758889",
                        "wikidata_label": "Golgi lumen",
                        "go_id": "GO:0005796",
                    },
                    {
                        "wikidata_id": "Q14866090",
                        "wikidata_label": "lysosomal lumen",
                        "go_id": "GO:0043202",
                    },
                    {
                        "wikidata_id": "Q14327652",
                        "wikidata_label": "integral component of membrane",
                        "go_id": "GO:0016021",
                    },
                    {
                        "wikidata_id": "Q14645596",
                        "wikidata_label": "extracellular region",
                        "go_id": "GO:0005576",
                    },
                    {
                        "wikidata_id": "Q14349455",
                        "wikidata_label": "membrane",
                        "go_id": "GO:0016020",
                    },
                ],
                [
                    {"wikidata_id": "Q30869", "wikidata_label": "nucleolus", "go_id": "GO:0005730"},
                    {"wikidata_id": "Q40260", "wikidata_label": "nucleus", "go_id": "GO:0005634"},
                    {"wikidata_id": "Q220599", "wikidata_label": "cytosol", "go_id": "GO:0005829"},
                    {
                        "wikidata_id": "Q154626",
                        "wikidata_label": "cytoskeleton",
                        "go_id": "GO:0005856",
                    },
                    {"wikidata_id": "Q79899", "wikidata_label": "cytoplasm", "go_id": "GO:0005737"},
                    {
                        "wikidata_id": "Q14817956",
                        "wikidata_label": "nucleoplasm",
                        "go_id": "GO:0005654",
                    },
                    {
                        "wikidata_id": "Q14888122",
                        "wikidata_label": "microtubule cytoskeleton",
                        "go_id": "GO:0015630",
                    },
                    {
                        "wikidata_id": "Q906755",
                        "wikidata_label": "nuclear matrix",
                        "go_id": "GO:0016363",
                    },
                ],
            ]
        )
        expected_data.name = WIKIDATA_CC_COL

        pd.testing.assert_series_equal(obtained_data[WIKIDATA_CC_COL], expected_data)
        self.assertIsInstance(metadata, dict)
