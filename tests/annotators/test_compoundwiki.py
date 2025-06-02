#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Wikidata annotator."""

import datetime
import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from pyBiodatafuse.annotators import compoundwiki
from pyBiodatafuse.annotators.compoundwiki import get_compound_annotations, get_version_compoundwiki
from pyBiodatafuse.constants import COMPOUNDWIKI_COL


class TestWikidata(unittest.TestCase):
    """Test the Compoundwiki class."""

    def test_version_compoundwiki(self):
        """Test the get_version_compoundwiki."""
        now = str(datetime.datetime.now())

        obtained_version = get_version_compoundwiki()

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
    def test_get_compound_annotation(self, mock_sparql_request):
        """Test the get_compound_annotation function."""
        compoundwiki.check_endpoint_compoundwiki = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["7858"],
                "identifier.source": ["PubChem-compound"],
                "target": ["7858"],
                "target.source": ["PubChem Compound"],
            }
        )

        obtained_data, metadata = get_compound_annotations(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                {
                    "AOP-Wiki Stressor ID": "9",
                    "CAS Registry Number": "107-18-6",
                    "ChEBI ID": "16605",
                    "ChEMBL ID": "CHEMBL234926",
                    "DSSTOX compound identifier": "DTXCID2044",
                    "EC number": "203-470-7",
                    "ECHA Substance Infocard ID": "100.003.156",
                    "InChI": "InChI=1S/C3H6O/c1-2-3-4/h2,4H,1,3H2",
                    "InChIKey": "XXROGKLTLUQVRX-UHFFFAOYSA-N",
                    "KEGG ID": "C02001",
                    "PubChem CID": "7858",
                    "SMILES (without stereochemistry)": "C=CCO",
                    "ToxBank Wiki": "https://wiki.toxbank.net/wiki/Allyl_alcohol",
                    "Wikidata Q identifier": "Q414553",
                    "chemical formula": "C₃H₆O",
                    "has role": "hepatotoxic agent",
                    "instance of": "chemical compound",
                    "part of": "AOP-Wiki Prototypical Stressors",
                }
            ]
        )
        expected_data.name = COMPOUNDWIKI_COL

        pd.testing.assert_series_equal(obtained_data[COMPOUNDWIKI_COL], expected_data)
        self.assertIsInstance(metadata, dict)
