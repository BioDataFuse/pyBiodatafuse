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

    @patch("pyBiodatafuse.annotators.compoundwiki.SPARQLWrapper.queryAndConvert")
    def test_get_compound_annotation(self, mock_sparql_request):
        """Test the get_compound_annotation function."""
        compoundwiki.check_endpoint_compoundwiki = Mock(return_value=True)

        mock_sparql_request.return_value = {
            "head": {
                "vars": ["compound", "compoundLabel", "propEntity", "propLabel", "val", "valLabel"]
            },
            "results": {
                "bindings": [
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P36",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "AOP-Wiki Stressor ID",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "9"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P23",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "CAS Registry Number",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "107-18-6"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P28",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "ChEBI ID"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "16605"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P41",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "ChEMBL ID"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "CHEMBL234926"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P22",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "DSSTOX compound identifier",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "DTXCID2044"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P43",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "EC number"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "203-470-7"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P44",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "ECHA Substance Infocard ID",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "100.003.156"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P9",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "InChI"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "InChI=1S/C3H6O/c1-2-3-4/h2,4H,1,3H2"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P10",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "InChIKey"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "XXROGKLTLUQVRX-UHFFFAOYSA-N"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P27",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "KEGG ID"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "C02001"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P13",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "PubChem CID"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "7858"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P12",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "SMILES (without stereochemistry)",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "C=CCO"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P4",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "ToxBank Wiki"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {
                            "type": "uri",
                            "value": "https://wiki.toxbank.net/wiki/Allyl_alcohol",
                        },
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P5",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "Wikidata Q identifier",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "Q414553"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P3",
                        },
                        "propLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "chemical formula",
                        },
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {"type": "literal", "value": "C₃H₆O"},
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P17",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "has role"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q42",
                        },
                        "valLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "hepatotoxic agent",
                        },
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P1",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "instance of"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q2",
                        },
                        "valLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "chemical compound",
                        },
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P21",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "part of"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q53",
                        },
                        "valLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "SEURAT-1 Gold Compounds",
                        },
                    },
                    {
                        "propEntity": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/P21",
                        },
                        "propLabel": {"xml:lang": "en", "type": "literal", "value": "part of"},
                        "compound": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q9",
                        },
                        "compoundLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "allyl alcohol",
                        },
                        "val": {
                            "type": "uri",
                            "value": "https://compoundcloud.wikibase.cloud/entity/Q4761",
                        },
                        "valLabel": {
                            "xml:lang": "en",
                            "type": "literal",
                            "value": "AOP-Wiki Prototypical Stressors",
                        },
                    },
                ]
            },
        }

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
                        "compound label": "allyl alcohol",
                    }
                ]
            ]
        )
        expected_data.name = COMPOUNDWIKI_COL

        pd.testing.assert_series_equal(obtained_data[COMPOUNDWIKI_COL], expected_data)
        self.assertIsInstance(metadata, dict)
