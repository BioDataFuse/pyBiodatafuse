#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the WikiPathways annotator."""
import unittest
from unittest.mock import patch

import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators.pubchem import get_protein_molecule_screened
from pyBiodatafuse.constants import PUBCHEM


class TestPubchem(unittest.TestCase):
    """Test the PubChem class."""

    # TODO after Pubchem update
    # @patch("pyBiodatafuse.annotators.pubchem.SPARQLWrapper.queryAndConvert")
    # def test_get_version_bgee(self, mock_sparql_request):
    #     """Test the get_version_bgee function."""
    #     version_data = '{"results": {"bindings": [{"date_modified": {"type": "literal", "value": "2023-11-01"}}]}}'
    #     mock_sparql_request.return_value = json.loads(version_data)

    #     obtained_version = get_version_bgee()

    #     expected_version = {"bgee_version": "2023-11-01"}

    #     self.assertEqual(obtained_version, expected_version)

    @patch("pyBiodatafuse.annotators.pubchem.SPARQLWrapper.queryAndConvert")
    def test_get_protein_molecule_screened(self, mock_sparql_request):
        """Test the get_protein_molecule_screened."""
        mock_sparql_request.side_effect = [
            {
                "head": {
                    "vars": [
                        "upProt",
                        "assay_type",
                        "outcome",
                        "compound_cid",
                        "ref_cit",
                        "compound_name",
                        "SMILES",
                        "InChI",
                    ]
                },
                "results": {
                    "bindings": [
                        {
                            "upProt": {
                                "type": "uri",
                                "value": "http://purl.uniprot.org/uniprot/P46098",
                            },
                            "assay_type": {
                                "type": "uri",
                                "value": "http://www.bioassayontology.org/bao#BAO_0000192",
                            },
                            "outcome": {
                                "type": "uri",
                                "value": "http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#active",
                            },
                            "compound_cid": {
                                "type": "uri",
                                "value": "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID9911844",
                            },
                            "ref_cit": {
                                "type": "literal",
                                "value": "Kikuchi C, Suzuki H, Hiranuma T, Koyama M.",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "compound_name": {
                                "type": "literal",
                                "value": "DR-4485 free base",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "SMILES": {
                                "type": "literal",
                                "value": "C1CC2=C(C=CC3=C2C(C1)(C(=O)N3)CCCCN4CCC(=CC4)C5=CC=C(C=C5)Cl)Cl",
                                "xml:lang": "en",
                            },
                            "InChI": {
                                "type": "literal",
                                "value": "InChI=1S/C26H28Cl2N2O/c27",
                                "xml:lang": "en",
                            },
                            "target_count": {
                                "type": "literal",
                                "value": "1",
                                "datatype": "http://www.w3.org/2001/XMLSchema#integer",
                            },
                        },
                        {
                            "upProt": {
                                "type": "uri",
                                "value": "http://purl.uniprot.org/uniprot/P00533",
                            },
                            "assay_type": {
                                "type": "uri",
                                "value": "http://www.bioassayontology.org/bao#BAO_0000186",
                            },
                            "outcome": {
                                "type": "uri",
                                "value": "http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#inconclusive",
                            },
                            "compound_cid": {
                                "type": "uri",
                                "value": "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID5151543",
                            },
                            "compound_name": {
                                "type": "literal",
                                "value": "1-(4-Methylpiperazin-1-yl)anthra-9,10-quinone",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "SMILES": {
                                "type": "literal",
                                "value": "CN1CCN(CC1)C2=CC=CC3=C2C(=O)C4=CC=CC=C4C3=O",
                                "xml:lang": "en",
                            },
                            "InChI": {
                                "type": "literal",
                                "value": "InChI=1S/C19H18N2O2/c1",
                                "xml:lang": "en",
                            },
                            "target_count": {
                                "type": "literal",
                                "value": "1",
                                "datatype": "http://www.w3.org/2001/XMLSchema#integer",
                            },
                        },
                        {
                            "upProt": {
                                "type": "uri",
                                "value": "http://purl.uniprot.org/uniprot/P00533",
                            },
                            "assay_type": {
                                "type": "uri",
                                "value": "http://www.bioassayontology.org/bao#BAO_0000186",
                            },
                            "outcome": {
                                "type": "uri",
                                "value": "http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#inconclusive",
                            },
                            "compound_cid": {
                                "type": "uri",
                                "value": "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID6706",
                            },
                            "compound_name": {
                                "type": "literal",
                                "value": "1-(Methylamino)anthraquinone",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "SMILES": {
                                "type": "literal",
                                "value": "CNC1=CC=CC2=C1C(=O)C3=CC=CC=C3C2=O",
                                "xml:lang": "en",
                            },
                            "InChI": {
                                "type": "literal",
                                "value": "InChI=1S/C15H11NO2/c1",
                                "xml:lang": "en",
                            },
                            "target_count": {
                                "type": "literal",
                                "value": "1",
                                "datatype": "http://www.w3.org/2001/XMLSchema#integer",
                            },
                        },
                    ]
                },
            }
        ]

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["EGFR", "EGFR", "EGFR", "HTR3A", "HTR3A"],
                "identifier.source": ["HGNC", "HGNC", "HGNC", "HGNC", "HGNC"],
                "target": ["C9JYS6", "P00533", "Q504U8", "A0A0B4J205", "P46098"],
                "target.source": [
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                ],
            }
        )

        obtained_data, metadata = get_protein_molecule_screened(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "assay_type": "AC50",
                        "outcome": "inconclusive",
                        "compound_cid": 5151543,
                        "compound_name": "1-(4-Methylpiperazin-1-yl)anthra-9,10-quinone",
                        "SMILES": "CN1CCN(CC1)C2=CC=CC3=C2C(=O)C4=CC=CC=C4C3=O",
                        "InChI": "InChI=1S/C19H18N2O2/c1",
                        "ref_cit": nan,
                    },
                    {
                        "assay_type": "AC50",
                        "outcome": "inconclusive",
                        "compound_cid": 6706,
                        "compound_name": "1-(Methylamino)anthraquinone",
                        "SMILES": "CNC1=CC=CC2=C1C(=O)C3=CC=CC=C3C2=O",
                        "InChI": "InChI=1S/C15H11NO2/c1",
                        "ref_cit": nan,
                    },
                ],
                [
                    {
                        "assay_type": "Ki",
                        "outcome": "active",
                        "compound_cid": 9911844,
                        "compound_name": "DR-4485 free base",
                        "SMILES": "C1CC2=C(C=CC3=C2C(C1)(C(=O)N3)CCCCN4CCC(=CC4)C5=CC=C(C=C5)Cl)Cl",
                        "InChI": "InChI=1S/C26H28Cl2N2O/c27",
                        "ref_cit": "Kikuchi C, Suzuki H, Hiranuma T, Koyama M.",
                    }
                ],
            ]
        )
        expected_data.name = f"{PUBCHEM}_assays"

        pd.testing.assert_series_equal(obtained_data[f"{PUBCHEM}_assays"], expected_data)
        self.assertIsInstance(metadata, dict)
