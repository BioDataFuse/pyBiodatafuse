#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the CompoundWiki annotator."""

import unittest
from unittest.mock import Mock

import pandas as pd

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.annotators import compoundwiki


class TestCompoundwiki(unittest.TestCase):
    """Test the CompoundWiki class."""

    def test_get_pubchem_compoundwiki_annotations(self):
        """Test the get_compound_annotations function when pubchem is the input."""
        compoundwiki.check_endpoint_compoundwiki = Mock(return_value=True)

        pubchem_dataframe = pd.DataFrame(
            {
                "identifier": ["KCNJ5"],
                "identifier.source": ["HGNC"],
                "target": ["P38398"],
                "target.source": ["Uniprot-TrEMBL"],
                "PubChem_assays": [
                    [
                        {
                            "pubchem_assay_id": "AID:504669",
                            "assay_type": "IC50",
                            "outcome": "active",
                            "compound_cid": "CID:2157",
                            "compound_name": "2-nitro-N-quinolin-8-ylbenzenesulfonamide",
                            "smiles": "C1=CC=C(C(=C1)[N+](=O)[O-])S(=O)(=O)NC2=CC=CC3=C2N=CC=C3",
                            "inchi": "InChI=1S/C15H11N3O4S/c19-18(20)13-8-1-2-9-14(13)23(21,22)17-12-7-3-5-11-6-4-10-16-15(11)12/h1-10,17H",
                        }
                    ]
                ],
            }
        )

        obtained_data, metadata = compoundwiki.get_compound_annotations(pubchem_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "pubchem_assay_id": "AID:504669",
                        "assay_type": "IC50",
                        "outcome": "active",
                        "compound_cid": "CID:2157",
                        "compound_name": "2-nitro-N-quinolin-8-ylbenzenesulfonamide",
                        "smiles": "C1=CC=C(C(=C1)[N+](=O)[O-])S(=O)(=O)NC2=CC=CC3=C2N=CC=C3",
                        "inchi": "InChI=1S/C15H11N3O4S/c19-18(20)13-8-1-2-9-14(13)23(21,22)17-12-7-3-5-11-6-4-10-16-15(11)12/h1-10,17H",
                        "CompoundWiki_compounds": [
                            {
                                "AOP-Wiki Stressor ID": "95",
                                "ArrayExpress identifier": "E-MTAB-665",
                                "CAS Registry Number": "1951-25-3",
                                "ChEBI ID": "2663",
                                "ChEMBL ID": "CHEMBL633",
                                "DSSTOX compound identifier": "DTXCID702592",
                                "EC number": "217-772-1",
                                "ECHA Substance Infocard ID": "100.016.157",
                                "InChI": "InChI=1S/C25H29I2NO3/c1-4-7-11-22-23(18-10-8-9-12-21(18)31-22)24(29)17-15-19(26)25(20(27)16-17)30-14-13-28(5-2)6-3/h8-10,12,15-16H,4-7,11,13-14H2,1-3H3",
                                "InChIKey": "IYIKLHRQXLHMJQ-UHFFFAOYSA-N",
                                "KEGG ID": "D02910",
                                "PubChem CID": "2157",
                                "SMILES (without stereochemistry)": "CCCCC1=C(C2=CC=CC=C2O1)C(=O)C3=CC(=C(C(=C3)I)OCCN(CC)CC)I",
                                "ToxBank Wiki": "https://wiki.toxbank.net/wiki/Amiodarone",
                                "Wikidata Q identifier": "Q410061",
                                "chemical formula": "C₂₅H₂₉I₂NO₃",
                                "has role": "hepatotoxic agent",
                                "has salt": "amiodarone hydrochloride",
                                "instance of": "chemical compound",
                                "mass": "645.3125",
                                "part of": "EFSA TXG-MAP collection",
                                "compound label": "amiodarone",
                                "input_identifier": "2157",
                            }
                        ],
                    }
                ]
            ]
        )
        expected_data.name = Cons.PUBCHEM_COMPOUND_ASSAYS_COL

        pd.testing.assert_series_equal(
            obtained_data[Cons.PUBCHEM_COMPOUND_ASSAYS_COL], expected_data
        )
        self.assertIsInstance(metadata, dict)
