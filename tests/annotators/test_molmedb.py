#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the MolMeDB annotator."""

import unittest
from unittest.mock import Mock, patch

import pandas as pd
from numpy import nan

from pyBiodatafuse.annotators import molmedb
from pyBiodatafuse.annotators.molmedb import (
    get_compound_gene_inhibitor,
    get_gene_compound_inhibitor,
)
from pyBiodatafuse.constants import MOLMEDB_COMPOUND_PROTEIN_COL, MOLMEDB_PROTEIN_COMPOUND_COL


class TestMolMeDb(unittest.TestCase):
    """Test the MolMeDB class."""

    # TODO after MolMeDB update
    # @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
    # def test_get_version_bgee(self, mock_sparql_request):
    #     """Test the get_version_bgee function."""
    #     version_data = '{"results": {"bindings": [{"date_modified": {"type": "literal", "value": "2023-11-01"}}]}}'
    #     mock_sparql_request.return_value = json.loads(version_data)

    #     obtained_version = get_version_bgee()

    #     expected_version = {"bgee_version": "2023-11-01"}

    #     self.assertEqual(obtained_version, expected_version)

    @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
    def test_get_gene_compound_inhibitor(self, mock_sparql_request):
        """Test the get_gene_compound_inhibitor."""
        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "identifier": ["SLC17A9", "SLC17A9", "SLC25A1", "SLC25A1", "KCNJ5", "KCNJ5"],
                "identifier.source": ["HGNC", "HGNC", "HGNC", "HGNC", "HGNC", "HGNC"],
                "target": ["Q5W197", "Q9BYT1", "P53007", "B4DP62", "P48544", "A0A5J6E2W8"],
                "target.source": [
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                    "Uniprot-TrEMBL",
                ],
            }
        )

        mock_sparql_request.side_effect = [
            {
                "head": {
                    "vars": [
                        "transporterID",
                        "label",
                        "InChIKey",
                        "SMILES",
                        "chembl_id",
                        "chebi_id",
                        "drugbank_id",
                        "pubchem_compound_id",
                        "molmedb_id",
                        "source_doi",
                        "source_pmid",
                    ]
                },
                "results": {
                    "bindings": [
                        {
                            "transporterID": {
                                "type": "literal",
                                "value": "P48544",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "label": {
                                "type": "literal",
                                "value": "Euphorbiaproliferin c",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "inchikey": {
                                "type": "literal",
                                "value": "MEMULCZBXUZFOZ-UHFFFAOYSA-N",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "smiles": {
                                "type": "literal",
                                "value": "CC(=O)OC1C2(C)OCC3(C(=O)C=CC(C(C)(C)OC(C)=O)C23)C(OC(C)=O)C2C(OC(=O)C(C)C)C(C)CC21OC(C)=O",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "molmedb_id": {
                                "type": "literal",
                                "value": "MM470852",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "source_pmid": {
                                "type": "literal",
                                "value": "27441737",
                                "datatype": "http://www.w3.org/2001/XMLSchema#int",
                            },
                        },
                        {
                            "transporterID": {
                                "type": "literal",
                                "value": "P48544",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "label": {
                                "type": "literal",
                                "value": "Euphornin",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "inchikey": {
                                "type": "literal",
                                "value": "BRVXVMOWTHQKHC-LVYIKVSWSA-N",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "smiles": {
                                "type": "literal",
                                "value": "CC(=O)OC1CC(OC(C)=O)C(C)(C)/C=C\\C(C)C(OC(C)=O)C2(O)CC(C)C(OC(=O)c3ccccc3)C2/C=C\\1C",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "chebi_id": {
                                "type": "literal",
                                "value": "140105",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "molmedb_id": {
                                "type": "literal",
                                "value": "MM470853",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "source_pmid": {
                                "type": "literal",
                                "value": "30411614",
                                "datatype": "http://www.w3.org/2001/XMLSchema#int",
                            },
                        },
                        {
                            "transporterID": {
                                "type": "literal",
                                "value": "Q9BYT1",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "label": {
                                "type": "literal",
                                "value": "[dichloro(phosphono)methyl]phosphonic acid",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "inchikey": {
                                "type": "literal",
                                "value": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "smiles": {
                                "type": "literal",
                                "value": "O=P(O)(O)C(Cl)(Cl)P(=O)(O)O",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "chebi_id": {
                                "type": "literal",
                                "value": "110423",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "drugbank_id": {
                                "type": "literal",
                                "value": "DB00720",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "pubchem_compound_id": {
                                "type": "literal",
                                "value": "25419",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "molmedb_id": {
                                "type": "literal",
                                "value": "MM17483",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "source_pmid": {
                                "type": "literal",
                                "value": "28720702",
                                "datatype": "http://www.w3.org/2001/XMLSchema#int",
                            },
                        },
                    ]
                },
            }
        ]
        molmedb.check_endpoint_molmedb = Mock(return_value=True)

        obtained_data, metadata = get_gene_compound_inhibitor(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "compound_name": nan,
                        "inchikey": nan,
                        "smiles": nan,
                        "compound_cid": nan,
                        "molmedb_id": nan,
                        "source_pmid": nan,
                        "chebi_id": nan,
                        "drugbank_id": nan,
                        "uniprot_trembl_id": nan,
                    }
                ],
                [
                    {
                        "compound_name": nan,
                        "inchikey": nan,
                        "smiles": nan,
                        "compound_cid": nan,
                        "molmedb_id": nan,
                        "source_pmid": nan,
                        "chebi_id": nan,
                        "drugbank_id": nan,
                        "uniprot_trembl_id": nan,
                    }
                ],
                [
                    {
                        "compound_name": nan,
                        "inchikey": nan,
                        "smiles": nan,
                        "compound_cid": nan,
                        "molmedb_id": nan,
                        "source_pmid": nan,
                        "chebi_id": nan,
                        "drugbank_id": nan,
                        "uniprot_trembl_id": nan,
                    },
                    {
                        "compound_name": "[dichloro(phosphono)methyl]phosphonic acid",
                        "inchikey": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                        "smiles": "O=P(O)(O)C(Cl)(Cl)P(=O)(O)O",
                        "compound_cid": "25419",
                        "molmedb_id": "MM17483",
                        "source_pmid": "28720702",
                        "chebi_id": "110423",
                        "drugbank_id": "DrugBank:DB00720",
                        "uniprot_trembl_id": "Q9BYT1",
                    },
                ],
                [
                    {
                        "compound_name": nan,
                        "inchikey": nan,
                        "smiles": nan,
                        "compound_cid": nan,
                        "molmedb_id": nan,
                        "source_pmid": nan,
                        "chebi_id": nan,
                        "drugbank_id": nan,
                        "uniprot_trembl_id": nan,
                    },
                    {
                        "compound_name": "[dichloro(phosphono)methyl]phosphonic acid",
                        "inchikey": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                        "smiles": "O=P(O)(O)C(Cl)(Cl)P(=O)(O)O",
                        "compound_cid": "25419",
                        "molmedb_id": "MM17483",
                        "source_pmid": "28720702",
                        "chebi_id": "110423",
                        "drugbank_id": "DrugBank:DB00720",
                        "uniprot_trembl_id": "Q9BYT1",
                    },
                ],
                [
                    {
                        "compound_name": nan,
                        "inchikey": nan,
                        "smiles": nan,
                        "compound_cid": nan,
                        "molmedb_id": nan,
                        "source_pmid": nan,
                        "chebi_id": nan,
                        "drugbank_id": nan,
                        "uniprot_trembl_id": nan,
                    }
                ],
                [
                    {
                        "compound_name": nan,
                        "inchikey": nan,
                        "smiles": nan,
                        "compound_cid": nan,
                        "molmedb_id": nan,
                        "source_pmid": nan,
                        "chebi_id": nan,
                        "drugbank_id": nan,
                        "uniprot_trembl_id": nan,
                    }
                ],
            ]
        )
        expected_data.name = MOLMEDB_PROTEIN_COMPOUND_COL

        pd.testing.assert_series_equal(obtained_data[MOLMEDB_PROTEIN_COMPOUND_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
    def test_get_compound_gene_inhibitor(self, mock_sparql_request):
        """Test the get_compound_gene_inhibitor."""
        mock_sparql_request.side_effect = [
            {
                "head": {
                    "vars": [
                        "inhibitorInChIKey",
                        "uniprot_trembl_id",
                        "hgcn_id",
                        "source_doi",
                        "source_pmid",
                    ]
                },
                "results": {
                    "bindings": [
                        {
                            "inhibitorInChIKey": {
                                "type": "literal",
                                "value": "OVVBIIBBRZVPAL-UHFFFAOYSA-N",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "uniprot_trembl_id": {
                                "type": "literal",
                                "value": "P23975",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "hgcn_id": {
                                "type": "literal",
                                "value": "SLC6A2",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "source_pmid": {
                                "type": "literal",
                                "value": "20223878",
                                "datatype": "http://www.w3.org/2001/XMLSchema#int",
                            },
                        },
                        {
                            "inhibitorInChIKey": {
                                "type": "literal",
                                "value": "OVVBIIBBRZVPAL-UHFFFAOYSA-N",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "uniprot_trembl_id": {
                                "type": "literal",
                                "value": "P31645",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "hgcn_id": {
                                "type": "literal",
                                "value": "SLC6A4",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "source_pmid": {
                                "type": "literal",
                                "value": "20223878",
                                "datatype": "http://www.w3.org/2001/XMLSchema#int",
                            },
                        },
                        {
                            "inhibitorInChIKey": {
                                "type": "literal",
                                "value": "LEJRLSZVESQKJK-UHFFFAOYSA-N",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "uniprot_trembl_id": {
                                "type": "literal",
                                "value": "Q01959",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "hgcn_id": {
                                "type": "literal",
                                "value": "SLC6A3",
                                "datatype": "http://www.w3.org/2001/XMLSchema#string",
                            },
                            "source_pmid": {
                                "type": "literal",
                                "value": "9703474",
                                "datatype": "http://www.w3.org/2001/XMLSchema#int",
                            },
                        },
                    ]
                },
            }
        ]

        bridgedb_dataframe_compounds = pd.DataFrame(
            {
                "identifier": ["10041551", "10025195", "2153"],
                "identifier.source": ["PubChem-compound", "PubChem-compound", "PubChem-compound"],
                "target": [
                    "OVVBIIBBRZVPAL-UHFFFAOYSA-N",
                    "LEJRLSZVESQKJK-UHFFFAOYSA-N",
                    "ZFXYFBGIUFBOJW-UHFFFAOYSA-N",
                ],
                "target.source": ["InChIKey", "InChIKey", "InChIKey"],
            }
        )

        obtained_data, metadata = get_compound_gene_inhibitor(bridgedb_dataframe_compounds)

        expected_data = pd.Series(
            [
                [
                    {
                        "uniprot_trembl_id": "Q01959",
                        "hgnc_symbol": "SLC6A3",
                        "source_pmid": "9703474",
                    }
                ],
                [
                    {
                        "uniprot_trembl_id": "P23975",
                        "hgnc_symbol": "SLC6A2",
                        "source_pmid": "20223878",
                    },
                    {
                        "uniprot_trembl_id": "P31645",
                        "hgnc_symbol": "SLC6A4",
                        "source_pmid": "20223878",
                    },
                ],
                [{"uniprot_trembl_id": nan, "hgnc_symbol": nan, "source_pmid": nan}],
            ]
        )
        expected_data.name = MOLMEDB_COMPOUND_PROTEIN_COL

        pd.testing.assert_series_equal(obtained_data[MOLMEDB_COMPOUND_PROTEIN_COL], expected_data)

        self.assertIsInstance(metadata, dict)
