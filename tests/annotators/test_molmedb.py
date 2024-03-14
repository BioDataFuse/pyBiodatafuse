#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the MolMeDB annotator."""

from unittest.mock import patch

import pandas as pd
import pytest
from numpy import nan

from pyBiodatafuse.annotators.molmedb import (
    get_gene_compound_inhibitor,
    get_compound_gene_inhibitor,
)

# MolMeDB still not versioned
# TODO after MolMeDB update
# @patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
# def test_molmedb_version(mock_sparql_request):
#     """Test the get_version_molmedb."""
#     mock_sparql_request.return_value = {
#         "head": {"link": [], "vars": ["title"]},
#         "results": {
#             "distinct": False,
#             "ordered": True,
#             "bindings": [{"title": {"type": "literal", "value": "XXX"}}],
#         },
#     }

#     obtained_version = get_version_molmedb()

#     expected_version = {"molmedb_version": "XXX"}

#     assert obtained_version == expected_version


@patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
def test_get_gene_compound_inhibitor(mock_sparql_request, bridgedb_dataframe_genes):
    """Test the get_gene_compound_inhibitor."""
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
                    "pdb_ligand_id",
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
                        "InChIKey": {
                            "type": "literal",
                            "value": "MEMULCZBXUZFOZ-UHFFFAOYSA-N",
                            "datatype": "http://www.w3.org/2001/XMLSchema#string",
                        },
                        "SMILES": {
                            "type": "literal",
                            "value": "CC(=O)OC1C2(C)OCC3(C(=O)C=CC(C(C)(C)OC",
                            "datatype": "http://www.w3.org/2001/XMLSchema#string",
                        },
                        "molmedb_id": {
                            "type": "literal",
                            "value": "MM470852",
                            "datatype": "http://www.w3.org/2001/XMLSchema#string",
                        },
                        "source_doi": {
                            "type": "literal",
                            "value": "10.1021/acs.jnatprod.6b00260",
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
                        "InChIKey": {
                            "type": "literal",
                            "value": "BRVXVMOWTHQKHC-LVYIKVSWSA-N",
                            "datatype": "http://www.w3.org/2001/XMLSchema#string",
                        },
                        "SMILES": {
                            "type": "literal",
                            "value": "CC(=O)OC1CC(OC(C)=O)C(C)(C)/C=C",
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
                        "source_doi": {
                            "type": "literal",
                            "value": "10.1021/acs.jnatprod.8b00500",
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
                        "InChIKey": {
                            "type": "literal",
                            "value": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                            "datatype": "http://www.w3.org/2001/XMLSchema#string",
                        },
                        "SMILES": {
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
                        "source_doi": {
                            "type": "literal",
                            "value": "10.1073/pnas.1704847114",
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

    obtained_data, metadata = get_gene_compound_inhibitor(bridgedb_dataframe_genes)

    expected_data = pd.Series(
        [
            [
                {
                    "compound_name": "Euphorbiaproliferin c",
                    "InChIKey": "MEMULCZBXUZFOZ-UHFFFAOYSA-N",
                    "SMILES": "CC(=O)OC1C2(C)OCC3(C(=O)C=CC(C(C)(C)OC",
                    "compound_cid": nan,
                    "molmedb_id": "MM470852",
                    "source_doi": "doi:10.1021/acs.jnatprod.6b00260",
                    "source_pmid": 27441737,
                    "chebi_id": nan,
                    "drugbank_id": nan,
                },
                {
                    "compound_name": "Euphornin",
                    "InChIKey": "BRVXVMOWTHQKHC-LVYIKVSWSA-N",
                    "SMILES": "CC(=O)OC1CC(OC(C)=O)C(C)(C)/C=C",
                    "compound_cid": nan,
                    "molmedb_id": "MM470853",
                    "source_doi": "doi:10.1021/acs.jnatprod.8b00500",
                    "source_pmid": 30411614,
                    "chebi_id": "140105",
                    "drugbank_id": nan,
                },
            ],
            [
                {
                    "compound_name": "[dichloro(phosphono)methyl]phosphonic acid",
                    "InChIKey": "ACSIXWWBWUQEHA-UHFFFAOYSA-N",
                    "SMILES": "O=P(O)(O)C(Cl)(Cl)P(=O)(O)O",
                    "compound_cid": 25419,
                    "molmedb_id": "MM17483",
                    "source_doi": "doi:10.1073/pnas.1704847114",
                    "source_pmid": 28720702,
                    "chebi_id": "110423",
                    "drugbank_id": "DB00720",
                }
            ],
            [
                {
                    "compound_name": nan,
                    "InChIKey": nan,
                    "SMILES": nan,
                    "compound_cid": nan,
                    "molmedb_id": nan,
                    "source_doi": nan,
                    "source_pmid": nan,
                    "chebi_id": nan,
                    "drugbank_id": nan,
                }
            ],
        ]
    )
    expected_data.name = "transporter_inhibitor"

    pd.testing.assert_series_equal(obtained_data["transporter_inhibitor"], expected_data)


@patch("pyBiodatafuse.annotators.molmedb.SPARQLWrapper.queryAndConvert")
def test_get_compound_gene_inhibitor(mock_sparql_request, bridgedb_dataframe_compounds):
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
                        "source_doi": {
                            "type": "literal",
                            "value": "10.1124/jpet.110.166264",
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
                        "source_doi": {
                            "type": "literal",
                            "value": "10.1124/jpet.110.166264",
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
                        "source_doi": {
                            "type": "literal",
                            "value": "10.1021/jm980066t",
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

    obtained_data, metadata = get_compound_gene_inhibitor(bridgedb_dataframe_compounds)

    expected_data = pd.Series(
        [
            [
                {
                    "uniprot_trembl_id": "Q01959",
                    "hgcn_id": "SLC6A3",
                    "source_doi": "doi:10.1021/jm980066t",
                    "source_pmid": "9703474",
                }
            ],
            [
                {
                    "uniprot_trembl_id": "P23975",
                    "hgcn_id": "SLC6A2",
                    "source_doi": "doi:10.1124/jpet.110.166264",
                    "source_pmid": "20223878",
                },
                {
                    "uniprot_trembl_id": "P31645",
                    "hgcn_id": "SLC6A4",
                    "source_doi": "doi:10.1124/jpet.110.166264",
                    "source_pmid": "20223878",
                },
            ],
            [{"uniprot_trembl_id": nan, "hgcn_id": nan, "source_doi": nan, "source_pmid": nan}],
        ]
    )
    expected_data.name = "transporter_inhibited"

    pd.testing.assert_series_equal(obtained_data["transporter_inhibited"], expected_data)


@pytest.fixture(scope="module")
def bridgedb_dataframe_genes():
    """Reusable sample Pandas DataFrame to be used as input for the tests with genes."""
    return pd.DataFrame(
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


@pytest.fixture(scope="module")
def bridgedb_dataframe_compounds():
    """Reusable sample Pandas DataFrame to be used as input for the tests with genes."""
    return pd.DataFrame(
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