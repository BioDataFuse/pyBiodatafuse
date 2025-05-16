#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the IntAct annotator."""

import unittest
from unittest.mock import Mock

import pandas as pd

from pyBiodatafuse.annotators import intact
from pyBiodatafuse.annotators.intact import (
    check_endpoint_intact,
    get_compound_interactions,
    get_gene_interactions,
)
from pyBiodatafuse.constants import INTACT_COMPOUND_INTERACT_COL, INTACT_INTERACT_COL


class TestIntact(unittest.TestCase):
    """Test the IntAct class."""

    def test_get_interactions(self):
        """Test the get_interactions function."""
        intact.check_endpoint_intact = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["DAG1"],
                "identifier.source": ["HGNC"],
                "target": ["ENSG00000173402"],
                "target.source": ["Ensembl"],
            }
        )

        obtained_data, metadata = get_gene_interactions(bridgedb_dataframe, interaction_type="both")

        expected_data = pd.Series(
            [
                [
                    {
                        "interaction_id": "EBI-7882257",
                        "interactor_id_A": "EBI-1755945",
                        "interactor_id_B": "EBI-1755945",
                        "binary_interaction_id": 13027532,
                        "confidence_values": ["intact-miscore:0.56"],
                        "score": 0.56,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "direct interaction",
                        "stoichiometry_A": "1-1",
                        "stoichiometry_B": "1-1",
                        "detection_method": "x-ray diffraction",
                        "detection_method_id": "MI:0114",
                        "host_organism": "In vitro",
                        "interactor_A_name": "dag1_human",
                        "interactor_B_name": "dag1_human",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "DAG1",
                        "molecule_B": "DAG1",
                        "id_A": "Q14118",
                        "id_B": "Q14118",
                        "pubmed_publication_id": "11423118",
                        "altIdsA": [
                            "Q969J9 (uniprotkb)",
                            "ENSP00000412067.2 (ensembl)",
                            "A8K6M7 (uniprotkb)",
                            "ENSP00000387859.2 (ensembl)",
                            "ENSP00000501140.1 (ensembl)",
                            "Q14118 (uniprotkb)",
                            "ENSP00000513218.1 (ensembl)",
                            "EBI-1755945 (intact)",
                            "ENSP00000405859.2 (ensembl)",
                            "ENSP00000401805.2 (ensembl)",
                            "ENSP00000501165.2 (ensembl)",
                            "ENSP00000401382.3 (ensembl)",
                            "ENSP00000415321.4 (ensembl)",
                            "ENSP00000388833.2 (ensembl)",
                            "ENSP00000410145.3 (ensembl)",
                            "ENSP00000513217.1 (ensembl)",
                            "ENSP00000513216.1 (ensembl)",
                            "ENSP00000312435.2 (ensembl)",
                        ],
                        "altIdsB": [
                            "Q969J9 (uniprotkb)",
                            "ENSP00000412067.2 (ensembl)",
                            "A8K6M7 (uniprotkb)",
                            "ENSP00000387859.2 (ensembl)",
                            "ENSP00000501140.1 (ensembl)",
                            "Q14118 (uniprotkb)",
                            "ENSP00000513218.1 (ensembl)",
                            "EBI-1755945 (intact)",
                            "ENSP00000405859.2 (ensembl)",
                            "ENSP00000401805.2 (ensembl)",
                            "ENSP00000501165.2 (ensembl)",
                            "ENSP00000401382.3 (ensembl)",
                            "ENSP00000415321.4 (ensembl)",
                            "ENSP00000388833.2 (ensembl)",
                            "ENSP00000410145.3 (ensembl)",
                            "ENSP00000513217.1 (ensembl)",
                            "ENSP00000513216.1 (ensembl)",
                            "ENSP00000312435.2 (ensembl)",
                        ],
                        "intact_link_to": "DAG1",
                    },
                    {
                        "interaction_id": "EBI-7882311",
                        "interactor_id_A": "EBI-1755945",
                        "interactor_id_B": "EBI-1755945",
                        "binary_interaction_id": 13027607,
                        "confidence_values": ["intact-miscore:0.56"],
                        "score": 0.56,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "direct interaction",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "elisa",
                        "detection_method_id": "MI:0411",
                        "host_organism": "In vitro",
                        "interactor_A_name": "dag1_human",
                        "interactor_B_name": "dag1_human",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "DAG1",
                        "molecule_B": "DAG1",
                        "id_A": "Q14118",
                        "id_B": "Q14118",
                        "pubmed_publication_id": "11423118",
                        "altIdsA": [
                            "Q969J9 (uniprotkb)",
                            "ENSP00000412067.2 (ensembl)",
                            "A8K6M7 (uniprotkb)",
                            "ENSP00000387859.2 (ensembl)",
                            "ENSP00000501140.1 (ensembl)",
                            "Q14118 (uniprotkb)",
                            "ENSP00000513218.1 (ensembl)",
                            "EBI-1755945 (intact)",
                            "ENSP00000405859.2 (ensembl)",
                            "ENSP00000401805.2 (ensembl)",
                            "ENSP00000501165.2 (ensembl)",
                            "ENSP00000401382.3 (ensembl)",
                            "ENSP00000415321.4 (ensembl)",
                            "ENSP00000388833.2 (ensembl)",
                            "ENSP00000410145.3 (ensembl)",
                            "ENSP00000513217.1 (ensembl)",
                            "ENSP00000513216.1 (ensembl)",
                            "ENSP00000312435.2 (ensembl)",
                        ],
                        "altIdsB": [
                            "Q969J9 (uniprotkb)",
                            "ENSP00000412067.2 (ensembl)",
                            "A8K6M7 (uniprotkb)",
                            "ENSP00000387859.2 (ensembl)",
                            "ENSP00000501140.1 (ensembl)",
                            "Q14118 (uniprotkb)",
                            "ENSP00000513218.1 (ensembl)",
                            "EBI-1755945 (intact)",
                            "ENSP00000405859.2 (ensembl)",
                            "ENSP00000401805.2 (ensembl)",
                            "ENSP00000501165.2 (ensembl)",
                            "ENSP00000401382.3 (ensembl)",
                            "ENSP00000415321.4 (ensembl)",
                            "ENSP00000388833.2 (ensembl)",
                            "ENSP00000410145.3 (ensembl)",
                            "ENSP00000513217.1 (ensembl)",
                            "ENSP00000513216.1 (ensembl)",
                            "ENSP00000312435.2 (ensembl)",
                        ],
                        "intact_link_to": "DAG1",
                    },
                    {
                        "interaction_id": "EBI-5327885",
                        "interactor_id_A": "EBI-5327879",
                        "interactor_id_B": "EBI-1755945",
                        "binary_interaction_id": 10324525,
                        "confidence_values": ["intact-miscore:0.4"],
                        "score": 0.4,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "physical association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "biophysical",
                        "detection_method_id": "MI:0013",
                        "host_organism": "Homo sapiens HeLa S3 epitheloid cervical carcinoma cell",
                        "interactor_A_name": "ganglioside_gm1",
                        "interactor_B_name": "dag1_human",
                        "interactor_A_species": "Chemical synthesis (Chemical synthesis)",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "ganglioside_gm1",
                        "molecule_B": "DAG1",
                        "id_A": "CHEBI:61048",
                        "id_B": "Q14118",
                        "pubmed_publication_id": "22106087",
                        "altIdsA": ["EBI-5327879 (intact)", "CHEBI:61048 (chebi)"],
                        "altIdsB": [
                            "Q969J9 (uniprotkb)",
                            "ENSP00000412067.2 (ensembl)",
                            "A8K6M7 (uniprotkb)",
                            "ENSP00000387859.2 (ensembl)",
                            "ENSP00000501140.1 (ensembl)",
                            "Q14118 (uniprotkb)",
                            "ENSP00000513218.1 (ensembl)",
                            "EBI-1755945 (intact)",
                            "ENSP00000405859.2 (ensembl)",
                            "ENSP00000401805.2 (ensembl)",
                            "ENSP00000501165.2 (ensembl)",
                            "ENSP00000401382.3 (ensembl)",
                            "ENSP00000415321.4 (ensembl)",
                            "ENSP00000388833.2 (ensembl)",
                            "ENSP00000410145.3 (ensembl)",
                            "ENSP00000513217.1 (ensembl)",
                            "ENSP00000513216.1 (ensembl)",
                            "ENSP00000312435.2 (ensembl)",
                        ],
                        "intact_link_to": None,
                    },
                ]
            ]
        )
        expected_data.name = INTACT_INTERACT_COL

        pd.testing.assert_series_equal(obtained_data[INTACT_INTERACT_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    def test_get_compound_interactions(self):
        """Test the get_compound_interactions function."""
        intact.check_endpoint_intact = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["15361"],
                "identifier.source": ["ChEBI"],
                "target": ["CHEBI:15361"],
                "target.source": ["ChEBI"],
            }
        )

        obtained_data, metadata = get_compound_interactions(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "interaction_id": "EBI-9301798",
                        "interactor_id_A": "EBI-9096",
                        "interactor_id_B": "EBI-6621808",
                        "binary_interaction_id": 13894862,
                        "confidence_values": ["intact-miscore:0.44"],
                        "score": 0.44,
                        "biological_role_A": "enzyme",
                        "biological_role_B": "enzyme target",
                        "type": "enzymatic reaction",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "enzymatic study",
                        "detection_method_id": "MI:0415",
                        "host_organism": "In vitro",
                        "interactor_A_name": "ilvb_yeast",
                        "interactor_B_name": "pyruvate",
                        "interactor_A_species": "Saccharomyces cerevisiae",
                        "interactor_B_species": "Chemical synthesis (Chemical synthesis)",
                        "molecule_A": "ILV2",
                        "molecule_B": "pyruvate",
                        "id_A": "P07342",
                        "id_B": "CHEBI:15361",
                        "pubmed_publication_id": "16390333",
                        "altIdsA": [
                            "P07342 (uniprotkb)",
                            "D6VZT1 (uniprotkb)",
                            "EBI-9096 (intact)",
                        ],
                        "altIdsB": ["CHEBI:15361 (chebi)", "EBI-6621808 (intact)"],
                        "intact_link_to": None,
                    },
                    {
                        "interaction_id": "EBI-6621805",
                        "interactor_id_A": "EBI-372327",
                        "interactor_id_B": "EBI-6621808",
                        "binary_interaction_id": 11900151,
                        "confidence_values": ["intact-miscore:0.44"],
                        "score": 0.44,
                        "biological_role_A": "enzyme",
                        "biological_role_B": "enzyme target",
                        "type": "enzymatic reaction",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "oxidoreduct assay",
                        "detection_method_id": "MI:0979",
                        "host_organism": "In vitro",
                        "interactor_A_name": "ldha_human",
                        "interactor_B_name": "pyruvate",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Chemical synthesis (Chemical synthesis)",
                        "molecule_A": "LDHA",
                        "molecule_B": "pyruvate",
                        "id_A": "P00338",
                        "id_B": "CHEBI:15361",
                        "pubmed_publication_id": "23523103",
                        "altIdsA": [
                            "ENSP00000395337.3 (ensembl)",
                            "Q53G53 (uniprotkb)",
                            "EBI-372327 (intact)",
                            "D3DQY3 (uniprotkb)",
                            "F8W819 (uniprotkb)",
                            "B7Z5E3 (uniprotkb)",
                            "ENSP00000445331.1 (ensembl)",
                            "ENSP00000500953.1 (ensembl)",
                            "Q9UDE9 (uniprotkb)",
                            "Q9UDE8 (uniprotkb)",
                            "Q6IBM7 (uniprotkb)",
                            "ENSP00000499898.1 (ensembl)",
                            "Q6ZNV1 (uniprotkb)",
                            "P00338 (uniprotkb)",
                            "B4DKQ2 (uniprotkb)",
                            "ENSP00000499977.1 (ensembl)",
                        ],
                        "altIdsB": ["CHEBI:15361 (chebi)", "EBI-6621808 (intact)"],
                        "intact_link_to": None,
                    },
                ]
            ]
        )
        expected_data.name = INTACT_COMPOUND_INTERACT_COL

        pd.testing.assert_series_equal(obtained_data[INTACT_COMPOUND_INTERACT_COL], expected_data)
        self.assertIsInstance(metadata, dict)
