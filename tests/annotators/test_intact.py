#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the IntAct annotator."""

import unittest
from unittest.mock import Mock

import pandas as pd

from pyBiodatafuse.annotators import intact
from pyBiodatafuse.constants import INTACT_COMPOUND_INTERACT_COL, INTACT_INTERACT_COL


class TestIntact(unittest.TestCase):
    """Test the IntAct class."""

    def test_get_interactions(self):
        """Test the get_interactions function with mocked API calls.

        Uses controlled mock data to verify:
        - gene-gene interactions where only one interactor is the input gene are kept (OR logic)
        - gene-compound interactions are kept under 'both' mode
        - interactions where neither interactor is the input gene are excluded
        """
        intact.check_endpoint_intact = Mock(return_value=True)
        # DAG1 AC is EBI-1755945
        intact.get_protein_intact_acs = Mock(return_value=["EBI-1755945"])

        mock_interactions = [
            # Gene-gene: both sides are DAG1 (self-interaction) – should be kept
            {
                "interaction_id": "EBI-SELF",
                "interactor_id_A": "EBI-1755945",
                "interactor_id_B": "EBI-1755945",
                "score": 0.56,
                "biological_role_A": "unspecified role",
                "biological_role_B": "unspecified role",
                "type": "direct interaction",
                "detection_method": "x-ray diffraction",
                "host_organism": "In vitro",
                "interactor_A_name": "dag1_human",
                "interactor_B_name": "dag1_human",
                "interactor_A_species": "Homo sapiens",
                "interactor_B_species": "Homo sapiens",
                "molecule_A": "DAG1",
                "molecule_B": "DAG1",
                "id_A": "uniprotkb:Q14118",
                "id_B": "uniprotkb:Q14118",
                "pubmed_publication_id": "11111111",
            },
            # Gene-gene: only B is DAG1 – should be kept (OR logic)
            {
                "interaction_id": "EBI-PARTNER",
                "interactor_id_A": "EBI-9999999",
                "interactor_id_B": "EBI-1755945",
                "score": 0.4,
                "biological_role_A": "unspecified role",
                "biological_role_B": "unspecified role",
                "type": "association",
                "detection_method": "anti tag coip",
                "host_organism": "Homo sapiens HEK293T",
                "interactor_A_name": "partner_human",
                "interactor_B_name": "dag1_human",
                "interactor_A_species": "Homo sapiens",
                "interactor_B_species": "Homo sapiens",
                "molecule_A": "PARTNER",
                "molecule_B": "DAG1",
                "id_A": "uniprotkb:P99999",
                "id_B": "uniprotkb:Q14118",
                "pubmed_publication_id": "22222222",
            },
            # Gene-compound: A is a ChEBI compound, B is DAG1 – kept under 'both'
            {
                "interaction_id": "EBI-COMPOUND",
                "interactor_id_A": "EBI-5327879",
                "interactor_id_B": "EBI-1755945",
                "score": 0.4,
                "biological_role_A": "unspecified role",
                "biological_role_B": "unspecified role",
                "type": "physical association",
                "detection_method": "biophysical",
                "host_organism": "In vitro",
                "interactor_A_name": "ganglioside_gm1",
                "interactor_B_name": "dag1_human",
                "interactor_A_species": "Chemical synthesis (Chemical synthesis)",
                "interactor_B_species": "Homo sapiens",
                "molecule_A": "ganglioside_gm1",
                "molecule_B": "DAG1",
                "id_A": "CHEBI:61048",
                "id_B": "uniprotkb:Q14118",
                "pubmed_publication_id": "33333333",
            },
            # Gene-gene: neither A nor B is DAG1 – should be excluded
            {
                "interaction_id": "EBI-UNRELATED",
                "interactor_id_A": "EBI-1111111",
                "interactor_id_B": "EBI-2222222",
                "score": 0.5,
                "biological_role_A": "unspecified role",
                "biological_role_B": "unspecified role",
                "type": "association",
                "detection_method": "anti tag coip",
                "host_organism": "In vitro",
                "interactor_A_name": "other_a",
                "interactor_B_name": "other_b",
                "interactor_A_species": "Homo sapiens",
                "interactor_B_species": "Homo sapiens",
                "molecule_A": "OTHERA",
                "molecule_B": "OTHERB",
                "id_A": "uniprotkb:P11111",
                "id_B": "uniprotkb:P22222",
                "pubmed_publication_id": "44444444",
            },
        ]
        intact.get_intact_interactions = Mock(return_value=mock_interactions)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["DAG1"],
                "identifier.source": ["HGNC"],
                "target": ["ENSG00000173402"],
                "target.source": ["Ensembl"],
            }
        )

        obtained_data, metadata = intact.get_gene_interactions(
            bridgedb_dataframe, interaction_type="both"
        )

        result = obtained_data[INTACT_INTERACT_COL].iloc[0]
        # 3 interactions should be kept (self, partner, compound); unrelated excluded
        self.assertEqual(len(result), 3)
        interaction_ids = {r["interaction_id"] for r in result}
        self.assertIn("EBI-SELF", interaction_ids)
        self.assertIn("EBI-PARTNER", interaction_ids)
        self.assertIn("EBI-COMPOUND", interaction_ids)
        self.assertNotIn("EBI-UNRELATED", interaction_ids)
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

        obtained_data, metadata = intact.get_compound_interactions(
            bridgedb_dataframe, interaction_type="both"
        )

        expected_data = pd.Series(
            [
                [
                    {
                        "interaction_id": "EBI-9301798",
                        "interactor_id_A": "EBI-9096",
                        "interactor_id_B": "EBI-6621808",
                        "score": 0.44,
                        "biological_role_A": "enzyme",
                        "biological_role_B": "enzyme target",
                        "type": "enzymatic reaction",
                        "detection_method": "enzymatic study",
                        "host_organism": "In vitro",
                        "interactor_A_name": "ilvb_yeast",
                        "interactor_B_name": "pyruvate",
                        "interactor_A_species": "Saccharomyces cerevisiae",
                        "interactor_B_species": "Chemical synthesis (Chemical synthesis)",
                        "molecule_A": "ILV2",
                        "molecule_B": "pyruvate",
                        "id_A": "uniprotkb:P07342",
                        "id_B": "CHEBI:15361",
                        "pubmed_publication_id": "16390333",
                        "intact_link_to": "EBI-9096",
                    },
                    {
                        "interaction_id": "EBI-6621805",
                        "interactor_id_A": "EBI-372327",
                        "interactor_id_B": "EBI-6621808",
                        "score": 0.44,
                        "biological_role_A": "enzyme",
                        "biological_role_B": "enzyme target",
                        "type": "enzymatic reaction",
                        "detection_method": "oxidoreduct assay",
                        "host_organism": "In vitro",
                        "interactor_A_name": "ldha_human",
                        "interactor_B_name": "pyruvate",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Chemical synthesis (Chemical synthesis)",
                        "molecule_A": "LDHA",
                        "molecule_B": "pyruvate",
                        "id_A": "uniprotkb:P00338",
                        "id_B": "CHEBI:15361",
                        "pubmed_publication_id": "23523103",
                        "intact_link_to": "EBI-372327",
                    },
                ]
            ]
        )
        expected_data.name = INTACT_COMPOUND_INTERACT_COL

        pd.testing.assert_series_equal(obtained_data[INTACT_COMPOUND_INTERACT_COL], expected_data)
        self.assertIsInstance(metadata, dict)
