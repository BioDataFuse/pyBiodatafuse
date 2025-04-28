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
    get_interactions,
)
from pyBiodatafuse.constants import INTACT_COMPOUND_INTERACT_COL, INTACT_INTERACT_COL


class TestIntact(unittest.TestCase):
    """Test the IntAct class."""

    def test_get_interactions(self):
        """Test the get_interactions function."""
        intact.check_endpoint_intact = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["TP53", "MDM2"],
                "identifier.source": ["HGNC", "HGNC"],
                "target": ["ENSG00000141510", "ENSG00000135679"],
                "target.source": ["Ensembl", "Ensembl"],
            }
        )

        obtained_data, metadata = get_interactions(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "interaction_id": "EBI-21879744",
                        "interactor_id_A": "EBI-3895853",
                        "interactor_id_B": "EBI-5279149",
                        "binary_interaction_id": 4493009,
                        "confidence_values": ["intact-miscore:0.35"],
                        "intact_score": 0.35,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "anti tag coip",
                        "detection_method_id": "MI:0007",
                        "host_organism": "Homo sapiens HEK293T embryonic kidney cell",
                        "interactor_A_name": "p04637-2",
                        "interactor_B_name": "q00987-11",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "TP53",
                        "molecule_B": "MDM2",
                        "id_A": "P04637-2",
                        "id_B": "Q00987-11",
                        "pubmed_publication_id": "28514442",
                        "ensembl": "ENSG00000141510",
                        "intact_link_to": "MDM2",
                    },
                ],
                [
                    {
                        "interaction_id": "EBI-1551596",
                        "interactor_id_A": "EBI-389668",
                        "interactor_id_B": "EBI-389668",
                        "binary_interaction_id": 2180886,
                        "confidence_values": ["intact-miscore:0.74"],
                        "intact_score": 0.74,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "physical association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "anti tag coip",
                        "detection_method_id": "MI:0007",
                        "host_organism": "Homo sapiens HEK293 embryonic kidney cell",
                        "interactor_A_name": "mdm2_human",
                        "interactor_B_name": "mdm2_human",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "MDM2",
                        "molecule_B": "MDM2",
                        "id_A": "Q00987",
                        "id_B": "Q00987",
                        "pubmed_publication_id": "17936559",
                        "ensembl": "ENSG00000135679",
                        "intact_link_to": "MDM2",
                    },
                    {
                        "interaction_id": "EBI-4307418",
                        "interactor_id_A": "EBI-389668",
                        "interactor_id_B": "EBI-389668",
                        "binary_interaction_id": 9859622,
                        "confidence_values": ["intact-miscore:0.74"],
                        "intact_score": 0.74,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "physical association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "two hybrid array",
                        "detection_method_id": "MI:0397",
                        "host_organism": "Saccharomyces cerevisiae (Baker's yeast)",
                        "interactor_A_name": "mdm2_human",
                        "interactor_B_name": "mdm2_human",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "MDM2",
                        "molecule_B": "MDM2",
                        "id_A": "Q00987",
                        "id_B": "Q00987",
                        "pubmed_publication_id": "22493164",
                        "ensembl": "ENSG00000135679",
                        "intact_link_to": "MDM2",
                    },
                    {
                        "interaction_id": "EBI-7156290",
                        "interactor_id_A": "EBI-389668",
                        "interactor_id_B": "EBI-389668",
                        "binary_interaction_id": 12259773,
                        "confidence_values": ["intact-miscore:0.74"],
                        "intact_score": 0.74,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "physical association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "anti tag coip",
                        "detection_method_id": "MI:0007",
                        "host_organism": "Homo sapiens HEK293 embryonic kidney cell",
                        "interactor_A_name": "mdm2_human",
                        "interactor_B_name": "mdm2_human",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "MDM2",
                        "molecule_B": "MDM2",
                        "id_A": "Q00987",
                        "id_B": "Q00987",
                        "pubmed_publication_id": "17159902",
                        "ensembl": "ENSG00000135679",
                        "intact_link_to": "MDM2",
                    },
                    {
                        "interaction_id": "EBI-7156376",
                        "interactor_id_A": "EBI-389668",
                        "interactor_id_B": "EBI-389668",
                        "binary_interaction_id": 12259846,
                        "confidence_values": ["intact-miscore:0.74"],
                        "intact_score": 0.74,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "physical association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "anti bait coip",
                        "detection_method_id": "MI:0006",
                        "host_organism": "Homo sapiens HEK293 embryonic kidney cell",
                        "interactor_A_name": "mdm2_human",
                        "interactor_B_name": "mdm2_human",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "MDM2",
                        "molecule_B": "MDM2",
                        "id_A": "Q00987",
                        "id_B": "Q00987",
                        "pubmed_publication_id": "17159902",
                        "ensembl": "ENSG00000135679",
                        "intact_link_to": "MDM2",
                    },
                    {
                        "interaction_id": "EBI-21879744",
                        "interactor_id_A": "EBI-3895853",
                        "interactor_id_B": "EBI-5279149",
                        "binary_interaction_id": 4493009,
                        "confidence_values": ["intact-miscore:0.35"],
                        "intact_score": 0.35,
                        "biological_role_A": "unspecified role",
                        "biological_role_B": "unspecified role",
                        "type": "association",
                        "stoichiometry_A": "0-0",
                        "stoichiometry_B": "0-0",
                        "detection_method": "anti tag coip",
                        "detection_method_id": "MI:0007",
                        "host_organism": "Homo sapiens HEK293T embryonic kidney cell",
                        "interactor_A_name": "p04637-2",
                        "interactor_B_name": "q00987-11",
                        "interactor_A_species": "Homo sapiens",
                        "interactor_B_species": "Homo sapiens",
                        "molecule_A": "TP53",
                        "molecule_B": "MDM2",
                        "id_A": "P04637-2",
                        "id_B": "Q00987-11",
                        "pubmed_publication_id": "28514442",
                        "ensembl": "ENSG00000135679",
                        "intact_link_to": "TP53",
                    },
                ],
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
                "identifier": ["DAG1"],
                "identifier.source": ["HGNC"],
                "target": ["ENSG00000173402"],
                "target.source": ["Ensembl"],
            }
        )

        obtained_data, metadata = get_compound_interactions(bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "interaction_id": "EBI-5327885",
                        "interactor_id_A": "EBI-5327879",
                        "interactor_id_B": "EBI-1755945",
                        "binary_interaction_id": 10324525,
                        "confidence_values": ["intact-miscore:0.4"],
                        "intact_score": 0.4,
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
                        "ensembl": "ENSG00000173402",
                    },
                ]
            ]
        )
        expected_data.name = INTACT_COMPOUND_INTERACT_COL

        pd.testing.assert_series_equal(obtained_data[INTACT_COMPOUND_INTERACT_COL], expected_data)
        self.assertIsInstance(metadata, dict)
