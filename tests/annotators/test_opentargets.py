#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenTargets annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
from numpy import nan

from pyBiodatafuse import id_mapper
from pyBiodatafuse.annotators import opentargets
from pyBiodatafuse.annotators.opentargets import (
    get_disease_compound_interactions,
    get_gene_compound_interactions,
    get_gene_go_process,
    get_gene_reactome_pathways,
    get_version_opentargets,
)
from pyBiodatafuse.constants import (
    OPENTARGETS_DISEASE_COMPOUND_COL,
    OPENTARGETS_GENE_COMPOUND_COL,
    OPENTARGETS_GO_COL,
    OPENTARGETS_REACTOME_COL,
)

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestOpentarget(unittest.TestCase):
    """Test the Opentarget class."""

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_version_opentargets(self, mock_post_version):
        """Test that the GraphQL endpoint returns the expected results."""
        mock_post_version.return_value.json.return_value = {
            "data": {
                "meta": {
                    "name": "Open Targets GraphQL & REST API Beta",
                    "apiVersion": {"x": "24", "y": "0", "z": "3"},
                    "dataVersion": {"year": "24", "month": "03"},
                }
            }
        }

        obtained_version = get_version_opentargets()

        expected_version = {
            "datasource": "Open Targets GraphQL & REST API Beta",
            "metadata": {
                "source_version": {"apiVersion": "24.0.3"},
                "data_version": "24-03",
            },
        }

        assert obtained_version == expected_version

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_gene_go_process(self, mock_post_go_process):
        """Test that the GraphQL endpoint returns the expected results."""
        opentargets.check_endpoint_opentargets = Mock(return_value=True)
        opentargets.get_version_opentargets = Mock(
            return_value={
                "data": {
                    "meta": {
                        "name": "Open Targets GraphQL & REST API Beta",
                        "apiVersion": {"x": "24", "y": "0", "z": "3"},
                        "dataVersion": {"year": "24", "month": "03"},
                    }
                }
            }
        )

        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "identifier": ["ALG14", "ALG2", "CHRNA1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["ENSG00000172339", "ENSG00000119523", "ENSG00000138435"],
                "target.source": ["Ensembl", "Ensembl", "Ensembl"],
            }
        )

        with open(os.path.join(data_file_folder, "opentargets_go_mock_data.json")) as f:
            mock_post_go_process.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_gene_go_process(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "go_id": "GO:0043541",
                        "go_name": "UDP-N-acetylglucosamine transferase complex",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0043495",
                        "go_name": "protein-membrane adaptor activity",
                        "go_type": "F",
                    },
                    {
                        "go_id": "GO:0098554",
                        "go_name": "cytoplasmic side of endoplasmic reticulum membrane",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0005789",
                        "go_name": "endoplasmic reticulum membrane",
                        "go_type": "C",
                    },
                    {"go_id": "GO:0005515", "go_name": "protein binding", "go_type": "F"},
                    {
                        "go_id": "GO:0006487",
                        "go_name": "protein N-linked glycosylation",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0006488",
                        "go_name": "dolichol-linked oligosaccharide biosynthetic process",
                        "go_type": "P",
                    },
                ],
                [
                    {"go_id": "GO:0012505", "go_name": "endomembrane system", "go_type": "C"},
                    {
                        "go_id": "GO:0000033",
                        "go_name": "alpha-1,3-mannosyltransferase activity",
                        "go_type": "F",
                    },
                    {
                        "go_id": "GO:0005789",
                        "go_name": "endoplasmic reticulum membrane",
                        "go_type": "C",
                    },
                    {"go_id": "GO:0006486", "go_name": "protein glycosylation", "go_type": "P"},
                    {
                        "go_id": "GO:0004378",
                        "go_name": "GDP-Man:Man1GlcNAc2-PP-Dol alpha-1,3-mannosyltransferase activity",
                        "go_type": "F",
                    },
                    {
                        "go_id": "GO:0006488",
                        "go_name": "dolichol-linked oligosaccharide biosynthetic process",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0102704",
                        "go_name": "GDP-Man:Man2GlcNAc2-PP-dolichol alpha-1,6-mannosyltransferase activity",
                        "go_type": "F",
                    },
                    {"go_id": "GO:0016020", "go_name": "membrane", "go_type": "C"},
                    {
                        "go_id": "GO:0098554",
                        "go_name": "cytoplasmic side of endoplasmic reticulum membrane",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0006487",
                        "go_name": "protein N-linked glycosylation",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0005515", "go_name": "protein binding", "go_type": "F"},
                ],
                [
                    {
                        "go_id": "GO:0007274",
                        "go_name": "neuromuscular synaptic transmission",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0050905", "go_name": "neuromuscular process", "go_type": "P"},
                    {"go_id": "GO:0009986", "go_name": "cell surface", "go_type": "C"},
                    {
                        "go_id": "GO:0070050",
                        "go_name": "neuron cellular homeostasis",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0043005", "go_name": "neuron projection", "go_type": "C"},
                    {
                        "go_id": "GO:0048630",
                        "go_name": "skeletal muscle tissue growth",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0035094", "go_name": "response to nicotine", "go_type": "P"},
                    {
                        "go_id": "GO:0042391",
                        "go_name": "regulation of membrane potential",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0005886", "go_name": "plasma membrane", "go_type": "C"},
                    {
                        "go_id": "GO:0005892",
                        "go_name": "acetylcholine-gated channel complex",
                        "go_type": "C",
                    },
                    {"go_id": "GO:0031594", "go_name": "neuromuscular junction", "go_type": "C"},
                    {"go_id": "GO:0051899", "go_name": "membrane depolarization", "go_type": "P"},
                    {"go_id": "GO:0050881", "go_name": "musculoskeletal movement", "go_type": "P"},
                    {
                        "go_id": "GO:0099634",
                        "go_name": "postsynaptic specialization membrane",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0022848",
                        "go_name": "acetylcholine-gated monoatomic cation-selective channel activity",
                        "go_type": "F",
                    },
                    {"go_id": "GO:0019228", "go_name": "neuronal action potential", "go_type": "P"},
                    {
                        "go_id": "GO:1904315",
                        "go_name": "transmitter-gated monoatomic ion channel activity involved in regulation of postsynaptic membrane potential",
                        "go_type": "F",
                    },
                    {
                        "go_id": "GO:0003009",
                        "go_name": "skeletal muscle contraction",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0045211", "go_name": "postsynaptic membrane", "go_type": "C"},
                    {
                        "go_id": "GO:0007271",
                        "go_name": "synaptic transmission, cholinergic",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0046716",
                        "go_name": "muscle cell cellular homeostasis",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0045202", "go_name": "synapse", "go_type": "C"},
                    {
                        "go_id": "GO:0060079",
                        "go_name": "excitatory postsynaptic potential",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0015464",
                        "go_name": "acetylcholine receptor activity",
                        "go_type": "F",
                    },
                    {"go_id": "GO:0042166", "go_name": "acetylcholine binding", "go_type": "F"},
                    {
                        "go_id": "GO:0006812",
                        "go_name": "monoatomic cation transport",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0007528",
                        "go_name": "neuromuscular junction development",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0095500",
                        "go_name": "acetylcholine receptor signaling pathway",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0034220",
                        "go_name": "monoatomic ion transmembrane transport",
                        "go_type": "P",
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_GO_COL

        pd.testing.assert_series_equal(obtained_data[OPENTARGETS_GO_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_gene_reactome_pathways(self, mock_post_reactome):
        """Test that the GraphQL endpoint returns the expected results."""
        opentargets.check_endpoint_opentargets = Mock(return_value=True)
        opentargets.get_version_opentargets = Mock(
            return_value={
                "data": {
                    "meta": {
                        "name": "Open Targets GraphQL & REST API Beta",
                        "apiVersion": {"x": "24", "y": "0", "z": "3"},
                        "dataVersion": {"year": "24", "month": "03"},
                    }
                }
            }
        )

        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "identifier": ["ALG14", "ALG2", "CHRNA1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["ENSG00000172339", "ENSG00000119523", "ENSG00000138435"],
                "target.source": ["Ensembl", "Ensembl", "Ensembl"],
            }
        )

        with open(os.path.join(data_file_folder, "opentargets_reactome_mock_data.json")) as f:
            mock_post_reactome.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_gene_reactome_pathways(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "pathway_label": "Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein",
                        "pathway_id": "Reactome:R-HSA-446193",
                    },
                    {
                        "pathway_label": "Defective ALG14 causes ALG14-CMS",
                        "pathway_id": "Reactome:R-HSA-5633231",
                    },
                ],
                [
                    {
                        "pathway_label": "Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein",
                        "pathway_id": "Reactome:R-HSA-446193",
                    },
                    {
                        "pathway_label": "Defective ALG2 causes CDG-1i",
                        "pathway_id": "Reactome:R-HSA-4549349",
                    },
                ],
                [
                    {
                        "pathway_label": "Highly calcium permeable nicotinic acetylcholine receptors",
                        "pathway_id": "Reactome:R-HSA-629597",
                    },
                    {
                        "pathway_label": "Highly calcium permeable postsynaptic nicotinic acetylcholine receptors",
                        "pathway_id": "Reactome:R-HSA-629594",
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_REACTOME_COL

        pd.testing.assert_series_equal(obtained_data[OPENTARGETS_REACTOME_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_gene_compound_interactions(self, mock_post_compound):
        """Test that the GraphQL endpoint returns the expected results."""
        opentargets.check_endpoint_opentargets = Mock(return_value=True)
        opentargets.get_version_opentargets = Mock(
            return_value={
                "data": {
                    "meta": {
                        "name": "Open Targets GraphQL & REST API Beta",
                        "apiVersion": {"x": "24", "y": "0", "z": "3"},
                        "dataVersion": {"year": "24", "month": "03"},
                    }
                }
            }
        )

        id_mapper.pubchem_xref = Mock(
            return_value=(
                pd.DataFrame(
                    {
                        "identifier": [
                            "CHEMBL1200641",
                            "CHEMBL1201244",
                            "CHEMBL1200648",
                            "CHEMBL1201248",
                            "CHEMBL703",
                            "CHEMBL1182833",
                            "CHEMBL984",
                        ],
                        "identifier.source": [
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                        ],
                        "target": [
                            "CID:62886",
                            "CID:441290",
                            "CID:441351",
                            "CID:62887",
                            "CID:5314",
                            "CID:5281042",
                            "CID:5281080",
                        ],
                        "target.source": [
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                        ],
                    }
                ),
                {},
            )
        )

        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "identifier": ["ALG14", "ALG2", "CHRNA1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["ENSG00000172339", "ENSG00000119523", "ENSG00000138435"],
                "target.source": ["Ensembl", "Ensembl", "Ensembl"],
            }
        )

        with open(os.path.join(data_file_folder, "opentargets_gene_compound_mock_data.json")) as f:
            mock_post_compound.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_gene_compound_interactions(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "chembl_id": nan,
                        "drugbank_id": nan,
                        "compound_cid": nan,
                        "compound_name": nan,
                        "clincal_trial_phase": nan,
                        "is_approved": nan,
                        "relation": nan,
                        "adverse_effect_count": nan,
                        "adverse_effect": nan,
                    }
                ],
                [
                    {
                        "chembl_id": nan,
                        "drugbank_id": nan,
                        "compound_cid": nan,
                        "compound_name": nan,
                        "clincal_trial_phase": nan,
                        "is_approved": nan,
                        "relation": nan,
                        "adverse_effect_count": nan,
                        "adverse_effect": nan,
                    }
                ],
                [
                    {
                        "chembl_id": "CHEMBL:CHEMBL1200641",
                        "drugbank_id": None,
                        "compound_cid": "CID:62886",
                        "compound_name": "CISATRACURIUM BESYLATE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 100.0,
                        "adverse_effect": [
                            {"name": "anaphylactic shock"},
                            {"name": "anaphylactic reaction"},
                            {"name": "feeding intolerance"},
                            {"name": "bronchospasm"},
                            {"name": "hyperthermia malignant"},
                            {"name": "circulatory collapse"},
                            {"name": "hypotension"},
                            {"name": "cardiac arrest"},
                            {"name": "hypertransaminasaemia"},
                            {"name": "cardio-respiratory arrest"},
                            {"name": "drug reaction with eosinophilia and systemic symptoms"},
                            {"name": "shock"},
                            {"name": "apnoea"},
                            {"name": "toxic epidermal necrolysis"},
                            {"name": "post procedural complication"},
                            {"name": "acute pulmonary oedema"},
                            {"name": "haemodynamic instability"},
                            {"name": "fatigue"},
                            {"name": "hepatocellular injury"},
                            {"name": "hyperbilirubinaemia"},
                            {"name": "bradycardia"},
                            {"name": "diarrhoea"},
                            {"name": "tachyphylaxis"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "headache"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201244",
                        "drugbank_id": "DrugBank:DB00728",
                        "compound_cid": "CID:441290",
                        "compound_name": "ROCURONIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 104.0,
                        "adverse_effect": [
                            {"name": "anaphylactic reaction"},
                            {"name": "anaphylactic shock"},
                            {"name": "cardiac arrest"},
                            {"name": "hypotension"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "hyperthermia malignant"},
                            {"name": "bronchospasm"},
                            {"name": "circulatory collapse"},
                            {"name": "tachycardia"},
                            {"name": "delayed recovery from anaesthesia"},
                            {"name": "neuromuscular blockade"},
                            {"name": "recurrence of neuromuscular blockade"},
                            {"name": "pulseless electrical activity"},
                            {"name": "bradycardia"},
                            {"name": "wrist surgery"},
                            {"name": "finger deformity"},
                            {"name": "tenosynovitis"},
                            {"name": "stress cardiomyopathy"},
                            {"name": "oxygen saturation decreased"},
                            {"name": "negative pressure pulmonary oedema"},
                            {"name": "airway peak pressure increased"},
                            {"name": "procedural hypotension"},
                            {"name": "hypoventilation"},
                            {"name": "fatigue"},
                            {"name": "hypoxia"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1200648",
                        "drugbank_id": None,
                        "compound_cid": "CID:441351",
                        "compound_name": "ROCURONIUM BROMIDE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 87.0,
                        "adverse_effect": [
                            {"name": "anaphylactic reaction"},
                            {"name": "anaphylactic shock"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "cardiac arrest"},
                            {"name": "hypotension"},
                            {"name": "bronchospasm"},
                            {"name": "hyperthermia malignant"},
                            {"name": "circulatory collapse"},
                            {"name": "tachycardia"},
                            {"name": "delayed recovery from anaesthesia"},
                            {"name": "oxygen saturation decreased"},
                            {"name": "recurrence of neuromuscular blockade"},
                            {"name": "neuromuscular blockade"},
                            {"name": "bradycardia"},
                            {"name": "pulseless electrical activity"},
                            {"name": "wrist surgery"},
                            {"name": "finger deformity"},
                            {"name": "endotracheal intubation complication"},
                            {"name": "tenosynovitis"},
                            {"name": "fatigue"},
                            {"name": "airway complication of anaesthesia"},
                            {"name": "stress cardiomyopathy"},
                            {"name": "procedural hypotension"},
                            {"name": "nausea"},
                            {"name": "haemodynamic instability"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201248",
                        "drugbank_id": "DrugBank:DB00565",
                        "compound_cid": "CID:62887",
                        "compound_name": "CISATRACURIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 96.0,
                        "adverse_effect": [
                            {"name": "anaphylactic shock"},
                            {"name": "anaphylactic reaction"},
                            {"name": "feeding intolerance"},
                            {"name": "circulatory collapse"},
                            {"name": "bronchospasm"},
                            {"name": "hyperthermia malignant"},
                            {"name": "hypotension"},
                            {"name": "cardiac arrest"},
                            {"name": "cardio-respiratory arrest"},
                            {"name": "apnoea"},
                            {"name": "shock"},
                            {"name": "post procedural complication"},
                            {"name": "haemodynamic instability"},
                            {"name": "bradycardia"},
                            {"name": "hepatocellular injury"},
                            {"name": "hypertransaminasaemia"},
                            {"name": "acute pulmonary oedema"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "fatigue"},
                            {"name": "hyperbilirubinaemia"},
                            {"name": "drug reaction with eosinophilia and systemic symptoms"},
                            {"name": "diarrhoea"},
                            {"name": "tachyphylaxis"},
                            {"name": "myocardial stunning"},
                            {"name": "rash maculo-papular"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL703",
                        "drugbank_id": "DrugBank:DB00202",
                        "compound_cid": "CID:5314",
                        "compound_name": "SUXAMETHONIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "activates",
                        "adverse_effect_count": 125.0,
                        "adverse_effect": [
                            {"name": "hyperthermia malignant"},
                            {"name": "anaphylactic shock"},
                            {"name": "cardiac arrest"},
                            {"name": "hypotension"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "anaphylactic reaction"},
                            {"name": "bronchospasm"},
                            {"name": "pseudocholinesterase deficiency"},
                            {"name": "premature baby"},
                            {"name": "anaesthetic complication"},
                            {"name": "rhabdomyolysis"},
                            {"name": "ventricular tachycardia"},
                            {"name": "pco2 increased"},
                            {"name": "circulatory collapse"},
                            {"name": "delayed recovery from anaesthesia"},
                            {"name": "apnoea"},
                            {"name": "pulseless electrical activity"},
                            {"name": "post procedural complication"},
                            {"name": "renal ischaemia"},
                            {"name": "diabetes insipidus"},
                            {"name": "cardio-respiratory arrest"},
                            {"name": "left ventricular dysfunction"},
                            {"name": "blood ph decreased"},
                            {"name": "foetal death"},
                            {"name": "muscle contractions involuntary"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1182833",
                        "drugbank_id": "DrugBank:DB01226",
                        "compound_cid": "CID:5281042",
                        "compound_name": "MIVACURIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 13.0,
                        "adverse_effect": [
                            {"name": "pseudocholinesterase deficiency"},
                            {"name": "blood pressure measurement"},
                            {"name": "performance status decreased"},
                            {"name": "weight bearing difficulty"},
                            {"name": "cardiovascular disorder"},
                            {"name": "oropharyngeal discomfort"},
                            {"name": "anaesthetic complication neurological"},
                            {"name": "hyperventilation"},
                            {"name": "impaired work ability"},
                            {"name": "pulseless electrical activity"},
                            {"name": "heart rate increased"},
                            {"name": "general physical health deterioration"},
                            {"name": "flatulence"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL984",
                        "drugbank_id": None,
                        "compound_cid": "CID:5281080",
                        "compound_name": "MIVACURIUM CHLORIDE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 12.0,
                        "adverse_effect": [
                            {"name": "pseudocholinesterase deficiency"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "blood pressure measurement"},
                            {"name": "acute generalised exanthematous pustulosis"},
                            {"name": "weight bearing difficulty"},
                            {"name": "anaesthetic complication neurological"},
                            {"name": "performance status decreased"},
                            {"name": "oropharyngeal discomfort"},
                            {"name": "hyperventilation"},
                            {"name": "impaired work ability"},
                            {"name": "anaphylactoid shock"},
                            {"name": "cardiovascular disorder"},
                        ],
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_GENE_COMPOUND_COL

        pd.testing.assert_series_equal(
            obtained_data[OPENTARGETS_GENE_COMPOUND_COL], expected_data, check_dtype=False
        )
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_disease_compound_interactions(self, mock_post_compound):
        """Test that the GraphQL endpoint returns the expected results."""
        opentargets.check_endpoint_opentargets = Mock(return_value=True)
        opentargets.get_version_opentargets = Mock(
            return_value={
                "data": {
                    "meta": {
                        "name": "Open Targets GraphQL & REST API Beta",
                        "apiVersion": {"x": "24", "y": "0", "z": "3"},
                        "dataVersion": {"year": "24", "month": "03"},
                    }
                }
            }
        )
        id_mapper.pubchem_xref = Mock(
            return_value=(
                pd.DataFrame(
                    {
                        "identifier": [
                            "CHEMBL384467",
                            "CHEMBL1504",
                            "CHEMBL655",
                            "CHEMBL778",
                            "CHEMBL277474",
                            "CHEMBL596",
                            "CHEMBL240163",
                            "CHEMBL1201231",
                            "CHEMBL1201014",
                            "CHEMBL1161",
                        ],
                        "identifier.source": [
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                            "name",
                        ],
                        "target": [
                            "CID:5743",
                            "CID:6436",
                            "CID:4192",
                            "CID:5311068",
                            "CID:2206",
                            "CID:3345",
                            "CID:65957",
                            "CID:72078",
                            "CID:441409",
                            "CID:441336",
                        ],
                        "target.source": [
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                            "PubChem Compound",
                        ],
                    }
                ),
                {},
            )
        )

        bridge_df = pd.DataFrame(
            {
                "identifier": ["A2ML1", "ABCA1", "A2ML1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": [
                    "EFO_0004992",
                    "MONDO:0005359",
                    "EFO_0700136",
                ],
                "target.source": [
                    "EFO",
                    "MONDO",
                    "EFO",
                ],
            }
        )

        with open(os.path.join(data_file_folder, "opentargets_disease_compound_mock.json")) as f:
            mock_post_compound.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_disease_compound_interactions(bridge_df)

        expected_data = pd.Series(
            [
                [
                    {
                        "chembl_id": "CHEMBL:CHEMBL384467",
                        "drugbank_id": "DrugBank:DB01234",
                        "compound_cid": "CID:5743",
                        "compound_name": "DEXAMETHASONE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 6.0,
                        "adverse_effect": [
                            {"name": "plasma cell myeloma"},
                            {"name": "febrile neutropenia"},
                            {"name": "neutropenia"},
                            {"name": "thrombocytopenia"},
                            {"name": "neuropathy peripheral"},
                            {"name": "febrile bone marrow aplasia"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1504",
                        "drugbank_id": None,
                        "compound_cid": "CID:6436",
                        "compound_name": "TRIAMCINOLONE ACETONIDE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 70.0,
                        "adverse_effect": [
                            {"name": "cushing's syndrome"},
                            {"name": "endophthalmitis"},
                            {"name": "non-infectious endophthalmitis"},
                            {"name": "adrenal insufficiency"},
                            {"name": "injection site atrophy"},
                            {"name": "macular degeneration"},
                            {"name": "intraocular pressure increased"},
                            {"name": "ocular hypertension"},
                            {"name": "skin atrophy"},
                            {"name": "visual acuity reduced"},
                            {"name": "eye infection toxoplasmal"},
                            {"name": "human antichimeric antibody positive"},
                            {"name": "joint destruction"},
                            {"name": "cytomegalovirus chorioretinitis"},
                            {"name": "wheelchair user"},
                            {"name": "body height below normal"},
                            {"name": "polyarthritis"},
                            {"name": "retinal artery occlusion"},
                            {"name": "eye inflammation"},
                            {"name": "maculopathy"},
                            {"name": "obesity"},
                            {"name": "adrenal suppression"},
                            {"name": "necrotising retinitis"},
                            {"name": "hypothalamic pituitary adrenal axis suppression"},
                            {"name": "glaucoma"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL655",
                        "drugbank_id": "DrugBank:DB00683",
                        "compound_cid": "CID:4192",
                        "compound_name": "MIDAZOLAM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 60.0,
                        "adverse_effect": [
                            {"name": "hypotension"},
                            {"name": "status epilepticus"},
                            {"name": "respiratory depression"},
                            {"name": "anaphylactic reaction"},
                            {"name": "delirium"},
                            {"name": "seizure"},
                            {"name": "oxygen saturation decreased"},
                            {"name": "sedation"},
                            {"name": "drug reaction with eosinophilia and systemic symptoms"},
                            {"name": "multiple-drug resistance"},
                            {"name": "tachycardia"},
                            {"name": "intestinal ischaemia"},
                            {"name": "cardiac arrest"},
                            {"name": "anaphylactic shock"},
                            {"name": "generalised tonic-clonic seizure"},
                            {"name": "anaesthetic complication neurological"},
                            {"name": "necrosis ischaemic"},
                            {"name": "hyperthermia malignant"},
                            {"name": "agitation"},
                            {"name": "unresponsive to stimuli"},
                            {"name": "systemic inflammatory response syndrome"},
                            {"name": "bradycardia"},
                            {"name": "apnoea"},
                            {"name": "upper airway obstruction"},
                            {"name": "myoclonus"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL778",
                        "drugbank_id": "DrugBank:DB00633",
                        "compound_cid": "CID:5311068",
                        "compound_name": "DEXMEDETOMIDINE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 80.0,
                        "adverse_effect": [
                            {"name": "delirium"},
                            {"name": "bradycardia"},
                            {"name": "hyperthermia"},
                            {"name": "agitation"},
                            {"name": "cardiac arrest"},
                            {"name": "withdrawal syndrome"},
                            {"name": "hypotension"},
                            {"name": "respiratory depression"},
                            {"name": "upper airway obstruction"},
                            {"name": "diabetes insipidus"},
                            {"name": "intestinal pseudo-obstruction"},
                            {"name": "acute motor axonal neuropathy"},
                            {"name": "drug withdrawal syndrome"},
                            {"name": "polyuria"},
                            {"name": "generalised tonic-clonic seizure"},
                            {"name": "tachyphylaxis"},
                            {"name": "myoclonus"},
                            {"name": "paradoxical drug reaction"},
                            {"name": "laryngospasm"},
                            {"name": "sedation complication"},
                            {"name": "confusional arousal"},
                            {"name": "quality of life decreased"},
                            {"name": "tachypnoea"},
                            {"name": "hyperthermia malignant"},
                            {"name": "tachycardia"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL277474",
                        "drugbank_id": "DrugBank:DB01435",
                        "compound_cid": "CID:2206",
                        "compound_name": "ANTIPYRINE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 3.0,
                        "adverse_effect": [
                            {"name": "miosis"},
                            {"name": "electrocardiogram q wave abnormal"},
                            {"name": "glucose-6-phosphate dehydrogenase deficiency"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL596",
                        "drugbank_id": "DrugBank:DB00813",
                        "compound_cid": "CID:3345",
                        "compound_name": "FENTANYL",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 46.0,
                        "adverse_effect": [
                            {"name": "serotonin syndrome"},
                            {"name": "respiratory depression"},
                            {"name": "somnolence"},
                            {"name": "drug withdrawal syndrome"},
                            {"name": "inadequate analgesia"},
                            {"name": "miosis"},
                            {"name": "unresponsive to stimuli"},
                            {"name": "drug dependence"},
                            {"name": "confusional state"},
                            {"name": "drug abuse"},
                            {"name": "withdrawal syndrome"},
                            {"name": "oxygen saturation decreased"},
                            {"name": "hypotension"},
                            {"name": "delirium"},
                            {"name": "apnoea"},
                            {"name": "hyperaesthesia"},
                            {"name": "toxicity to various agents"},
                            {"name": "depressed level of consciousness"},
                            {"name": "breakthrough pain"},
                            {"name": "presbyacusis"},
                            {"name": "pharmaceutical product complaint"},
                            {"name": "anaphylactic reaction"},
                            {"name": "drug interaction"},
                            {"name": "coma"},
                            {"name": "bradycardia"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL240163",
                        "drugbank_id": "DrugBank:DB11774",
                        "compound_cid": "CID:65957",
                        "compound_name": "PAZUFLOXACIN",
                        "clincal_trial_phase": 3.0,
                        "is_approved": False,
                        "relation": "treats",
                        "adverse_effect_count": 2.0,
                        "adverse_effect": [
                            {"name": "oculomucocutaneous syndrome"},
                            {"name": "chlamydia serology positive"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201231",
                        "drugbank_id": "DrugBank:DB14631",
                        "compound_cid": "CID:72078",
                        "compound_name": "PREDNISOLONE PHOSPHORIC ACID",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 0.0,
                        "adverse_effect": nan,
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201014",
                        "drugbank_id": None,
                        "compound_cid": "CID:441409",
                        "compound_name": "PREDNISOLONE SODIUM PHOSPHATE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 60.0,
                        "adverse_effect": [
                            {"name": "choroiditis"},
                            {"name": "intraocular pressure increased"},
                            {"name": "endophthalmitis"},
                            {"name": "uveitis"},
                            {"name": "hepatitis b reactivation"},
                            {"name": "visual acuity reduced"},
                            {"name": "chorioretinopathy"},
                            {"name": "eye pain"},
                            {"name": "osteonecrosis"},
                            {"name": "bedridden"},
                            {"name": "condition aggravated"},
                            {"name": "lower limb fracture"},
                            {"name": "osteoporosis"},
                            {"name": "loss of personal independence in daily activities"},
                            {"name": "hypogammaglobulinaemia"},
                            {"name": "hepatitis e"},
                            {"name": "meningomyelitis herpes"},
                            {"name": "foot fracture"},
                            {"name": "epstein-barr virus infection"},
                            {"name": "neutropenia"},
                            {"name": "corneal deposits"},
                            {"name": "cytomegalovirus infection"},
                            {"name": "post transplant lymphoproliferative disorder"},
                            {"name": "diarrhoea"},
                            {"name": "meningitis cryptococcal"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1161",
                        "drugbank_id": "DrugBank:DB14512",
                        "compound_cid": "CID:441336",
                        "compound_name": "MOMETASONE FUROATE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 40.0,
                        "adverse_effect": [
                            {"name": "poor quality device used"},
                            {"name": "dyspnoea"},
                            {"name": "asthma"},
                            {"name": "wheezing"},
                            {"name": "cough"},
                            {"name": "product dose omission issue"},
                            {"name": "gastrooesophageal reflux disease"},
                            {"name": "pneumonia aspiration"},
                            {"name": "wrong technique in device usage process"},
                            {"name": "steroid withdrawal syndrome"},
                            {"name": "poor quality product administered"},
                            {"name": "coma"},
                            {"name": "chest discomfort"},
                            {"name": "nasal congestion"},
                            {"name": "dysphonia"},
                            {"name": "poor quality drug administered"},
                            {"name": "nasal polyps"},
                            {"name": "sleep disorder due to a general medical condition"},
                            {"name": "rhinitis allergic"},
                            {"name": "anosmia"},
                            {"name": "productive cough"},
                            {"name": "diarrhoea"},
                            {"name": "product availability issue"},
                            {"name": "nasal discomfort"},
                            {"name": "product odour abnormal"},
                        ],
                    },
                ],
                [
                    {
                        "chembl_id": nan,
                        "drugbank_id": nan,
                        "compound_cid": nan,
                        "compound_name": nan,
                        "clincal_trial_phase": nan,
                        "is_approved": nan,
                        "relation": nan,
                        "adverse_effect_count": nan,
                        "adverse_effect": nan,
                    }
                ],
                [
                    {
                        "chembl_id": nan,
                        "drugbank_id": nan,
                        "compound_cid": nan,
                        "compound_name": nan,
                        "clincal_trial_phase": nan,
                        "is_approved": nan,
                        "relation": nan,
                        "adverse_effect_count": nan,
                        "adverse_effect": nan,
                    }
                ],
            ]
        )
        expected_data.name = OPENTARGETS_DISEASE_COMPOUND_COL

        pd.testing.assert_series_equal(
            obtained_data[OPENTARGETS_DISEASE_COMPOUND_COL],
            expected_data,
            check_dtype=False,  # to avoid dict checks
        )
        self.assertIsInstance(metadata, dict)
