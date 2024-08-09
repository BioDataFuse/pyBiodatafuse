#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenTargets annotator.

Run as: python -m unittest tests/annotators/test_opentargets.py"""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
import numpy as np

from pyBiodatafuse.annotators import opentargets
from pyBiodatafuse.annotators.opentargets import (
    get_compound_disease_interactions,
    get_gene_compound_interactions,
    get_gene_go_process,
    get_gene_reactome_pathways,
    get_version_opentargets,
)
from pyBiodatafuse.constants import (
    OPENTARGETS_COMPOUND_COL,
    OPENTARGETS_DISEASE_COL,
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
                "source_version": {"apiVersion": {"x": "24", "y": "0", "z": "3"}},
                "data_version": {"dataVersion": {"year": "24", "month": "03"}},
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
                    {"go_id": "GO:0031965", "go_name": "nuclear membrane", "go_type": "C"},
                    {
                        "go_id": "GO:0005789",
                        "go_name": "endoplasmic reticulum membrane",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0043541",
                        "go_name": "UDP-N-acetylglucosamine transferase complex",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0004577",
                        "go_name": "N-acetylglucosaminyldiphosphodolichol N-acetylglucosaminyltransferase activity",
                        "go_type": "F",
                    },
                    {
                        "go_id": "GO:0006488",
                        "go_name": "dolichol-linked oligosaccharide biosynthetic process",
                        "go_type": "P",
                    },
                ],
                [
                    {
                        "go_id": "GO:0046982",
                        "go_name": "protein heterodimerization activity",
                        "go_type": "F",
                    },
                    {"go_id": "GO:0016020", "go_name": "membrane", "go_type": "C"},
                    {
                        "go_id": "GO:0000033",
                        "go_name": "alpha-1,3-mannosyltransferase activity",
                        "go_type": "F",
                    },
                    {"go_id": "GO:0005515", "go_name": "protein binding", "go_type": "F"},
                    {
                        "go_id": "GO:0006490",
                        "go_name": "oligosaccharide-lipid intermediate biosynthetic process",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0005634", "go_name": "nucleus", "go_type": "C"},
                    {
                        "go_id": "GO:0048306",
                        "go_name": "calcium-dependent protein binding",
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
                        "go_id": "GO:0102704",
                        "go_name": "GDP-Man:Man2GlcNAc2-PP-dolichol alpha-1,6-mannosyltransferase activity",
                        "go_type": "F",
                    },
                    {
                        "go_id": "GO:0006488",
                        "go_name": "dolichol-linked oligosaccharide biosynthetic process",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0005737", "go_name": "cytoplasm", "go_type": "C"},
                    {
                        "go_id": "GO:0051592",
                        "go_name": "response to calcium ion",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0048471",
                        "go_name": "perinuclear region of cytoplasm",
                        "go_type": "C",
                    },
                ],
                [
                    {"go_id": "GO:0050905", "go_name": "neuromuscular process", "go_type": "P"},
                    {"go_id": "GO:0009986", "go_name": "cell surface", "go_type": "C"},
                    {
                        "go_id": "GO:0070050",
                        "go_name": "neuron cellular homeostasis",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0043005", "go_name": "neuron projection", "go_type": "C"},
                    {
                        "go_id": "GO:0034220",
                        "go_name": "monoatomic ion transmembrane transport",
                        "go_type": "P",
                    },
                    {"go_id": "GO:0035094", "go_name": "response to nicotine", "go_type": "P"},
                    {
                        "go_id": "GO:0048630",
                        "go_name": "skeletal muscle tissue growth",
                        "go_type": "P",
                    },
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
                    {
                        "go_id": "GO:0031594",
                        "go_name": "neuromuscular junction",
                        "go_type": "C",
                    },
                    {
                        "go_id": "GO:0051899",
                        "go_name": "membrane depolarization",
                        "go_type": "P",
                    },
                    {
                        "go_id": "GO:0050881",
                        "go_name": "musculoskeletal movement",
                        "go_type": "P",
                    },
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
                    {
                        "go_id": "GO:0019228",
                        "go_name": "neuronal action potential",
                        "go_type": "P",
                    },
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
                    {
                        "go_id": "GO:0060079",
                        "go_name": "excitatory postsynaptic potential",
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
                        "go_id": "GO:0007274",
                        "go_name": "neuromuscular synaptic transmission",
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
                        "pathway_id": "R-HSA-446193",
                    },
                    {
                        "pathway_label": "Defective ALG14 causes ALG14-CMS",
                        "pathway_id": "R-HSA-5633231",
                    },
                ],
                [
                    {
                        "pathway_label": "Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein",
                        "pathway_id": "R-HSA-446193",
                    },
                    {
                        "pathway_label": "Defective ALG2 causes CDG-1i",
                        "pathway_id": "R-HSA-4549349",
                    },
                ],
                [
                    {
                        "pathway_label": "Highly calcium permeable nicotinic acetylcholine receptors",
                        "pathway_id": "R-HSA-629597",
                    },
                    {
                        "pathway_label": "Highly calcium permeable postsynaptic nicotinic acetylcholine receptors",
                        "pathway_id": "R-HSA-629594",
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

        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "identifier": ["ALG14", "ALG2", "CHRNA1"],
                "identifier.source": ["HGNC", "HGNC", "HGNC"],
                "target": ["ENSG00000172339", "ENSG00000119523", "ENSG00000138435"],
                "target.source": ["Ensembl", "Ensembl", "Ensembl"],
            }
        )

        with open(os.path.join(data_file_folder, "opentargets_compound_mock_data.json")) as f:
            mock_post_compound.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_gene_compound_interactions(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "chembl_id": np.nan,
                        "drugbank_id": np.nan,
                        "compound_cid": np.nan,
                        "compound_name": np.nan,
                        "is_approved": np.nan,
                        "relation": np.nan,
                        "adverse_effect_count": np.nan,
                        "adverse_effect": np.nan,
                    }
                ],
                [
                    {
                        "chembl_id": np.nan,
                        "drugbank_id": np.nan,
                        "compound_cid": np.nan,
                        "compound_name": np.nan,
                        "is_approved": np.nan,
                        "relation": np.nan,
                        "adverse_effect_count": np.nan,
                        "adverse_effect": np.nan,
                    }
                ],
                [
                    {
                        "chembl_id": "CHEMBL1201248",
                        "drugbank_id": "DB00565",
                        "compound_cid": "62887",
                        "compound_name": "CISATRACURIUM",
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 103.0,
                        "adverse_effect": [
                            {"name": "anaphylactic shock"},
                            {"name": "anaphylactic reaction"},
                            {"name": "feeding intolerance"},
                            {"name": "circulatory collapse"},
                            {"name": "bronchospasm"},
                            {"name": "hyperthermia malignant"},
                            {"name": "hypotension"},
                            {"name": "cardiac arrest"},
                            {"name": "apnoea"},
                            {"name": "cardio-respiratory arrest"},
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
                        "chembl_id": "CHEMBL1200641",
                        "drugbank_id": np.nan,
                        "compound_cid": "62886",
                        "compound_name": "CISATRACURIUM BESYLATE",
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 88.0,
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
                            {"name": "acute pulmonary oedema"},
                            {"name": "post procedural complication"},
                            {"name": "haemodynamic instability"},
                            {"name": "hepatocellular injury"},
                            {"name": "fatigue"},
                            {"name": "hyperbilirubinaemia"},
                            {"name": "bradycardia"},
                            {"name": "diarrhoea"},
                            {"name": "tachyphylaxis"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "headache"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL703",
                        "drugbank_id": "DB00202",
                        "compound_cid": "441290",
                        "compound_name": "SUXAMETHONIUM",
                        "is_approved": True,
                        "relation": "activates",
                        "adverse_effect_count": 127.0,
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
                            {"name": "left ventricular dysfunction"},
                            {"name": "blood ph decreased"},
                            {"name": "foetal death"},
                            {"name": "muscle contractions involuntary"},
                            {"name": "bradycardia"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL1201244",
                        "drugbank_id": "DB00728",
                        "compound_cid": np.nan,
                        "compound_name": "ROCURONIUM",
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 94.0,
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
                            {"name": "pulseless electrical activity"},
                            {"name": "recurrence of neuromuscular blockade"},
                            {"name": "bradycardia"},
                            {"name": "finger deformity"},
                            {"name": "wrist surgery"},
                            {"name": "tenosynovitis"},
                            {"name": "stress cardiomyopathy"},
                            {"name": "oxygen saturation decreased"},
                            {"name": "negative pressure pulmonary oedema"},
                            {"name": "airway peak pressure increased"},
                            {"name": "procedural hypotension"},
                            {"name": "fatigue"},
                            {"name": "hypoventilation"},
                            {"name": "limb operation"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL1201219",
                        "drugbank_id": "DB01339",
                        "compound_cid": np.nan,
                        "compound_name": "VECURONIUM",
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 29.0,
                        "adverse_effect": [
                            {"name": "feeding intolerance"},
                            {"name": "hyperthermia malignant"},
                            {"name": "post procedural complication"},
                            {"name": "carcinoid crisis"},
                            {"name": "bradycardia foetal"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "therapeutic product cross-reactivity"},
                            {"name": "anaphylactic reaction"},
                            {"name": "paralysis"},
                            {"name": "vasoplegia syndrome"},
                            {"name": "hypotension"},
                            {"name": "cardiac arrest"},
                            {"name": "wound infection"},
                            {"name": "atelectasis"},
                            {"name": "neonatal hypoxia"},
                            {"name": "negative pressure pulmonary oedema"},
                            {"name": "urinary retention"},
                            {"name": "non-cardiogenic pulmonary oedema"},
                            {"name": "product packaging confusion"},
                            {"name": "premature baby"},
                            {"name": "haemodynamic instability"},
                            {"name": "nonreassuring foetal heart rate pattern"},
                            {"name": "bradyarrhythmia"},
                            {"name": "drug resistance"},
                            {"name": "intervertebral discitis"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL983",
                        "drugbank_id": np.nan,
                        "compound_cid": np.nan,
                        "compound_name": "SUCCINYLCHOLINE CHLORIDE",
                        "is_approved": True,
                        "relation": "activates",
                        "adverse_effect_count": 130.0,
                        "adverse_effect": [
                            {"name": "hyperthermia malignant"},
                            {"name": "anaphylactic shock"},
                            {"name": "cardiac arrest"},
                            {"name": "neuromuscular block prolonged"},
                            {"name": "hypotension"},
                            {"name": "pseudocholinesterase deficiency"},
                            {"name": "anaphylactic reaction"},
                            {"name": "rhabdomyolysis"},
                            {"name": "anaesthetic complication"},
                            {"name": "ventricular tachycardia"},
                            {"name": "bronchospasm"},
                            {"name": "pco2 increased"},
                            {"name": "pulseless electrical activity"},
                            {"name": "ventricular fibrillation"},
                            {"name": "diabetes insipidus"},
                            {"name": "renal ischaemia"},
                            {"name": "post procedural complication"},
                            {"name": "muscle contractions involuntary"},
                            {"name": "premature baby"},
                            {"name": "left ventricular dysfunction"},
                            {"name": "apnoea"},
                            {"name": "blood ph decreased"},
                            {"name": "foetal death"},
                            {"name": "bradycardia"},
                            {"name": "po2 increased"},
                        ],
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_COMPOUND_COL

        pd.testing.assert_series_equal(obtained_data[OPENTARGETS_COMPOUND_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_compound_disease_interactions(self, mock_post_disease):
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

        bridgedb_dataframe_compounds = pd.DataFrame(
            {
                "identifier": ["100208", "10040286", "10041551", "5291"],
                "identifier.source": [
                    "PubChem-compound",
                    "PubChem-compound",
                    "PubChem-compound",
                    "PubChem-compound",
                ],
                "target": ["100208", "10040286", "10041551", "5291"],
                "target.source": [
                    "PubChem Compound",
                    "PubChem Compound",
                    "PubChem Compound",
                    "PubChem Compound",
                ],
            }
        )

        with open(
            os.path.join(data_file_folder, "opentargets_compound_disease_mock_data.json")
        ) as f:
            mock_post_disease.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_compound_disease_interactions(bridgedb_dataframe_compounds)

        expected_data = pd.Series(
            [
                [
                    {
                        "disease_id": np.nan,
                        "disease_name": np.nan,
                        "therapeutic_areas": np.nan,
                        "disease_xrefs": np.nan,
                    }
                ],
                [
                    {
                        "disease_id": np.nan,
                        "disease_name": np.nan,
                        "therapeutic_areas": np.nan,
                        "disease_xrefs": np.nan,
                    }
                ],
                [
                    {
                        "disease_id": np.nan,
                        "disease_name": np.nan,
                        "therapeutic_areas": np.nan,
                        "disease_xrefs": np.nan,
                    }
                ],
                [
                    {
                        "disease_id": "umls:C0009375",
                        "disease_name": "colonic neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "EFO_0004288",
                            "NCI_C2953",
                            "MESH_D003110",
                            "UMLS_C0009375",
                            "MONDO_0005401",
                        ],
                    },
                    {
                        "disease_id": "umls:C0238033, umls:C0242787, umls:C0242788",
                        "disease_name": "male breast carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "UMLS_C0238033",
                            "DO_1614",
                            "MESH_D018567",
                            "EFO_0006861",
                            "UMLS_C0242787",
                            "MONDO_0005628",
                            "NCI_C3862",
                            "UMLS_C0242788",
                        ],
                    },
                    {
                        "disease_id": "umls:C0007115",
                        "disease_name": "thyroid cancer",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["DO_1781", "UMLS_C0007115", "NCI_C7510"],
                    },
                    {
                        "disease_id": "umls:C0036421",
                        "disease_name": "systemic scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0009690:urinary system disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "NCI_C72070",
                            "MONDO_0005100",
                            "DO_418",
                            "NCI_C72070",
                            "MESH_D012595",
                            "UMLS_C0036421",
                            "MESH_D012595",
                            "ORDO_90291",
                        ],
                    },
                    {
                        "disease_id": "umls:C0751690",
                        "disease_name": "malignant peripheral nerve sheath tumor",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "NCI_C3798",
                            "MONDO_0017827",
                            "MESH_D018317",
                            "DO_5940",
                            "MESH_D018319",
                            "NCI_C3798",
                            "UMLS_C0751690",
                            "ORDO_3148",
                        ],
                    },
                    {
                        "disease_id": "umls:C0036420",
                        "disease_name": "localised scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "NCI_C72069",
                            "DO_8472",
                            "MONDO_0019562",
                            "MESH_D012594",
                            "ORDO_90289",
                            "UMLS_C0036420",
                            "DO_419",
                            "MESH_D012594",
                        ],
                    },
                    {
                        "disease_id": "umls:C1527249",
                        "disease_name": "colorectal cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "ORDO_466667",
                            "DO_9256",
                            "DO_5672",
                            "UMLS_C1527249",
                            "EFO_0005842",
                            "NCI_C4978",
                            "OMIM_114500",
                        ],
                    },
                    {
                        "disease_id": "umls:C0206180",
                        "disease_name": "anaplastic large cell lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "MESH_D017728",
                            "UMLS_C0206180",
                            "NCI_C3720",
                            "MONDO_0020325",
                            "EFO_0003032",
                            "DO_0050744",
                            "ORDO_98841",
                        ],
                    },
                    {
                        "disease_id": "umls:C0032285",
                        "disease_name": "pneumonia",
                        "therapeutic_areas": "EFO_0005741:infectious disease, OTAR_0000010:respiratory or thoracic disease",
                        "disease_xrefs": [
                            "UMLS_C0032285",
                            "MESH_D011014",
                            "MONDO_0005249",
                            "NCI_C3333",
                            "DO_552",
                            "NCI_C3333",
                            "MESH_D011014",
                        ],
                    },
                    {
                        "disease_id": "umls:C2973725",
                        "disease_name": "pulmonary arterial hypertension",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "ORDO_182090",
                            "OMIM_615371",
                            "NCI_C3120",
                            "DO_6432",
                            "UMLS_C2973725",
                            "MESH_D006976",
                            "MONDO_0015924",
                            "MESH_D000081029",
                        ],
                    },
                    {
                        "disease_id": "umls:C3693482",
                        "disease_name": "dermatofibrosarcoma protuberans",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "UMLS_C3693482",
                            "OMIM_607907",
                            "NCI_C4683",
                            "MESH_D018223",
                            "DO_3507",
                            "ORDO_31112",
                        ],
                    },
                    {
                        "disease_id": "umls:CN971653, umls:C0025202",
                        "disease_name": "melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_CN971653",
                            "UMLS_C0025202",
                            "OMIM_155755",
                            "MESH_D008545",
                            "ORDO_411533",
                            "EFO_0000756",
                            "HPO_HP:0002861",
                            "DO_1909",
                            "OMIM_155600",
                            "MONDO_0005105",
                            "NCI_C3224",
                        ],
                    },
                    {
                        "disease_id": "umls:C0004096",
                        "disease_name": "asthma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease",
                        "disease_xrefs": [
                            "UMLS_C0004096",
                            "DO_2841",
                            "EFO_0000270",
                            "HPO_HP:0002099",
                            "MESH_D001249",
                            "NCI_C28397",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023470",
                        "disease_name": "myeloid leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": ["MESH_D007951", "DO_8692", "UMLS_C0023470", "NCI_C3172"],
                    },
                    {
                        "disease_id": "EFO_0007328",
                        "disease_name": "influenza",
                        "therapeutic_areas": "EFO_0005741:infectious disease, OTAR_0000010:respiratory or thoracic disease",
                        "disease_xrefs": [
                            "MESH_D009976",
                            "DO_8469",
                            "NCI_C53482",
                            "MESH_D007251",
                            "MESH_D007251",
                            "MONDO_0005812",
                            "NCI_C53482",
                            "EFO_0007411",
                        ],
                    },
                    {
                        "disease_id": "umls:C1168401",
                        "disease_name": "head and neck squamous cell carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "EFO_0000181",
                            "ORDO_67037",
                            "OMIM_275355",
                            "NCI_C34447",
                            "MESH_D000077195",
                            "MONDO_0010150",
                            "DO_5520",
                            "MESH_C535575",
                            "UMLS_C1168401",
                        ],
                    },
                    {
                        "disease_id": "umls:C0003873",
                        "disease_name": "rheumatoid arthritis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "NCI_C2884",
                            "ORDO_284130",
                            "MONDO_0008383",
                            "MESH_D001172",
                            "NCI_C2884",
                            "UMLS_C0003873",
                            "OMIM_180300",
                            "OMIM_604302",
                            "MESH_D001172",
                            "DO_7148",
                            "HPO_HP:0001370",
                        ],
                    },
                    {
                        "disease_id": "umls:C0008487",
                        "disease_name": "chordoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "HPO_HP:0010762",
                            "MESH_D002817",
                            "OMIM_215400",
                            "NCI_C2947",
                            "ORDO_178",
                            "DO_3302",
                            "UMLS_C0008487",
                        ],
                    },
                    {
                        "disease_id": "umls:C0037286",
                        "disease_name": "skin neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "UMLS_C0037286",
                            "DO_3165",
                            "NCI_C3372",
                            "MESH_D012878",
                            "MONDO_0002531",
                            "DO_4159",
                            "MESH_D012878",
                            "NCI_C3372",
                        ],
                    },
                    {
                        "disease_id": "umls:C0036421",
                        "disease_name": "systemic scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0009690:urinary system disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "NCI_C72070",
                            "MONDO_0005100",
                            "DO_418",
                            "NCI_C72070",
                            "MESH_D012595",
                            "UMLS_C0036421",
                            "MESH_D012595",
                            "ORDO_90291",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000220",
                        "disease_name": "acute lymphoblastic leukemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "MONDO_0004967",
                            "ORDO_513",
                            "EFO_0000220",
                            "DO_1037",
                            "NCI_C3167",
                            "HPO_HP:0006721",
                            "DO_9952",
                        ],
                    },
                    {
                        "disease_id": "umls:C0026769",
                        "disease_name": "multiple sclerosis",
                        "therapeutic_areas": "EFO_0000540:immune system disease, EFO_0000618:nervous system disease",
                        "disease_xrefs": [
                            "DO_2377",
                            "NCI_C3243",
                            "UMLS_C0026769",
                            "MESH_D009103",
                            "ORDO_802",
                            "EFO_0003885",
                        ],
                    },
                    {
                        "disease_id": "umls:C0027831",
                        "disease_name": "neurofibromatosis type 1",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                        "disease_xrefs": [
                            "MESH_D009456",
                            "NCI_C3273",
                            "DO_0111253",
                            "MESH_C538607",
                            "OMIM_162200",
                            "ORDO_636",
                            "UMLS_C0027831",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0005147",
                        "disease_name": "type 1 diabetes mellitus",
                        "therapeutic_areas": "OTAR_0000020:nutritional or metabolic disease, EFO_0001379:endocrine system disease, EFO_0000540:immune system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "DO_9744",
                            "EFO_0001359",
                            "NCI_C2986",
                            "ORDO_243377",
                            "MESH_D003922",
                            "OMIM_222100",
                        ],
                    },
                    {
                        "disease_id": "umls:C0206726",
                        "disease_name": "gliosarcoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "ORDO_251576",
                            "EFO_1001465",
                            "MONDO_0016681",
                            "NCI_C3796",
                            "DO_3071",
                            "UMLS_C0206726",
                            "MESH_D018316",
                        ],
                    },
                    {
                        "disease_id": "umls:C0007115",
                        "disease_name": "thyroid cancer",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["DO_1781", "UMLS_C0007115", "NCI_C7510"],
                    },
                    {
                        "disease_id": "umls:CN244903, umls:C0262584, umls:C0149925",
                        "disease_name": "small cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                        "disease_xrefs": [
                            "MONDO_0008433",
                            "DO_5411",
                            "UMLS_CN244903",
                            "EFO_0000702",
                            "NCI_C4917",
                            "UMLS_C0262584",
                            "UMLS_C0149925",
                            "OMIM_182280",
                            "DO_5409",
                            "MESH_D055752",
                            "ORDO_70573",
                        ],
                    },
                    {
                        "disease_id": "umls:C0032463",
                        "disease_name": "polycythemia vera",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "OMIM_263300",
                            "NCI_C3336",
                            "MESH_D011087",
                            "UMLS_C0032463",
                            "EFO_0002429",
                            "MONDO_0009891",
                            "DO_8997",
                            "ORDO_729",
                        ],
                    },
                    {
                        "disease_id": "umls:C0014457",
                        "disease_name": "Eosinophilia",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                        "disease_xrefs": ["UMLS_C0014457", "MESH_D004802"],
                    },
                    {
                        "disease_id": "umls:C0018133",
                        "disease_name": "graft versus host disease",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": [
                            "DO_0081267",
                            "MESH_D006086",
                            "ORDO_39812",
                            "UMLS_C0018133",
                            "NCI_C3063",
                        ],
                    },
                    {
                        "disease_id": "umls:C0025500",
                        "disease_name": "mesothelioma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_C0025500",
                            "MESH_D008654",
                            "NCI_C3234",
                            "OMIM_156240",
                            "MONDO_0005065",
                            "EFO_0000588",
                        ],
                    },
                    {
                        "disease_id": "umls:C0079218, umls:C1851124, umls:CN072436",
                        "disease_name": "Desmoid-type fibromatosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_C0079218",
                            "DO_0080366",
                            "ORDO_873",
                            "NCI_C9182",
                            "UMLS_C1851124",
                            "OMIM_135290",
                            "MONDO_0007608",
                            "UMLS_CN072436",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023968",
                        "disease_name": "loiasis",
                        "therapeutic_areas": "EFO_0005741:infectious disease, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "UMLS_C0023968",
                            "MESH_D008118",
                            "MONDO_0016566",
                            "DO_13523",
                            "ORDO_2404",
                            "NCI_C34784",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0100096",
                        "disease_name": "COVID-19",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                        "disease_xrefs": ["MESH_C000657245", "MESH_D000086382", "DO_0080600"],
                    },
                    {
                        "disease_id": "umls:C0023473",
                        "disease_name": "chronic myelogenous leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "NCI_C3174",
                            "UMLS_C0023473",
                            "EFO_0000339",
                            "OMIM_608232",
                            "MONDO_0011996",
                            "NCI_C3177",
                            "ORDO_521",
                            "DO_8552",
                            "DO_0081088",
                        ],
                    },
                    {
                        "disease_id": "umls:C0006826",
                        "disease_name": "cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_0050687",
                            "DO_0050686",
                            "UMLS_C0006826",
                            "NCI_C9305",
                            "EFO_0000311",
                            "DO_162",
                        ],
                    },
                    {
                        "disease_id": "umls:C0278883",
                        "disease_name": "metastatic melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "MONDO_0005191",
                            "EFO_0002617",
                            "UMLS_C0278883",
                            "NCI_C8925",
                        ],
                    },
                    {
                        "disease_id": "umls:C0007102",
                        "disease_name": "malignant colon neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": ["NCI_C9242", "UMLS_C0007102", "DO_219"],
                    },
                    {
                        "disease_id": "MONDO_0005149",
                        "disease_name": "pulmonary hypertension",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                        "disease_xrefs": ["DO_6432", "EFO_0001361", "MESH_D006976"],
                    },
                    {
                        "disease_id": "MONDO_0011962",
                        "disease_name": "endometrial cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["EFO_0004230", "OMIM_608089", "DO_1380", "NCI_C27815"],
                    },
                    {
                        "disease_id": "umls:C0034091",
                        "disease_name": "pulmonary venoocclusive disease",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "NCI_C85039",
                            "ORDO_31837",
                            "MESH_D011668",
                            "DO_5453",
                            "UMLS_C0034091",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000691",
                        "disease_name": "sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_1115",
                            "MONDO_0005089",
                            "NCI_C9118",
                            "EFO_0000691",
                            "MESH_D012509",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0008170",
                        "disease_name": "ovarian cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                        "disease_xrefs": [
                            "ORDO_213500",
                            "OMIM_167000",
                            "MESH_D010051",
                            "DO_2394",
                            "NCI_C7431",
                        ],
                    },
                    {
                        "disease_id": "EFO_1000657",
                        "disease_name": "rectum cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": ["EFO_1000657", "DO_1993", "MONDO_0006519", "NCI_C7418"],
                    },
                    {
                        "disease_id": "umls:C0281508",
                        "disease_name": "desmoplastic small round cell tumor",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "MONDO_0019373",
                            "NCI_C8300",
                            "EFO_1000895",
                            "UMLS_C0281508",
                            "ORDO_83469",
                            "MESH_D058405",
                            "DO_6785",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0008903",
                        "disease_name": "lung cancer",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["OMIM_211980", "NCI_C7377", "DO_1324"],
                    },
                    {
                        "disease_id": "MONDO_0007254",
                        "disease_name": "breast cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": ["DO_1612", "NCI_C9335"],
                    },
                    {
                        "disease_id": "umls:C0238463",
                        "disease_name": "papillary thyroid carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                        "disease_xrefs": [
                            "NCI_C4035",
                            "UMLS_C0238463",
                            "DO_3969",
                            "HPO_HP:0002895",
                            "EFO_0000641",
                            "OMIM_188550",
                            "MONDO_0005075",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0100096",
                        "disease_name": "COVID-19",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                        "disease_xrefs": ["MESH_C000657245", "MESH_D000086382", "DO_0080600"],
                    },
                    {
                        "disease_id": "umls:C0235974",
                        "disease_name": "pancreatic carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "UMLS_C0235974",
                            "DO_4905",
                            "OMIM_260350",
                            "NCI_C3850",
                            "MONDO_0005192",
                            "EFO_0002618",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0018364",
                        "disease_name": "malignant epithelial tumor of ovary",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                        "disease_xrefs": ["DO_2151", "MESH_C538090", "ORDO_398934", "NCI_C40026"],
                    },
                    {
                        "disease_id": "MONDO_0002715",
                        "disease_name": "uterine cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["NCI_C3552", "MESH_D014594", "DO_363"],
                    },
                    {
                        "disease_id": "umls:C0007097",
                        "disease_name": "carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "NCI_C2916",
                            "EFO_0000313",
                            "DO_305",
                            "UMLS_C0007097",
                            "MESH_D002277",
                            "MONDO_0004993",
                        ],
                    },
                    {
                        "disease_id": "umls:C0079218, umls:C1851124, umls:CN072436",
                        "disease_name": "Desmoid-type fibromatosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_C0079218",
                            "DO_0080366",
                            "ORDO_873",
                            "NCI_C9182",
                            "UMLS_C1851124",
                            "OMIM_135290",
                            "MONDO_0007608",
                            "UMLS_CN072436",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023493",
                        "disease_name": "T-cell acute lymphoblastic leukemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "MONDO_0004963",
                            "EFO_0000209",
                            "DO_0050523",
                            "DO_5603",
                            "NCI_C3183",
                            "UMLS_C0023493",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0007254",
                        "disease_name": "breast cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": ["DO_1612", "NCI_C9335"],
                    },
                    {
                        "disease_id": "umls:C0280630",
                        "disease_name": "Uterine Carcinosarcoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "MESH_D012192",
                            "NCI_C42700",
                            "DO_6171",
                            "MONDO_0006485",
                            "UMLS_C0280630",
                            "EFO_1000613",
                        ],
                    },
                    {
                        "disease_id": "umls:CN971653, umls:C0025202",
                        "disease_name": "melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_CN971653",
                            "UMLS_C0025202",
                            "OMIM_155755",
                            "MESH_D008545",
                            "ORDO_411533",
                            "EFO_0000756",
                            "HPO_HP:0002861",
                            "DO_1909",
                            "OMIM_155600",
                            "MONDO_0005105",
                            "NCI_C3224",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0001056",
                        "disease_name": "gastric cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": ["DO_10534", "ORDO_63443", "NCI_C9331", "OMIM_613659"],
                    },
                    {
                        "disease_id": "EFO_1001471",
                        "disease_name": "Merkel cell skin cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": ["MESH_D015266", "NCI_C9231"],
                    },
                    {
                        "disease_id": "umls:C0041296",
                        "disease_name": "tuberculosis",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                        "disease_xrefs": [
                            "DO_399",
                            "UMLS_C0041296",
                            "MESH_D014376",
                            "NCI_C3423",
                            "ORDO_3389",
                        ],
                    },
                    {
                        "disease_id": "EFO_1000359",
                        "disease_name": "Malignant Pancreatic Neoplasm",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "NCI_C9005",
                            "OMIM_618680",
                            "MONDO_0009831",
                            "DO_1793",
                            "NCI_C9005",
                        ],
                    },
                    {
                        "disease_id": "umls:C0153535, umls:C0151779, umls:C0153536",
                        "disease_name": "cutaneous melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                        "disease_xrefs": [
                            "UMLS_C0153535",
                            "UMLS_C0151779",
                            "DO_8923",
                            "MONDO_0005012",
                            "OMIM_615848",
                            "EFO_0000389",
                            "UMLS_C0153536",
                            "NCI_C3510",
                            "OMIM_155600",
                            "HPO_HP:0012056",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0008903",
                        "disease_name": "lung cancer",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["OMIM_211980", "NCI_C7377", "DO_1324"],
                    },
                    {
                        "disease_id": "umls:C0008479",
                        "disease_name": "chondrosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_C0008479",
                            "MESH_D002813",
                            "HPO_HP:0006765",
                            "MONDO_0008977",
                            "EFO_0000333",
                            "OMIM_215300",
                            "ORDO_55880",
                            "DO_3371",
                            "NCI_C2946",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000768",
                        "disease_name": "idiopathic pulmonary fibrosis",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "MESH_D054990",
                            "ORDO_2032",
                            "EFO_0000768",
                            "DO_0050156",
                            "OMIM_178500",
                            "NCI_C35716",
                        ],
                    },
                    {
                        "disease_id": "umls:C0007131",
                        "disease_name": "non-small cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "MONDO_0005233",
                            "DO_3908",
                            "NCI_C2926",
                            "EFO_0003060",
                            "UMLS_C0007131",
                            "MESH_D002289",
                            "HPO_HP:0030358",
                            "ORDO_488201",
                            "OMIM_211980",
                        ],
                    },
                    {
                        "disease_id": "umls:C2237660",
                        "disease_name": "wet macular degeneration",
                        "therapeutic_areas": "MONDO_0002025:psychiatric disorder, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000618:nervous system disease",
                        "disease_xrefs": [
                            "MONDO_0005417",
                            "DO_10873",
                            "MESH_D057135",
                            "UMLS_C2237660",
                            "MESH_D057135",
                        ],
                    },
                    {
                        "disease_id": "umls:C2973725",
                        "disease_name": "pulmonary arterial hypertension",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "ORDO_182090",
                            "OMIM_615371",
                            "NCI_C3120",
                            "DO_6432",
                            "UMLS_C2973725",
                            "MESH_D006976",
                            "MONDO_0015924",
                            "MESH_D000081029",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0000870",
                        "disease_name": "childhood acute lymphoblastic leukemia",
                        "therapeutic_areas": "EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["DO_0080144", "NCI_C3168"],
                    },
                    {
                        "disease_id": "umls:C0278996",
                        "disease_name": "head and neck malignant neoplasia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_11934",
                            "NCI_C4013",
                            "EFO_0006859",
                            "UMLS_C0278996",
                            "MONDO_0005627",
                        ],
                    },
                    {
                        "disease_id": "umls:C1621958, umls:C0017636, umls:CN227279",
                        "disease_name": "glioblastoma multiforme",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                        "disease_xrefs": [
                            "UMLS_C1621958",
                            "MESH_D005909",
                            "NCI_C3058",
                            "UMLS_C0017636",
                            "ORDO_360",
                            "HPO_HP:0100843",
                            "UMLS_CN227279",
                            "HPO_HP:0012174",
                            "MONDO_0018177",
                            "DO_3068",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000497",
                        "disease_name": "fibromatosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "EFO_0000497",
                            "MONDO_0005031",
                            "OMIM_135290",
                            "NCI_C3042",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000094",
                        "disease_name": "B-cell acute lymphoblastic leukemia",
                        "therapeutic_areas": "EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "DO_7061",
                            "MONDO_0004947",
                            "EFO_0000094",
                            "DO_0080630",
                            "NCI_C8936",
                        ],
                    },
                    {
                        "disease_id": "umls:C0024305",
                        "disease_name": "non-Hodgkins lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "ORDO_547",
                            "UMLS_C0024305",
                            "DO_0060060",
                            "MONDO_0018908",
                            "NCI_C3211",
                            "MESH_D008228",
                            "EFO_0005952",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023890",
                        "disease_name": "cirrhosis of liver",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "NCI_C2951",
                            "DO_5082",
                            "UMLS_C0023890",
                            "MESH_D008103",
                            "OMIM_215600",
                            "MESH_D008103",
                            "NCI_C2951",
                            "MONDO_0005155",
                        ],
                    },
                    {
                        "disease_id": "umls:C0027651, umls:CN236628",
                        "disease_name": "neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_14566",
                            "NCI_C3262",
                            "UMLS_C0027651",
                            "HPO_HP:0002664",
                            "MESH_D009369",
                            "MONDO_0005070",
                            "UMLS_CN236628",
                            "EFO_0000616",
                        ],
                    },
                    {
                        "disease_id": "umls:C0019034, umls:C0002895",
                        "disease_name": "sickle cell anemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "EFO_1001797",
                            "OMIM_603903",
                            "UMLS_C0019034",
                            "NCI_C34383",
                            "DO_10923",
                            "UMLS_C0002895",
                            "ORDO_232",
                            "MESH_D000755",
                        ],
                    },
                    {
                        "disease_id": "umls:C0238198",
                        "disease_name": "gastrointestinal stromal tumor",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "ORDO_44890",
                            "MESH_D046152",
                            "UMLS_C0238198",
                            "OMIM_606764",
                            "DO_9253",
                            "NCI_C3868",
                        ],
                    },
                    {
                        "disease_id": "umls:C0016057",
                        "disease_name": "fibrosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "HPO_HP:0100244",
                            "MESH_D005354",
                            "MONDO_0005164",
                            "UMLS_C0016057",
                            "NCI_C3043",
                            "EFO_0002087",
                            "ORDO_2030",
                            "DO_3355",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023473",
                        "disease_name": "chronic myelogenous leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "NCI_C3174",
                            "UMLS_C0023473",
                            "EFO_0000339",
                            "OMIM_608232",
                            "MONDO_0011996",
                            "NCI_C3177",
                            "ORDO_521",
                            "DO_8552",
                            "DO_0081088",
                        ],
                    },
                    {
                        "disease_id": "umls:C0024535",
                        "disease_name": "Plasmodium falciparum malaria",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, EFO_0005741:infectious disease, EFO_0005803:hematologic disease, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "UMLS_C0024535",
                            "MONDO_0005920",
                            "MESH_D016778",
                            "MESH_D016778",
                            "DO_14067",
                            "NCI_C34798",
                        ],
                    },
                    {
                        "disease_id": "umls:C0376358",
                        "disease_name": "prostate cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["MESH_D011471", "UMLS_C0376358", "NCI_C7378", "DO_10283"],
                    },
                    {
                        "disease_id": "umls:C0006826",
                        "disease_name": "cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_0050687",
                            "DO_0050686",
                            "UMLS_C0006826",
                            "NCI_C9305",
                            "EFO_0000311",
                            "DO_162",
                        ],
                    },
                    {
                        "disease_id": "umls:C0740457",
                        "disease_name": "kidney cancer",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["MESH_D007680", "NCI_C7548", "UMLS_C0740457", "DO_263"],
                    },
                    {
                        "disease_id": "umls:C0376544",
                        "disease_name": "hematopoietic and lymphoid cell neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": ["UMLS_C0376544", "NCI_C27134"],
                    },
                    {
                        "disease_id": "EFO_1000158",
                        "disease_name": "Central Nervous System Neoplasm",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["NCI_C9293", "MONDO_0006130", "EFO_1000158"],
                    },
                    {
                        "disease_id": "umls:C0018133",
                        "disease_name": "graft versus host disease",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": [
                            "DO_0081267",
                            "MESH_D006086",
                            "ORDO_39812",
                            "UMLS_C0018133",
                            "NCI_C3063",
                        ],
                    },
                    {
                        "disease_id": "umls:C0009402, umls:CN221574",
                        "disease_name": "colorectal carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "DO_0080199",
                            "MONDO_0024331",
                            "UMLS_C0009402",
                            "UMLS_CN221574",
                            "NCI_C2955",
                        ],
                    },
                    {
                        "disease_id": "umls:C0153467",
                        "disease_name": "peritoneum cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": ["DO_1725", "NCI_C3538", "UMLS_C0153467"],
                    },
                    {
                        "disease_id": "umls:C0079748",
                        "disease_name": "lymphoblastic lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": ["NCI_C9360", "UMLS_C0079748", "DO_0080147"],
                    },
                    {
                        "disease_id": "MONDO_0005149",
                        "disease_name": "pulmonary hypertension",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                        "disease_xrefs": ["DO_6432", "EFO_0001361", "MESH_D006976"],
                    },
                    {
                        "disease_id": "umls:C1621958, umls:C0017636, umls:CN227279",
                        "disease_name": "glioblastoma multiforme",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                        "disease_xrefs": [
                            "UMLS_C1621958",
                            "MESH_D005909",
                            "NCI_C3058",
                            "UMLS_C0017636",
                            "ORDO_360",
                            "HPO_HP:0100843",
                            "UMLS_CN227279",
                            "HPO_HP:0012174",
                            "MONDO_0018177",
                            "DO_3068",
                        ],
                    },
                    {
                        "disease_id": "umls:C0238198",
                        "disease_name": "gastrointestinal stromal tumor",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "ORDO_44890",
                            "MESH_D046152",
                            "UMLS_C0238198",
                            "OMIM_606764",
                            "DO_9253",
                            "NCI_C3868",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0011705",
                        "disease_name": "lymphangioleiomyomatosis",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["MESH_D018192", "OMIM_606690", "NCI_C3725"],
                    },
                    {
                        "disease_id": "MONDO_0002158",
                        "disease_name": "fallopian tube cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["DO_1964", "ORDO_180242", "NCI_C7480"],
                    },
                    {
                        "disease_id": "umls:C0278701",
                        "disease_name": "gastric adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "UMLS_C0278701",
                            "DO_3717",
                            "ORDO_464463",
                            "EFO_0000503",
                            "NCI_C4004",
                            "MONDO_0005036",
                        ],
                    },
                    {
                        "disease_id": "umls:C1540912",
                        "disease_name": "Hypereosinophilic syndrome",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, EFO_0000540:immune system disease, EFO_0005803:hematologic disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "NCI_C27038",
                            "ORDO_Orphanet_168956",
                            "ORDO_168956",
                            "MESH_D017681",
                            "MESH_D017681",
                            "MONDO_0015691",
                            "DO_999",
                            "UMLS_C1540912",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000712",
                        "disease_name": "stroke",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "MESH_D020521",
                            "OMIM_601367",
                            "MESH_D020521",
                            "NCI_C3390",
                            "MONDO_0005098",
                            "NCI_C3390",
                            "HPO_HP:0001297",
                        ],
                    },
                    {
                        "disease_id": "EFO_1001814",
                        "disease_name": "nephrogenic fibrosing dermopathy",
                        "therapeutic_areas": "EFO_0010285:integumentary system disease",
                        "disease_xrefs": ["ORDO_Orphanet_137617", "NCI_C84920", "MESH_D054989"],
                    },
                    {
                        "disease_id": "umls:C1540912",
                        "disease_name": "Hypereosinophilic syndrome",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, EFO_0000540:immune system disease, EFO_0005803:hematologic disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "NCI_C27038",
                            "ORDO_Orphanet_168956",
                            "ORDO_168956",
                            "MESH_D017681",
                            "MESH_D017681",
                            "MONDO_0015691",
                            "DO_999",
                            "UMLS_C1540912",
                        ],
                    },
                    {
                        "disease_id": "umls:C0007112",
                        "disease_name": "prostate adenocarcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "NCI_C2919",
                            "OMIM_176807",
                            "DO_2526",
                            "MONDO_0005082",
                            "EFO_0000673",
                            "UMLS_C0007112",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000658",
                        "disease_name": "plexiform neurofibroma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["MESH_D018318", "NCI_C3797", "DO_5151"],
                    },
                    {
                        "disease_id": "umls:C1322286, umls:C0205969, umls:CN207411",
                        "disease_name": "Thymic Carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "EFO_1000576",
                            "DO_3284",
                            "MONDO_0006451",
                            "UMLS_C1322286",
                            "UMLS_C0205969",
                            "ORDO_99868",
                            "NCI_C7569",
                            "DO_4554",
                            "UMLS_CN207411",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000691",
                        "disease_name": "sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_1115",
                            "MONDO_0005089",
                            "NCI_C9118",
                            "EFO_0000691",
                            "MESH_D012509",
                        ],
                    },
                    {
                        "disease_id": "umls:C0011644",
                        "disease_name": "scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "NCI_C0011644",
                            "NCI_C26746",
                            "MESH_D012594",
                            "MONDO_0019340",
                            "DO_419",
                            "ORDO_801",
                            "UMLS_C0011644",
                            "HPO_HP:0100324",
                        ],
                    },
                    {
                        "disease_id": "umls:C0948968, umls:C0001815, umls:C2355576",
                        "disease_name": "primary myelofibrosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "MONDO_0009692",
                            "UMLS_C0948968",
                            "DO_4971",
                            "MESH_D055728",
                            "OMIM_254450",
                            "ORDO_824",
                            "UMLS_C0001815",
                            "NCI_C2862",
                            "UMLS_C2355576",
                            "EFO_0002430",
                        ],
                    },
                    {
                        "disease_id": "umls:C0699790, umls:C1319315",
                        "disease_name": "colorectal adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "UMLS_C0699790",
                            "MONDO_0005008",
                            "DO_0050913",
                            "NCI_C5105",
                            "EFO_0000365",
                            "DO_0050861",
                            "UMLS_C1319315",
                            "OMIM_114500",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023467",
                        "disease_name": "acute myeloid leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "ORDO_519",
                            "UMLS_C0023467",
                            "OMIM_601626",
                            "MESH_D015470",
                            "DO_9119",
                            "EFO_0000222",
                            "MONDO_0018874",
                            "NCI_C3171",
                        ],
                    },
                    {
                        "disease_id": "EFO_0004610",
                        "disease_name": "acute lung injury",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": [
                            "MONDO_0015796",
                            "ORDO_178320",
                            "MESH_D055371",
                            "NCI_C155766",
                            "MESH_D055371",
                        ],
                    },
                    {
                        "disease_id": "umls:C0039101",
                        "disease_name": "synovial sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "UMLS_C0039101",
                            "DO_5485",
                            "EFO_0001376",
                            "ORDO_3273",
                            "MONDO_0010434",
                            "OMIM_300813",
                            "NCI_C3400",
                            "MESH_D013584",
                            "HPO_HP:0012570",
                        ],
                    },
                    {
                        "disease_id": "EFO_0003956",
                        "disease_name": "seasonal allergic rhinitis",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000010:respiratory or thoracic disease",
                        "disease_xrefs": [
                            "NCI_C92188",
                            "MONDO_0005324",
                            "MESH_D006255",
                            "NCI_C92188",
                            "MESH_D006255",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0044917",
                        "disease_name": "T-lymphoblastic lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": ["NCI_C6919", "EFO_1001830"],
                    },
                    {
                        "disease_id": "umls:C1292778",
                        "disease_name": "chronic myeloproliferative disorder",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "OMIM_131440",
                            "NCI_C103126",
                            "MONDO_0020076",
                            "EFO_0004251",
                            "ORDO_98274",
                            "UMLS_C1292778",
                            "NCI_C4345",
                            "DO_2226",
                        ],
                    },
                    {
                        "disease_id": "umls:C0027651, umls:CN236628",
                        "disease_name": "neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_14566",
                            "NCI_C3262",
                            "UMLS_C0027651",
                            "HPO_HP:0002664",
                            "MESH_D009369",
                            "MONDO_0005070",
                            "UMLS_CN236628",
                            "EFO_0000616",
                        ],
                    },
                    {
                        "disease_id": "umls:C0867389",
                        "disease_name": "chronic graft versus host disease",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": ["NCI_C4981", "ORDO_99921", "UMLS_C0867389"],
                    },
                    {
                        "disease_id": "umls:C0023418",
                        "disease_name": "leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "NCI_C3161",
                            "UMLS_C0023418",
                            "DO_1240",
                            "MESH_D007938",
                            "EFO_0000565",
                            "HPO_HP:0001909",
                            "MONDO_0005059",
                        ],
                    },
                    {
                        "disease_id": "umls:C0024305",
                        "disease_name": "non-Hodgkins lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "ORDO_547",
                            "UMLS_C0024305",
                            "DO_0060060",
                            "MONDO_0018908",
                            "NCI_C3211",
                            "MESH_D008228",
                            "EFO_0005952",
                        ],
                    },
                    {
                        "disease_id": "EFO_1001919",
                        "disease_name": "Spinal cord injury",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": ["MESH_D013119"],
                    },
                    {
                        "disease_id": "umls:C0152271, umls:C0023448",
                        "disease_name": "lymphoid leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "MESH_D007945",
                            "HPO_HP:0005526",
                            "EFO_0004289",
                            "UMLS_C0152271",
                            "MONDO_0005402",
                            "UMLS_C0023448",
                            "DO_10747",
                            "DO_1037",
                            "NCI_C7539",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000220",
                        "disease_name": "acute lymphoblastic leukemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "MONDO_0004967",
                            "ORDO_513",
                            "EFO_0000220",
                            "DO_1037",
                            "NCI_C3167",
                            "HPO_HP:0006721",
                            "DO_9952",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000349",
                        "disease_name": "clear cell renal carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": [
                            "DO_4467",
                            "NCI_C4033",
                            "OMIM_144700",
                            "ORDO_319276",
                            "EFO_0000349",
                            "MONDO_0005005",
                        ],
                    },
                    {
                        "disease_id": "umls:C0079748",
                        "disease_name": "lymphoblastic lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": ["NCI_C9360", "UMLS_C0079748", "DO_0080147"],
                    },
                    {
                        "disease_id": "umls:C0023418",
                        "disease_name": "leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "NCI_C3161",
                            "UMLS_C0023418",
                            "DO_1240",
                            "MESH_D007938",
                            "EFO_0000565",
                            "HPO_HP:0001909",
                            "MONDO_0005059",
                        ],
                    },
                    {
                        "disease_id": "umls:C0024299",
                        "disease_name": "lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                        "disease_xrefs": [
                            "EFO_0000574",
                            "UMLS_C0024299",
                            "MESH_D008223",
                            "ORDO_223735",
                            "OMIM_605027",
                            "DO_0060058",
                            "MONDO_0005062",
                            "NCI_C3208",
                        ],
                    },
                    {
                        "disease_id": "umls:C0038356",
                        "disease_name": "stomach neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "NCI_C3387",
                            "MESH_D013274",
                            "UMLS_C0038356",
                            "MESH_D013274",
                            "NCI_C3387",
                            "MONDO_0021085",
                            "DO_10534",
                        ],
                    },
                    {
                        "disease_id": "umls:C0023467",
                        "disease_name": "acute myeloid leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                        "disease_xrefs": [
                            "ORDO_519",
                            "UMLS_C0023467",
                            "OMIM_601626",
                            "MESH_D015470",
                            "DO_9119",
                            "EFO_0000222",
                            "MONDO_0018874",
                            "NCI_C3171",
                        ],
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_DISEASE_COL

        pd.testing.assert_series_equal(obtained_data[OPENTARGETS_DISEASE_COL], expected_data)
        self.assertIsInstance(metadata, dict)

    @patch("pyBiodatafuse.annotators.opentargets.requests.post")
    def test_get_gene_disease_interactions(self, mock_post_disease):
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

        # Mocking functions
        opentargets.get_gene_compound_interactions = Mock(
            return_value=(
                pd.DataFrame(
                    {
                        "identifier": ["ALG14", "CHRNA1"],
                        "identifier.source": ["HGNC", "HGNC"],
                        "target": ["ENSG00000172339", "ENSG00000138435"],
                        "target.source": ["Ensembl", "Ensembl"],
                        OPENTARGETS_COMPOUND_COL: [
                            [
                                {
                                    "chembl_id": np.nan,
                                    "drugbank_id": np.nan,
                                    "compound_cid": np.nan,
                                    "compound_name": np.nan,
                                    "is_approved": np.nan,
                                    "relation": np.nan,
                                    "adverse_effect_count": np.nan,
                                    "adverse_effect": np.nan,
                                }
                            ],
                            [
                                {
                                    "chembl_id": "CHEMBL1201248",
                                    "drugbank_id": "DB00565",
                                    "compound_cid": "62887",
                                    "compound_name": "CISATRACURIUM",
                                    "is_approved": True,
                                    "relation": "inhibits",
                                    "adverse_effect_count": 103.0,
                                    "adverse_effect": [
                                        {"name": "anaphylactic shock"},
                                        {"name": "anaphylactic reaction"},
                                        {"name": "feeding intolerance"},
                                        {"name": "circulatory collapse"},
                                        {"name": "bronchospasm"},
                                        {"name": "hyperthermia malignant"},
                                        {"name": "hypotension"},
                                        {"name": "cardiac arrest"},
                                        {"name": "apnoea"},
                                        {"name": "cardio-respiratory arrest"},
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
                                        {
                                            "name": "drug reaction with eosinophilia and systemic symptoms"
                                        },
                                        {"name": "diarrhoea"},
                                        {"name": "tachyphylaxis"},
                                        {"name": "myocardial stunning"},
                                        {"name": "rash maculo-papular"},
                                    ],
                                },
                                {
                                    "chembl_id": "CHEMBL1200641",
                                    "drugbank_id": None,
                                    "compound_cid": "62886",
                                    "compound_name": "CISATRACURIUM BESYLATE",
                                    "is_approved": True,
                                    "relation": "inhibits",
                                    "adverse_effect_count": 88.0,
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
                                        {
                                            "name": "drug reaction with eosinophilia and systemic symptoms"
                                        },
                                        {"name": "shock"},
                                        {"name": "apnoea"},
                                        {"name": "toxic epidermal necrolysis"},
                                        {"name": "acute pulmonary oedema"},
                                        {"name": "post procedural complication"},
                                        {"name": "haemodynamic instability"},
                                        {"name": "hepatocellular injury"},
                                        {"name": "fatigue"},
                                        {"name": "hyperbilirubinaemia"},
                                        {"name": "bradycardia"},
                                        {"name": "diarrhoea"},
                                        {"name": "tachyphylaxis"},
                                        {"name": "neuromuscular block prolonged"},
                                        {"name": "headache"},
                                    ],
                                },
                                {
                                    "chembl_id": "CHEMBL703",
                                    "drugbank_id": "DB00202",
                                    "compound_cid": "441290",
                                    "compound_name": "SUXAMETHONIUM",
                                    "is_approved": True,
                                    "relation": "activates",
                                    "adverse_effect_count": 127.0,
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
                                        {"name": "left ventricular dysfunction"},
                                        {"name": "blood ph decreased"},
                                        {"name": "foetal death"},
                                        {"name": "muscle contractions involuntary"},
                                        {"name": "bradycardia"},
                                    ],
                                },
                                {
                                    "chembl_id": "CHEMBL1201244",
                                    "drugbank_id": "DB00728",
                                    "compound_cid": np.nan,
                                    "compound_name": "ROCURONIUM",
                                    "is_approved": True,
                                    "relation": "inhibits",
                                    "adverse_effect_count": 94.0,
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
                                        {"name": "pulseless electrical activity"},
                                        {"name": "recurrence of neuromuscular blockade"},
                                        {"name": "bradycardia"},
                                        {"name": "finger deformity"},
                                        {"name": "wrist surgery"},
                                        {"name": "tenosynovitis"},
                                        {"name": "stress cardiomyopathy"},
                                        {"name": "oxygen saturation decreased"},
                                        {"name": "negative pressure pulmonary oedema"},
                                        {"name": "airway peak pressure increased"},
                                        {"name": "procedural hypotension"},
                                        {"name": "fatigue"},
                                        {"name": "hypoventilation"},
                                        {"name": "limb operation"},
                                    ],
                                },
                                {
                                    "chembl_id": "CHEMBL1201219",
                                    "drugbank_id": "DB01339",
                                    "compound_cid": np.nan,
                                    "compound_name": "VECURONIUM",
                                    "is_approved": True,
                                    "relation": "inhibits",
                                    "adverse_effect_count": 29.0,
                                    "adverse_effect": [
                                        {"name": "feeding intolerance"},
                                        {"name": "hyperthermia malignant"},
                                        {"name": "post procedural complication"},
                                        {"name": "carcinoid crisis"},
                                        {"name": "bradycardia foetal"},
                                        {"name": "neuromuscular block prolonged"},
                                        {"name": "therapeutic product cross-reactivity"},
                                        {"name": "anaphylactic reaction"},
                                        {"name": "paralysis"},
                                        {"name": "vasoplegia syndrome"},
                                        {"name": "hypotension"},
                                        {"name": "cardiac arrest"},
                                        {"name": "wound infection"},
                                        {"name": "atelectasis"},
                                        {"name": "neonatal hypoxia"},
                                        {"name": "negative pressure pulmonary oedema"},
                                        {"name": "urinary retention"},
                                        {"name": "non-cardiogenic pulmonary oedema"},
                                        {"name": "product packaging confusion"},
                                        {"name": "premature baby"},
                                        {"name": "haemodynamic instability"},
                                        {"name": "nonreassuring foetal heart rate pattern"},
                                        {"name": "bradyarrhythmia"},
                                        {"name": "drug resistance"},
                                        {"name": "intervertebral discitis"},
                                    ],
                                },
                                {
                                    "chembl_id": "CHEMBL983",
                                    "drugbank_id": None,
                                    "compound_cid": np.nan,
                                    "compound_name": "SUCCINYLCHOLINE CHLORIDE",
                                    "is_approved": True,
                                    "relation": "activates",
                                    "adverse_effect_count": 130.0,
                                    "adverse_effect": [
                                        {"name": "hyperthermia malignant"},
                                        {"name": "anaphylactic shock"},
                                        {"name": "cardiac arrest"},
                                        {"name": "neuromuscular block prolonged"},
                                        {"name": "hypotension"},
                                        {"name": "pseudocholinesterase deficiency"},
                                        {"name": "anaphylactic reaction"},
                                        {"name": "rhabdomyolysis"},
                                        {"name": "anaesthetic complication"},
                                        {"name": "ventricular tachycardia"},
                                        {"name": "bronchospasm"},
                                        {"name": "pco2 increased"},
                                        {"name": "pulseless electrical activity"},
                                        {"name": "ventricular fibrillation"},
                                        {"name": "diabetes insipidus"},
                                        {"name": "renal ischaemia"},
                                        {"name": "post procedural complication"},
                                        {"name": "muscle contractions involuntary"},
                                        {"name": "premature baby"},
                                        {"name": "left ventricular dysfunction"},
                                        {"name": "apnoea"},
                                        {"name": "blood ph decreased"},
                                        {"name": "foetal death"},
                                        {"name": "bradycardia"},
                                        {"name": "po2 increased"},
                                    ],
                                },
                            ],
                        ],
                    }
                ),
                {
                    "datasource": "Open Targets GraphQL & REST API Beta",
                    "metadata": {
                        "source_version": {"apiVersion": {"x": "24", "y": "1", "z": "4"}},
                        "data_version": {"dataVersion": {"year": "24", "month": "06"}},
                    },
                    "query": {
                        "size": 2,
                        "input_type": "Ensembl",
                        "number_of_added_nodes": 6,
                        "number_of_added_edges": 6,
                        "time": "0:00:00.541977",
                        "date": "2024-08-09 16:44:14",
                        "url": "https://api.platform.opentargets.org/api/v4/graphql",
                    },
                },
            )
        )
        opentargets.get_gene_disease_interactions = Mock(
            return_value=(
                pd.read_csv(os.path.join(data_file_folder, "opentargets_gene_disease_test.csv")),
                {
                    "datasource": "Open Targets GraphQL & REST API Beta",
                    "metadata": {
                        "source_version": {"apiVersion": {"x": "24", "y": "1", "z": "4"}},
                        "data_version": {"dataVersion": {"year": "24", "month": "06"}},
                    },
                    "query": {
                        "size": 6,
                        "number_of_added_nodes": 14,
                        "number_of_added_edges": 17,
                        "time": "0:00:00.162167",
                        "date": "2024-08-09 16:49:26",
                        "url": "https://api.platform.opentargets.org/api/v4/graphql",
                        "input_type": "chembl_id",
                    },
                },
            )
        )

        with open(
            os.path.join(data_file_folder, "opentargets_compound_disease_mock_data.json")
        ) as f:
            mock_post_disease.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_compound_disease_interactions(bridgedb_dataframe_genes)

        expected_data = pd.Series(
            [
                [
                    {
                        "disease_id": np.nan,
                        "disease_name": np.nan,
                        "therapeutic_areas": np.nan,
                        "disease_xrefs": np.nan,
                    }
                ],
                [
                    {
                        "disease_id": "EFO_1000637",
                        "disease_name": "acute respiratory distress syndrome",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease",
                        "disease_xrefs": [
                            "MESH_D012128",
                            "NCI_C3353",
                            "NCI_C3353",
                            "MONDO_0006502",
                        ],
                    },
                    {
                        "disease_id": "EFO_0000712",
                        "disease_name": "stroke",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "MESH_D020521",
                            "OMIM_601367",
                            "MESH_D020521",
                            "NCI_C3390",
                            "MONDO_0005098",
                            "NCI_C3390",
                            "HPO_HP:0001297",
                        ],
                    },
                    {
                        "disease_id": "umls:C0428977",
                        "disease_name": "Bradycardia",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                        "disease_xrefs": ["UMLS_C0428977", "MESH_D001919"],
                    },
                    {
                        "disease_id": "umls:C0231528, umls:C1963177",
                        "disease_name": "Myalgia",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                        "disease_xrefs": [
                            "NCI_C27009",
                            "UMLS_C0231528",
                            "MESH_D063806",
                            "UMLS_C1963177",
                        ],
                    },
                    {
                        "disease_id": "umls:C0158288",
                        "disease_name": "spinal stenosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease",
                        "disease_xrefs": [
                            "NCI_C177444",
                            "HPO_HP:0003416",
                            "MESH_D013130",
                            "UMLS_C0158288",
                            "DO_6725",
                            "MONDO_0005965",
                            "MESH_D013130",
                        ],
                    },
                    {
                        "disease_id": "umls:C0018790, umls:C0444720",
                        "disease_name": "cardiac arrest",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "UMLS_C0018790",
                            "MONDO_0000745",
                            "MESH_D006323",
                            "UMLS_C0444720",
                            "HPO_HP:0001695",
                        ],
                    },
                    {
                        "disease_id": "EFO_0003843",
                        "disease_name": "pain",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                        "disease_xrefs": ["NCI_C3303", "MESH_D010146"],
                    },
                    {
                        "disease_id": "EFO_0010581",
                        "disease_name": "organophosphate poisoning",
                        "therapeutic_areas": "OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": ["MESH_68062025"],
                    },
                    {
                        "disease_id": "EFO_0003931",
                        "disease_name": "bone fracture",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000009:injury, poisoning or other complication",
                        "disease_xrefs": [
                            "MONDO_0005315",
                            "MESH_D050723",
                            "MESH_D050723",
                            "NCI_C3046",
                        ],
                    },
                    {
                        "disease_id": "umls:C0018790, umls:C0444720",
                        "disease_name": "cardiac arrest",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                        "disease_xrefs": [
                            "UMLS_C0018790",
                            "MONDO_0000745",
                            "MESH_D006323",
                            "UMLS_C0444720",
                            "HPO_HP:0001695",
                        ],
                    },
                    {
                        "disease_id": "umls:C0231528, umls:C1963177",
                        "disease_name": "Myalgia",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                        "disease_xrefs": [
                            "NCI_C27009",
                            "UMLS_C0231528",
                            "MESH_D063806",
                            "UMLS_C1963177",
                        ],
                    },
                    {
                        "disease_id": "umls:C0022658",
                        "disease_name": "kidney disease",
                        "therapeutic_areas": "EFO_0009690:urinary system disease",
                        "disease_xrefs": [
                            "MONDO_0005240",
                            "NCI_C3149",
                            "NCI_C34843",
                            "MESH_D007674",
                            "NCI_C3149",
                            "UMLS_C0022658",
                            "MESH_D007674",
                            "DO_557",
                        ],
                    },
                    {
                        "disease_id": "umls:C0004604",
                        "disease_name": "Back pain",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                        "disease_xrefs": ["UMLS_C0004604", "MESH_D001416"],
                    },
                    {
                        "disease_id": "umls:C0376358",
                        "disease_name": "prostate cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["MESH_D011471", "UMLS_C0376358", "NCI_C7378", "DO_10283"],
                    },
                    {
                        "disease_id": "umls:C0235974",
                        "disease_name": "pancreatic carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                        "disease_xrefs": [
                            "UMLS_C0235974",
                            "DO_4905",
                            "OMIM_260350",
                            "NCI_C3850",
                            "MONDO_0005192",
                            "EFO_0002618",
                        ],
                    },
                    {
                        "disease_id": "MONDO_0021117",
                        "disease_name": "lung neoplasm",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                        "disease_xrefs": ["NCI_C3200", "MESH_D008175"],
                    },
                    {
                        "disease_id": "EFO_1000637",
                        "disease_name": "acute respiratory distress syndrome",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease",
                        "disease_xrefs": [
                            "MESH_D012128",
                            "NCI_C3353",
                            "NCI_C3353",
                            "MONDO_0006502",
                        ],
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_DISEASE_COL

        pd.testing.assert_series_equal(obtained_data[OPENTARGETS_DISEASE_COL], expected_data)
        self.assertIsInstance(metadata, dict)
