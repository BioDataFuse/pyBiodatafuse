#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the OpenTargets annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd

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
                        "chembl_id": None,
                        "drugbank_id": None,
                        "compound_cid": None,
                        "compound_name": None,
                        "is_approved": None,
                        "relation": None,
                        "adverse_effect_count": None,
                        "adverse_effect": None,
                    }
                ],
                [
                    {
                        "chembl_id": None,
                        "drugbank_id": None,
                        "compound_cid": None,
                        "compound_name": None,
                        "is_approved": None,
                        "relation": None,
                        "adverse_effect_count": None,
                        "adverse_effect": None,
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
                        "compound_cid": None,
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
                        "compound_cid": None,
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
                        "compound_cid": None,
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

        bridgedb_dataframe_genes = pd.DataFrame(
            {
                "target": ["CHEMBL1201583", "CHEMBL941"],
                "target.source": ["chembl", "chembl"],
            }
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
                        "disease_id": "EFO_0000326",
                        "disease_name": "central nervous system cancer",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1335302",
                        "disease_name": "pancreatic ductal adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0345958",
                        "disease_name": "large cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0280630",
                        "disease_name": "Uterine Carcinosarcoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0684337, umls:C3489398, umls:C0877849, umls:C0553580",
                        "disease_name": "Ewing sarcoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0206686",
                        "disease_name": "adrenal cortex carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, EFO_0000319:cardiovascular disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "umls:C0278996",
                        "disease_name": "head and neck malignant neoplasia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0037930",
                        "disease_name": "spinal cord neoplasm",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0002158",
                        "disease_name": "fallopian tube cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0008903",
                        "disease_name": "lung cancer",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0027859",
                        "disease_name": "Vestibular schwannoma",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "umls:C0006826",
                        "disease_name": "cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1266044",
                        "disease_name": "collecting duct carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1515289, umls:C0206724",
                        "disease_name": "sex cord-stromal tumor",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0030354",
                        "disease_name": "papilloma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0007224",
                        "disease_name": "coronavirus infectious disease",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                    },
                    {
                        "disease_id": "MONDO_0018364",
                        "disease_name": "malignant epithelial tumor of ovary",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "umls:C0302592",
                        "disease_name": "cervical carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0008497",
                        "disease_name": "choriocarcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0555198",
                        "disease_name": "malignant glioma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0235974",
                        "disease_name": "pancreatic carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0011882",
                        "disease_name": "diabetic neuropathy",
                        "therapeutic_areas": "EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "EFO_1001100",
                        "disease_name": "peritoneal neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C1458155",
                        "disease_name": "breast neoplasm",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "EFO_0009708",
                        "disease_name": "metastasis",
                        "therapeutic_areas": "GO_0008150:biological_process",
                    },
                    {
                        "disease_id": "umls:C0032285",
                        "disease_name": "pneumonia",
                        "therapeutic_areas": "EFO_0005741:infectious disease, OTAR_0000010:respiratory or thoracic disease",
                    },
                    {
                        "disease_id": "umls:C1561643",
                        "disease_name": "chronic kidney disease",
                        "therapeutic_areas": "EFO_0009690:urinary system disease",
                    },
                    {
                        "disease_id": "EFO_0009637",
                        "disease_name": "pleural effusion",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "EFO_0000182",
                        "disease_name": "hepatocellular carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "EFO_0004243",
                        "disease_name": "carcinoid tumor",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C2145472",
                        "disease_name": "urothelial carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0750952",
                        "disease_name": "biliary tract cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0037661",
                        "disease_name": "somatostatinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0007131",
                        "disease_name": "non-small cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0005684",
                        "disease_name": "urinary bladder cancer",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:CN971653, umls:C0025202",
                        "disease_name": "melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0206696",
                        "disease_name": "signet ring cell carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_1000158",
                        "disease_name": "Central Nervous System Neoplasm",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000349",
                        "disease_name": "clear cell renal carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0008526",
                        "disease_name": "choroiditis",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0238198",
                        "disease_name": "gastrointestinal stromal tumor",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0020634",
                        "disease_name": "grade III meningioma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0271051",
                        "disease_name": "macular retinal edema",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0812413",
                        "disease_name": "malignant pleural mesothelioma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0024437",
                        "disease_name": "macular degeneration",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "EFO_1000657",
                        "disease_name": "rectum cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0280793",
                        "disease_name": "oligoastrocytoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0042164",
                        "disease_name": "uveitis",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0152013",
                        "disease_name": "lung adenocarcinoma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C3539878, umls:C4722518",
                        "disease_name": "triple-negative breast cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0033999",
                        "disease_name": "pterygium",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system, MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0003962",
                        "disease_name": "Ascites",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "umls:C0206726",
                        "disease_name": "gliosarcoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023434, umls:C0855095",
                        "disease_name": "chronic lymphocytic leukemia",
                        "therapeutic_areas": "EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0281508",
                        "disease_name": "desmoplastic small round cell tumor",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0024677",
                        "disease_name": "pancreatic insulinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0278601",
                        "disease_name": "inflammatory breast carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "EFO_0005577",
                        "disease_name": "pharynx cancer",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0678222",
                        "disease_name": "breast carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0016057",
                        "disease_name": "fibrosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000632",
                        "disease_name": "oligodendroglioma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "EFO_0001416",
                        "disease_name": "cervical adenocarcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0549473",
                        "disease_name": "thyroid carcinoma",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0039446",
                        "disease_name": "telangiectasis",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C3463824, umls:C0033027",
                        "disease_name": "myelodysplastic syndrome",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0040100",
                        "disease_name": "Thymoma",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, EFO_0000540:immune system disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C1175175",
                        "disease_name": "severe acute respiratory syndrome",
                        "therapeutic_areas": "EFO_0005741:infectious disease, OTAR_0000010:respiratory or thoracic disease",
                    },
                    {
                        "disease_id": "umls:C0007112",
                        "disease_name": "prostate adenocarcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023269",
                        "disease_name": "leiomyosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005741:infectious disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0024301",
                        "disease_name": "follicular lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0034091",
                        "disease_name": "pulmonary venoocclusive disease",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0025500",
                        "disease_name": "mesothelioma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0011962",
                        "disease_name": "endometrial cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0009402, umls:CN221574",
                        "disease_name": "colorectal carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0044903",
                        "disease_name": "myelofibrosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0007137",
                        "disease_name": "squamous cell carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0376358",
                        "disease_name": "prostate cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0699790, umls:C1319315",
                        "disease_name": "colorectal adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0338106",
                        "disease_name": "colon adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0034069",
                        "disease_name": "pulmonary fibrosis",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0153467",
                        "disease_name": "peritoneum cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0259783",
                        "disease_name": "mixed glioma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023903, umls:C0345904",
                        "disease_name": "liver cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "EFO_0003833",
                        "disease_name": "brain neoplasm",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0017638",
                        "disease_name": "glioma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0278701",
                        "disease_name": "gastric adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0334579",
                        "disease_name": "anaplastic astrocytoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "EFO_0004212",
                        "disease_name": "Keloid",
                        "therapeutic_areas": "EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "EFO_0004284",
                        "disease_name": "upper aerodigestive tract neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0278883",
                        "disease_name": "metastatic melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0206666",
                        "disease_name": "placental site trophoblastic tumor",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, OTAR_0000014:pregnancy or perinatal disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0005753",
                        "disease_name": "ocular vascular disease",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0042373",
                        "disease_name": "vascular disease",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "MONDO_0044937",
                        "disease_name": "rectal carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MP_0001914",
                        "disease_name": "hemorrhage",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "umls:C1266042, umls:C3887514",
                        "disease_name": "chromophobe renal cell carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1527249",
                        "disease_name": "colorectal cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:CN205129, umls:C1306837, umls:C1336078",
                        "disease_name": "papillary renal cell carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0206644",
                        "disease_name": "benign fibrous histiocytoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0007115",
                        "disease_name": "thyroid cancer",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0007532",
                        "disease_name": "uterine corpus cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0345905",
                        "disease_name": "intrahepatic cholangiocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0392788",
                        "disease_name": "extranodal nasal NK/T cell lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, EFO_0005741:infectious disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "MONDO_0007254",
                        "disease_name": "breast cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C1621958, umls:C0017636, umls:CN227279",
                        "disease_name": "glioblastoma multiforme",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0024299",
                        "disease_name": "lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0039101",
                        "disease_name": "synovial sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0206655",
                        "disease_name": "alveolar rhabdomyosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0242383",
                        "disease_name": "age-related macular degeneration",
                        "therapeutic_areas": "MONDO_0002025:psychiatric disorder, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0279663",
                        "disease_name": "ovarian serous cystadenocarcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "umls:C0027819",
                        "disease_name": "neuroblastoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0598798",
                        "disease_name": "lymphoid neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_0001065",
                        "disease_name": "endometriosis",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease",
                    },
                    {
                        "disease_id": "umls:C0206695",
                        "disease_name": "neuroendocrine carcinoma",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C2237660",
                        "disease_name": "wet macular degeneration",
                        "therapeutic_areas": "MONDO_0002025:psychiatric disorder, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "EFO_0000183",
                        "disease_name": "Hodgkins lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, EFO_0005741:infectious disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0206681",
                        "disease_name": "clear cell adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0152018",
                        "disease_name": "esophageal carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0035412",
                        "disease_name": "rhabdomyosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0205698",
                        "disease_name": "undifferentiated carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023473",
                        "disease_name": "chronic myelogenous leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0023827",
                        "disease_name": "liposarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1333407, umls:C0280788",
                        "disease_name": "anaplastic ependymoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0007097",
                        "disease_name": "carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_1001471",
                        "disease_name": "Merkel cell skin cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0035305",
                        "disease_name": "retinal detachment",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0025149, umls:C1334410",
                        "disease_name": "medulloblastoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0024305",
                        "disease_name": "non-Hodgkins lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0751711",
                        "disease_name": "anterior ischemic optic neuropathy",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C1334633",
                        "disease_name": "neoplasm of mature B-cells",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0206657, umls:C0279544",
                        "disease_name": "alveolar soft part sarcoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0206698",
                        "disease_name": "cholangiocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MP_0001845",
                        "disease_name": "inflammation",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "umls:C1328479",
                        "disease_name": "pancreatic endocrine carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0004114",
                        "disease_name": "astrocytoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0035344",
                        "disease_name": "retinopathy of prematurity",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "EFO_0001075",
                        "disease_name": "ovarian carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "EFO_0009784",
                        "disease_name": "central serous retinopathy",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0035328",
                        "disease_name": "retinal vein occlusion",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "EFO_1001469",
                        "disease_name": "Mantle cell lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:CN204398, umls:C4551687",
                        "disease_name": "soft tissue sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0205898",
                        "disease_name": "Pineoblastoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C2986658",
                        "disease_name": "diffuse intrinsic pontine glioma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0205851",
                        "disease_name": "germ cell tumor",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0001430",
                        "disease_name": "adenoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0008170",
                        "disease_name": "ovarian cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "umls:C0027439",
                        "disease_name": "nasopharyngeal neoplasm",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0035309",
                        "disease_name": "retinopathy",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "MONDO_0001657",
                        "disease_name": "brain cancer",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0005567",
                        "disease_name": "malignant peritoneal mesothelioma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0021632",
                        "disease_name": "primary brain neoplasm",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0153535, umls:C0151779, umls:C0153536",
                        "disease_name": "cutaneous melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0025286",
                        "disease_name": "meningioma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000681",
                        "disease_name": "renal cell carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:CN201941, umls:C0014474",
                        "disease_name": "ependymoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0017601",
                        "disease_name": "glaucoma",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0001418",
                        "disease_name": "adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000558",
                        "disease_name": "Kaposi's sarcoma",
                        "therapeutic_areas": "EFO_0005741:infectious disease, MONDO_0045024:cancer or benign tumor, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0079744",
                        "disease_name": "diffuse large B-cell lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C1292778",
                        "disease_name": "myeloproliferative disorder",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0025290",
                        "disease_name": "aseptic meningitis",
                        "therapeutic_areas": "EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0684337",
                        "disease_name": "peripheral primitive neuroectodermal tumor",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "EFO_0010282",
                        "disease_name": "gastrointestinal disease",
                        "therapeutic_areas": "EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0100096",
                        "disease_name": "COVID-19",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                    },
                    {
                        "disease_id": "umls:C0206687, umls:C1569637",
                        "disease_name": "endometrioid carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0042170",
                        "disease_name": "Vogt-Koyanagi-Harada disease",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0019562",
                        "disease_name": "von Hippel-Lindau disease",
                        "therapeutic_areas": "EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0026764, umls:C0268381",
                        "disease_name": "multiple myeloma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0238033, umls:C0242787, umls:C0242788",
                        "disease_name": "male breast carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0206754",
                        "disease_name": "neuroendocrine neoplasm",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0919267",
                        "disease_name": "ovarian neoplasm",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0011884",
                        "disease_name": "diabetic retinopathy",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, OTAR_0000020:nutritional or metabolic disease, MONDO_0024458:disorder of visual system, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0206687, umls:C0476089",
                        "disease_name": "endometrial carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_1001513",
                        "disease_name": "liver neoplasm",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "EFO_0000691",
                        "disease_name": "sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0238122",
                        "disease_name": "Fallopian Tube Carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023418",
                        "disease_name": "leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0546837, umls:C0496775",
                        "disease_name": "esophageal cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "EFO_0004142",
                        "disease_name": "colorectal neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0153420, umls:C0699791",
                        "disease_name": "gastric carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0009405",
                        "disease_name": "Lynch syndrome",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0024236",
                        "disease_name": "lymphedema",
                        "therapeutic_areas": "EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_0003860",
                        "disease_name": "pancreatic neoplasm",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0024440",
                        "disease_name": "cystoid macular edema",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C0027873",
                        "disease_name": "neuromyelitis optica",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0007102",
                        "disease_name": "malignant colon neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0027092",
                        "disease_name": "Myopia",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "umls:C0018923",
                        "disease_name": "angiosarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "EFO_0000196",
                        "disease_name": "metastatic prostate cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0002974",
                        "disease_name": "cervical cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0039445",
                        "disease_name": "hereditary hemorrhagic telangiectasia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0751690",
                        "disease_name": "malignant peripheral nerve sheath tumor",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:CN244903, umls:C0262584, umls:C0149925",
                        "disease_name": "small cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "EFO_0008626",
                        "disease_name": "vitreous hemorrhage",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "umls:C1168401",
                        "disease_name": "head and neck squamous cell carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0740457",
                        "disease_name": "kidney cancer",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0153452",
                        "disease_name": "gallbladder cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0009375",
                        "disease_name": "colonic neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0334590",
                        "disease_name": "anaplastic oligodendroglioma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0027651, umls:CN236628",
                        "disease_name": "neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0043510",
                        "disease_name": "brain injury",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, OTAR_0000009:injury, poisoning or other complication",
                    },
                    {
                        "disease_id": "umls:C0085109",
                        "disease_name": "corneal neovascularization",
                        "therapeutic_areas": "MONDO_0024458:disorder of visual system",
                    },
                    {
                        "disease_id": "MONDO_0001056",
                        "disease_name": "gastric cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                ],
                [
                    {
                        "disease_id": "umls:C0009375",
                        "disease_name": "colonic neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0238033, umls:C0242787, umls:C0242788",
                        "disease_name": "male breast carcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0007115",
                        "disease_name": "thyroid cancer",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0036421",
                        "disease_name": "systemic scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0009690:urinary system disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0751690",
                        "disease_name": "malignant peripheral nerve sheath tumor",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0036420",
                        "disease_name": "localised scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C1527249",
                        "disease_name": "colorectal cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0206180",
                        "disease_name": "anaplastic large cell lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0032285",
                        "disease_name": "pneumonia",
                        "therapeutic_areas": "EFO_0005741:infectious disease, OTAR_0000010:respiratory or thoracic disease",
                    },
                    {
                        "disease_id": "umls:C2973725",
                        "disease_name": "pulmonary arterial hypertension",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C3693482",
                        "disease_name": "dermatofibrosarcoma protuberans",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:CN971653, umls:C0025202",
                        "disease_name": "melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0004096",
                        "disease_name": "asthma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease",
                    },
                    {
                        "disease_id": "umls:C0023470",
                        "disease_name": "myeloid leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "EFO_0007328",
                        "disease_name": "influenza",
                        "therapeutic_areas": "EFO_0005741:infectious disease, OTAR_0000010:respiratory or thoracic disease",
                    },
                    {
                        "disease_id": "umls:C1168401",
                        "disease_name": "head and neck squamous cell carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0003873",
                        "disease_name": "rheumatoid arthritis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0008487",
                        "disease_name": "chordoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0037286",
                        "disease_name": "skin neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0036421",
                        "disease_name": "systemic scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0009690:urinary system disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "EFO_0000220",
                        "disease_name": "acute lymphoblastic leukemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0026769",
                        "disease_name": "multiple sclerosis",
                        "therapeutic_areas": "EFO_0000540:immune system disease, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0027831",
                        "disease_name": "neurofibromatosis type 1",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "MONDO_0005147",
                        "disease_name": "type 1 diabetes mellitus",
                        "therapeutic_areas": "OTAR_0000020:nutritional or metabolic disease, EFO_0001379:endocrine system disease, EFO_0000540:immune system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0206726",
                        "disease_name": "gliosarcoma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0007115",
                        "disease_name": "thyroid cancer",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:CN244903, umls:C0262584, umls:C0149925",
                        "disease_name": "small cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "umls:C0032463",
                        "disease_name": "polycythemia vera",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0014457",
                        "disease_name": "Eosinophilia",
                        "therapeutic_areas": "EFO_0000651:phenotype",
                    },
                    {
                        "disease_id": "umls:C0018133",
                        "disease_name": "graft versus host disease",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000009:injury, poisoning or other complication",
                    },
                    {
                        "disease_id": "umls:C0025500",
                        "disease_name": "mesothelioma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0079218, umls:C1851124, umls:CN072436",
                        "disease_name": "Desmoid-type fibromatosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023968",
                        "disease_name": "loiasis",
                        "therapeutic_areas": "EFO_0005741:infectious disease, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "MONDO_0100096",
                        "disease_name": "COVID-19",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                    },
                    {
                        "disease_id": "umls:C0023473",
                        "disease_name": "chronic myelogenous leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0006826",
                        "disease_name": "cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0278883",
                        "disease_name": "metastatic melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0007102",
                        "disease_name": "malignant colon neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0005149",
                        "disease_name": "pulmonary hypertension",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "MONDO_0011962",
                        "disease_name": "endometrial cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0034091",
                        "disease_name": "pulmonary venoocclusive disease",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "EFO_0000691",
                        "disease_name": "sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0008170",
                        "disease_name": "ovarian cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "EFO_1000657",
                        "disease_name": "rectum cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0281508",
                        "disease_name": "desmoplastic small round cell tumor",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0008903",
                        "disease_name": "lung cancer",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0007254",
                        "disease_name": "breast cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0238463",
                        "disease_name": "papillary thyroid carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "MONDO_0100096",
                        "disease_name": "COVID-19",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                    },
                    {
                        "disease_id": "umls:C0235974",
                        "disease_name": "pancreatic carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0018364",
                        "disease_name": "malignant epithelial tumor of ovary",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease",
                    },
                    {
                        "disease_id": "MONDO_0002715",
                        "disease_name": "uterine cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0007097",
                        "disease_name": "carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0079218, umls:C1851124, umls:CN072436",
                        "disease_name": "Desmoid-type fibromatosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023493",
                        "disease_name": "T-cell acute lymphoblastic leukemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "MONDO_0007254",
                        "disease_name": "breast cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0280630",
                        "disease_name": "Uterine Carcinosarcoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:CN971653, umls:C0025202",
                        "disease_name": "melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0001056",
                        "disease_name": "gastric cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "EFO_1001471",
                        "disease_name": "Merkel cell skin cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C0041296",
                        "disease_name": "tuberculosis",
                        "therapeutic_areas": "EFO_0005741:infectious disease",
                    },
                    {
                        "disease_id": "EFO_1000359",
                        "disease_name": "Malignant Pancreatic Neoplasm",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, EFO_0009605:pancreas disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0153535, umls:C0151779, umls:C0153536",
                        "disease_name": "cutaneous melanoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "MONDO_0008903",
                        "disease_name": "lung cancer",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0008479",
                        "disease_name": "chondrosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000768",
                        "disease_name": "idiopathic pulmonary fibrosis",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0007131",
                        "disease_name": "non-small cell lung carcinoma",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C2237660",
                        "disease_name": "wet macular degeneration",
                        "therapeutic_areas": "MONDO_0002025:psychiatric disorder, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C2973725",
                        "disease_name": "pulmonary arterial hypertension",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "MONDO_0000870",
                        "disease_name": "childhood acute lymphoblastic leukemia",
                        "therapeutic_areas": "EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0278996",
                        "disease_name": "head and neck malignant neoplasia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1621958, umls:C0017636, umls:CN227279",
                        "disease_name": "glioblastoma multiforme",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "EFO_0000497",
                        "disease_name": "fibromatosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000094",
                        "disease_name": "B-cell acute lymphoblastic leukemia",
                        "therapeutic_areas": "EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0024305",
                        "disease_name": "non-Hodgkins lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0023890",
                        "disease_name": "cirrhosis of liver",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0027651, umls:CN236628",
                        "disease_name": "neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0019034, umls:C0002895",
                        "disease_name": "sickle cell anemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0238198",
                        "disease_name": "gastrointestinal stromal tumor",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0016057",
                        "disease_name": "fibrosarcoma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0023473",
                        "disease_name": "chronic myelogenous leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0024535",
                        "disease_name": "Plasmodium falciparum malaria",
                        "therapeutic_areas": "EFO_0001379:endocrine system disease, EFO_0005741:infectious disease, EFO_0005803:hematologic disease, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0376358",
                        "disease_name": "prostate cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0006826",
                        "disease_name": "cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0740457",
                        "disease_name": "kidney cancer",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0376544",
                        "disease_name": "hematopoietic and lymphoid cell neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_1000158",
                        "disease_name": "Central Nervous System Neoplasm",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0018133",
                        "disease_name": "graft versus host disease",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000009:injury, poisoning or other complication",
                    },
                    {
                        "disease_id": "umls:C0009402, umls:CN221574",
                        "disease_name": "colorectal carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0153467",
                        "disease_name": "peritoneum cancer",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0079748",
                        "disease_name": "lymphoblastic lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "MONDO_0005149",
                        "disease_name": "pulmonary hypertension",
                        "therapeutic_areas": "EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C1621958, umls:C0017636, umls:CN227279",
                        "disease_name": "glioblastoma multiforme",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0000618:nervous system disease",
                    },
                    {
                        "disease_id": "umls:C0238198",
                        "disease_name": "gastrointestinal stromal tumor",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "MONDO_0011705",
                        "disease_name": "lymphangioleiomyomatosis",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "MONDO_0002158",
                        "disease_name": "fallopian tube cancer",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0278701",
                        "disease_name": "gastric adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C1540912",
                        "disease_name": "Hypereosinophilic syndrome",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, EFO_0000540:immune system disease, EFO_0005803:hematologic disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "EFO_0000712",
                        "disease_name": "stroke",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "EFO_1001814",
                        "disease_name": "nephrogenic fibrosing dermopathy",
                        "therapeutic_areas": "EFO_0010285:integumentary system disease",
                    },
                    {
                        "disease_id": "umls:C1540912",
                        "disease_name": "Hypereosinophilic syndrome",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, EFO_0000540:immune system disease, EFO_0005803:hematologic disease, EFO_0000319:cardiovascular disease",
                    },
                    {
                        "disease_id": "umls:C0007112",
                        "disease_name": "prostate adenocarcinoma",
                        "therapeutic_areas": "OTAR_0000017:reproductive system or breast disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0000658",
                        "disease_name": "plexiform neurofibroma",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C1322286, umls:C0205969, umls:CN207411",
                        "disease_name": "Thymic Carcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0001379:endocrine system disease, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_0000691",
                        "disease_name": "sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0011644",
                        "disease_name": "scleroderma",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0024458:disorder of visual system, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0948968, umls:C0001815, umls:C2355576",
                        "disease_name": "primary myelofibrosis",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "umls:C0699790, umls:C1319315",
                        "disease_name": "colorectal adenocarcinoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0023467",
                        "disease_name": "acute myeloid leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                    {
                        "disease_id": "EFO_0004610",
                        "disease_name": "acute lung injury",
                        "therapeutic_areas": "OTAR_0000010:respiratory or thoracic disease, OTAR_0000009:injury, poisoning or other complication",
                    },
                    {
                        "disease_id": "umls:C0039101",
                        "disease_name": "synovial sarcoma",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "EFO_0003956",
                        "disease_name": "seasonal allergic rhinitis",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000010:respiratory or thoracic disease",
                    },
                    {
                        "disease_id": "MONDO_0044917",
                        "disease_name": "T-lymphoblastic lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C1292778",
                        "disease_name": "chronic myeloproliferative disorder",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0027651, umls:CN236628",
                        "disease_name": "neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0867389",
                        "disease_name": "chronic graft versus host disease",
                        "therapeutic_areas": "EFO_0000540:immune system disease, OTAR_0000009:injury, poisoning or other complication",
                    },
                    {
                        "disease_id": "umls:C0023418",
                        "disease_name": "leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0024305",
                        "disease_name": "non-Hodgkins lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_1001919",
                        "disease_name": "Spinal cord injury",
                        "therapeutic_areas": "EFO_0000618:nervous system disease, OTAR_0000009:injury, poisoning or other complication",
                    },
                    {
                        "disease_id": "umls:C0152271, umls:C0023448",
                        "disease_name": "lymphoid leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_0000220",
                        "disease_name": "acute lymphoblastic leukemia",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "EFO_0000349",
                        "disease_name": "clear cell renal carcinoma",
                        "therapeutic_areas": "EFO_0009690:urinary system disease, MONDO_0045024:cancer or benign tumor",
                    },
                    {
                        "disease_id": "umls:C0079748",
                        "disease_name": "lymphoblastic lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0023418",
                        "disease_name": "leukemia",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0024299",
                        "disease_name": "lymphoma",
                        "therapeutic_areas": "OTAR_0000018:genetic, familial or congenital disease, MONDO_0045024:cancer or benign tumor, EFO_0005803:hematologic disease",
                    },
                    {
                        "disease_id": "umls:C0038356",
                        "disease_name": "stomach neoplasm",
                        "therapeutic_areas": "MONDO_0045024:cancer or benign tumor, EFO_0010282:gastrointestinal disease",
                    },
                    {
                        "disease_id": "umls:C0023467",
                        "disease_name": "acute myeloid leukemia",
                        "therapeutic_areas": "OTAR_0000006:musculoskeletal or connective tissue disease, EFO_0005803:hematologic disease, MONDO_0045024:cancer or benign tumor, EFO_0000540:immune system disease",
                    },
                ],
            ]
        )
        expected_data.name = OPENTARGETS_DISEASE_COL

        pd.testing.assert_series_equal(obtained_data[OPENTARGETS_DISEASE_COL], expected_data)
        self.assertIsInstance(metadata, dict)
