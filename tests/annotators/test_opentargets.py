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
from pyBiodatafuse.annotators.opentargets import (  # get_compound_disease_interactions,
    get_disease_compound_interactions,
    get_gene_compound_interactions,
    get_gene_go_process,
    get_gene_reactome_pathways,
    get_version_opentargets,
)
from pyBiodatafuse.constants import (  # OPENTARGETS_DISEASE_COL,
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
                            "CHEMBL1201248",
                            "CHEMBL1200641",
                            "CHEMBL703",
                            "CHEMBL1201244",
                            "CHEMBL1201219",
                            "CHEMBL983",
                            "CHEMBL1200648",
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
                            "pubchem.compound:62887",
                            "pubchem.compound:62886",
                            "pubchem.compound:5314",
                            "pubchem.compound:441290",
                            "pubchem.compound:39765",
                            "pubchem.compound:22475",
                            "pubchem.compound:441351",
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
                        "chembl_id": "CHEMBL:CHEMBL1201248",
                        "drugbank_id": "DrugBank:DB00565",
                        "compound_cid": "pubchem.compound:62887",
                        "compound_name": "CISATRACURIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 102.0,
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
                            {"name": "diarrhoea"},
                            {"name": "drug reaction with eosinophilia and systemic symptoms"},
                            {"name": "tachyphylaxis"},
                            {"name": "myocardial stunning"},
                            {"name": "rash maculo-papular"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1200641",
                        "drugbank_id": None,
                        "compound_cid": "pubchem.compound:62886",
                        "compound_name": "CISATRACURIUM BESYLATE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 83.0,
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
                        "chembl_id": "CHEMBL:CHEMBL703",
                        "drugbank_id": "DrugBank:DB00202",
                        "compound_cid": "pubchem.compound:5314",
                        "compound_name": "SUXAMETHONIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "activates",
                        "adverse_effect_count": 119.0,
                        "adverse_effect": [
                            {"name": "hyperthermia malignant"},
                            {"name": "anaphylactic shock"},
                            {"name": "cardiac arrest"},
                            {"name": "hypotension"},
                            {"name": "anaphylactic reaction"},
                            {"name": "neuromuscular block prolonged"},
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
                            {"name": "renal ischaemia"},
                            {"name": "post procedural complication"},
                            {"name": "diabetes insipidus"},
                            {"name": "cardio-respiratory arrest"},
                            {"name": "left ventricular dysfunction"},
                            {"name": "blood ph decreased"},
                            {"name": "foetal death"},
                            {"name": "muscle contractions involuntary"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201244",
                        "drugbank_id": "DrugBank:DB00728",
                        "compound_cid": "pubchem.compound:441290",
                        "compound_name": "ROCURONIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 98.0,
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
                            {"name": "airway peak pressure increased"},
                            {"name": "negative pressure pulmonary oedema"},
                            {"name": "procedural hypotension"},
                            {"name": "hypoventilation"},
                            {"name": "fatigue"},
                            {"name": "hypoxia"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201219",
                        "drugbank_id": "DrugBank:DB01339",
                        "compound_cid": "pubchem.compound:39765",
                        "compound_name": "VECURONIUM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 27.0,
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
                            {"name": "intervertebral discitis"},
                            {"name": "drug resistance"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL983",
                        "drugbank_id": None,
                        "compound_cid": "pubchem.compound:22475",
                        "compound_name": "SUCCINYLCHOLINE CHLORIDE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "activates",
                        "adverse_effect_count": 121.0,
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
                            {"name": "left ventricular dysfunction"},
                            {"name": "premature baby"},
                            {"name": "apnoea"},
                            {"name": "blood ph decreased"},
                            {"name": "foetal death"},
                            {"name": "bradycardia"},
                            {"name": "po2 increased"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1200648",
                        "drugbank_id": None,
                        "compound_cid": "pubchem.compound:441351",
                        "compound_name": "ROCURONIUM BROMIDE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "inhibits",
                        "adverse_effect_count": 94.0,
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
                            {"name": "diarrhoea"},
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
                            "CHEMBL762",
                            "CHEMBL1161",
                            "CHEMBL787",
                            "CHEMBL70",
                            "CHEMBL480",
                            "CHEMBL596",
                            "CHEMBL1201191",
                            "CHEMBL131",
                            "CHEMBL384467",
                            "CHEMBL1504",
                            "CHEMBL655",
                            "CHEMBL778",
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
                            "name",
                            "name",
                        ],
                        "target": [
                            "pubchem.compound:4636",
                            "pubchem.compound:441336",
                            "pubchem.compound:5281040",
                            "pubchem.compound:5288826",
                            "pubchem.compound:3883",
                            "pubchem.compound:3345",
                            "pubchem.compound:1549000",
                            "pubchem.compound:5755",
                            "pubchem.compound:5743",
                            "pubchem.compound:6436",
                            "pubchem.compound:4192",
                            "pubchem.compound:5311068",
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
                "identifier": ["A2ML1", "ABCA1"],
                "identifier.source": ["HGNC", "HGNC"],
                "target": ["EFO_0004992", "EFO_0700136"],
                "target.source": ["EFO", "EFO"],
            }
        )

        with open(os.path.join(data_file_folder, "opentargets_disease_compound_mock.json")) as f:
            mock_post_compound.return_value.json.return_value = json.load(f)

        obtained_data, metadata = get_disease_compound_interactions(bridge_df)

        expected_data = pd.Series(
            [
                [
                    {
                        "chembl_id": "CHEMBL:CHEMBL762",
                        "drugbank_id": "DrugBank:DB00935",
                        "compound_cid": "pubchem.compound:4636",
                        "compound_name": "OXYMETAZOLINE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 5.0,
                        "adverse_effect": [
                            {"name": "ductus arteriosus stenosis foetal"},
                            {"name": "cerebral vasoconstriction"},
                            {"name": "tricuspid valve incompetence"},
                            {"name": "pupils unequal"},
                            {"name": "mydriasis"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1161",
                        "drugbank_id": "DrugBank:DB14512",
                        "compound_cid": "pubchem.compound:441336",
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
                            {"name": "nasal polyps"},
                            {"name": "poor quality drug administered"},
                            {"name": "rhinitis allergic"},
                            {"name": "sleep disorder due to a general medical condition"},
                            {"name": "anosmia"},
                            {"name": "productive cough"},
                            {"name": "diarrhoea"},
                            {"name": "product availability issue"},
                            {"name": "total lung capacity increased"},
                            {"name": "nasal discomfort"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL787",
                        "drugbank_id": "DrugBank:DB00471",
                        "compound_cid": "pubchem.compound:5281040",
                        "compound_name": "MONTELUKAST",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 65.0,
                        "adverse_effect": [
                            {"name": "asthma"},
                            {"name": "wheezing"},
                            {"name": "allergic granulomatous angiitis"},
                            {"name": "obstructive airways disorder"},
                            {"name": "therapeutic product effect incomplete"},
                            {"name": "suicidal ideation"},
                            {"name": "blood count abnormal"},
                            {"name": "sleep disorder due to a general medical condition"},
                            {"name": "dyspnoea"},
                            {"name": "aggression"},
                            {"name": "nightmare"},
                            {"name": "loss of personal independence in daily activities"},
                            {"name": "abnormal behaviour"},
                            {"name": "cough"},
                            {"name": "sleep terror"},
                            {"name": "productive cough"},
                            {"name": "eosinophilic granulomatosis with polyangiitis"},
                            {"name": "depression"},
                            {"name": "anxiety"},
                            {"name": "sputum discoloured"},
                            {"name": "chest discomfort"},
                            {"name": "nasal polyps"},
                            {"name": "anger"},
                            {"name": "temperature regulation disorder"},
                            {"name": "food allergy"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL70",
                        "drugbank_id": "DrugBank:DB00295",
                        "compound_cid": "pubchem.compound:5288826",
                        "compound_name": "MORPHINE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 42.0,
                        "adverse_effect": [
                            {"name": "respiratory depression"},
                            {"name": "sedation"},
                            {"name": "somnolence"},
                            {"name": "pain"},
                            {"name": "confusional state"},
                            {"name": "drug withdrawal syndrome"},
                            {"name": "hyperaesthesia"},
                            {"name": "mental status changes"},
                            {"name": "drug hypersensitivity"},
                            {"name": "delirium"},
                            {"name": "constipation"},
                            {"name": "unresponsive to stimuli"},
                            {"name": "depressed level of consciousness"},
                            {"name": "hypotension"},
                            {"name": "miosis"},
                            {"name": "contrast media allergy"},
                            {"name": "toxicity to various agents"},
                            {"name": "lethargy"},
                            {"name": "respiratory rate decreased"},
                            {"name": "oxygen saturation decreased"},
                            {"name": "withdrawal syndrome"},
                            {"name": "drug abuse"},
                            {"name": "coma"},
                            {"name": "nausea"},
                            {"name": "myoclonus"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL480",
                        "drugbank_id": "DrugBank:DB00448",
                        "compound_cid": "pubchem.compound:3883",
                        "compound_name": "LANSOPRAZOLE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 22.0,
                        "adverse_effect": [
                            {"name": "chronic kidney disease"},
                            {"name": "end stage renal disease"},
                            {"name": "renal injury"},
                            {"name": "acute kidney injury"},
                            {"name": "hyponatraemia"},
                            {"name": "renal failure"},
                            {"name": "colitis microscopic"},
                            {"name": "hypomagnesaemia"},
                            {"name": "tubulointerstitial nephritis"},
                            {"name": "rebound acid hypersecretion"},
                            {"name": "subacute cutaneous lupus erythematosus"},
                            {"name": "nephrogenic anaemia"},
                            {"name": "hypocalcaemia"},
                            {"name": "gastrooesophageal reflux disease"},
                            {"name": "renal haemangioma"},
                            {"name": "drug reaction with eosinophilia and systemic symptoms"},
                            {"name": "agranulocytosis"},
                            {"name": "hypertensive nephropathy"},
                            {"name": "hepatocellular injury"},
                            {"name": "nephropathy"},
                            {"name": "renal impairment"},
                            {"name": "injection site pain"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL596",
                        "drugbank_id": "DrugBank:DB00813",
                        "compound_cid": "pubchem.compound:3345",
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
                            {"name": "delirium"},
                            {"name": "hypotension"},
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
                            {"name": "cardiac arrest"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL1201191",
                        "drugbank_id": "DrugBank:DB06282",
                        "compound_cid": "pubchem.compound:1549000",
                        "compound_name": "LEVOCETIRIZINE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 25.0,
                        "adverse_effect": [
                            {"name": "drug-induced liver injury"},
                            {"name": "febrile convulsion"},
                            {"name": "alopecia areata"},
                            {"name": "dermatitis contact"},
                            {"name": "allergic reaction to excipient"},
                            {"name": "urticaria"},
                            {"name": "somnolence"},
                            {"name": "angioedema"},
                            {"name": "erythema multiforme"},
                            {"name": "fixed eruption"},
                            {"name": "face oedema"},
                            {"name": "loss of consciousness"},
                            {"name": "dermatitis bullous"},
                            {"name": "congenital inguinal hernia"},
                            {"name": "hypogammaglobulinaemia"},
                            {"name": "anaphylactic shock"},
                            {"name": "swelling face"},
                            {"name": "eosinophilia"},
                            {"name": "altered state of consciousness"},
                            {"name": "bronchopulmonary dysplasia"},
                            {"name": "anaphylactic reaction"},
                            {"name": "haemorrhage subcutaneous"},
                            {"name": "arrhythmia"},
                            {"name": "eye oedema"},
                            {"name": "neonatal infection"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL131",
                        "drugbank_id": "DrugBank:DB00860",
                        "compound_cid": "pubchem.compound:5755",
                        "compound_name": "PREDNISOLONE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 11.0,
                        "adverse_effect": [
                            {"name": "cytomegalovirus infection"},
                            {"name": "pneumocystis jirovecii pneumonia"},
                            {"name": "transplant rejection"},
                            {"name": "febrile neutropenia"},
                            {"name": "polyomavirus-associated nephropathy"},
                            {"name": "post transplant lymphoproliferative disorder"},
                            {"name": "nocardiosis"},
                            {"name": "epstein-barr virus infection"},
                            {"name": "cytomegalovirus viraemia"},
                            {"name": "hepatitis e"},
                            {"name": "bk virus infection"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL384467",
                        "drugbank_id": "DrugBank:DB01234",
                        "compound_cid": "pubchem.compound:5743",
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
                        "compound_cid": "pubchem.compound:6436",
                        "compound_name": "TRIAMCINOLONE ACETONIDE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 69.0,
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
                            {"name": "cytomegalovirus chorioretinitis"},
                            {"name": "wheelchair user"},
                            {"name": "joint destruction"},
                            {"name": "body height below normal"},
                            {"name": "retinal artery occlusion"},
                            {"name": "eye inflammation"},
                            {"name": "polyarthritis"},
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
                        "compound_cid": "pubchem.compound:4192",
                        "compound_name": "MIDAZOLAM",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 57.0,
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
                            {"name": "unresponsive to stimuli"},
                            {"name": "agitation"},
                            {"name": "systemic inflammatory response syndrome"},
                            {"name": "bradycardia"},
                            {"name": "upper airway obstruction"},
                            {"name": "apnoea"},
                            {"name": "myoclonus"},
                        ],
                    },
                    {
                        "chembl_id": "CHEMBL:CHEMBL778",
                        "drugbank_id": "DrugBank:DB00633",
                        "compound_cid": "pubchem.compound:5311068",
                        "compound_name": "DEXMEDETOMIDINE",
                        "clincal_trial_phase": 4.0,
                        "is_approved": True,
                        "relation": "treats",
                        "adverse_effect_count": 94.0,
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
