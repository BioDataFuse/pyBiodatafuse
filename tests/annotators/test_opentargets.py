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
    get_gene_go_process,
    get_gene_compound_interactions,
    get_gene_disease_associations,
    get_gene_reactome_pathways,
    get_version_opentargets,
)
from pyBiodatafuse.constants import OPENTARGETS_GO_COL, OPENTARGETS_REACTOME_COL

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
    def test_get_gene_reactome_pathways(self, mock_post_go_reactome):
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
            mock_post_go_reactome.return_value.json.return_value = json.load(f)

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
    def test_get_gene_compound_interactions(self, mock_post_go_reactome):
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

        # with open(os.path.join(data_file_folder, "opentargets_reactome_mock_data.json")) as f:
        #     mock_post_go_reactome.return_value.json.return_value = json.load(f)

        # obtained_data, metadata = get_gene_reactome_pathways(bridgedb_dataframe_genes)

        # expected_data = pd.Series(
        #     [
        #         [
        #             {
        #                 "pathway_label": "Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein",
        #                 "pathway_id": "R-HSA-446193",
        #             },
        #             {
        #                 "pathway_label": "Defective ALG14 causes ALG14-CMS",
        #                 "pathway_id": "R-HSA-5633231",
        #             },
        #         ],
        #         [
        #             {
        #                 "pathway_label": "Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein",
        #                 "pathway_id": "R-HSA-446193",
        #             },
        #             {
        #                 "pathway_label": "Defective ALG2 causes CDG-1i",
        #                 "pathway_id": "R-HSA-4549349",
        #             },
        #         ],
        #         [
        #             {
        #                 "pathway_label": "Highly calcium permeable nicotinic acetylcholine receptors",
        #                 "pathway_id": "R-HSA-629597",
        #             },
        #             {
        #                 "pathway_label": "Highly calcium permeable postsynaptic nicotinic acetylcholine receptors",
        #                 "pathway_id": "R-HSA-629594",
        #             },
        #         ],
        #     ]
        # )
        # expected_data.name = OPENTARGETS_REACTOME_COL

        # pd.testing.assert_series_equal(obtained_data[OPENTARGETS_REACTOME_COL], expected_data)
        # self.assertIsInstance(metadata, dict)
