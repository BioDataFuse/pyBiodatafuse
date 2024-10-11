#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the DISGENET annotator."""

import json
import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
from requests import Session

from pyBiodatafuse.annotators import disgenet
from pyBiodatafuse.annotators.disgenet import get_gene_disease, get_version_disgenet
from pyBiodatafuse.constants import DISGENET_DISEASE_COL

data_file_folder = os.path.join(os.path.dirname(__file__), "data")


class TestDisgenet(unittest.TestCase):
    """Test the DISGENET class."""

    @patch.object(Session, "get")
    def test_get_version_disgenet(self, mock_sparql_version):
        """Test that the API endpoint returns the expected get_version_bgee results."""
        mock_sparql_version.side_effect = {
            "status": "OK",
            "payload": {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"},
        }

        obtained_version = get_version_disgenet(api_key="test")

        expected_version = {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"}

        assert obtained_version == expected_version

    @patch.object(Session, "get")
    def test_get_gene_disease(self, mock_post_gene_disease):
        """Test the get_gene_disease function."""
        disgenet.get_version_disgenet = Mock(
            return_value={
                "status": "OK",
                "payload": {"lastUpdate": "10 Jul 2024", "version": "DISGENET v24.2"},
            }
        )  # Mock the version call
        disgenet.check_endpoint_disgenet = Mock(return_value=True)

        bridgedb_dataframe = pd.DataFrame(
            {
                "identifier": ["ALG14"],
                "identifier.source": ["HGNC"],
                "target": ["199857"],
                "target.source": ["NCBI Gene"],
            }
        )

        mock_post_gene_disease.return_value.ok = True
        mock_post_gene_disease.return_value.text = '{"status":"OK","paging":{"pageSize":100,"totalElements":6,"totalElementsInPage":6,"currentPageNumber":0},"warnings":["gene_ensembl_id > The parameter value is empty / unspecified.","gene_symbol > The parameter value is empty / unspecified.","disease > The parameter value is empty / unspecified.","chemical > The parameter value is empty / unspecified.","dis_type > The parameter value is empty / unspecified.","dis_class_list > The parameter value is empty / unspecified.","min_score > The parameter value is unspecified or not a double number.","max_score > The parameter value is unspecified or not a double number.","min_ei > The parameter value is unspecified or not a double number.","max_ei > The parameter value is unspecified or not a double number.","min_dsi > The parameter value is unspecified or not a double number.","max_dsi > The parameter value is unspecified or not a double number.","min_dpi > The parameter value is unspecified or not a double number.","max_dpi > The parameter value is unspecified or not a double number.","min_pli > The parameter value is unspecified or not a double number.","max_pli > The parameter value is unspecified or not a double number."],"requestpar":{"gene_ncbi_id":"199857","source":"CURATED","page_number":0},"userinfo":{"profile":"ACADEMIC","profileDescription":"Your DISGENET profile is ACADEMIC; you can access DISGENET REST endpoints with restrictions."},"payload":[{"assocID":15912171,"symbolOfGene":"ALG14","geneNcbiID":199857,"geneEnsemblIDs":["ENSG00000172339"],"geneNcbiType":"protein-coding","geneDSI":0.622,"geneDPI":0.542,"genepLI":6.9786E-5,"geneProteinStrIDs":["Q96F25"],"geneProteinClassIDs":[],"geneProteinClassNames":[],"diseaseVocabularies":["MESH_D018981","NCI_C84615","ORDO_137","DO_5212","MONDO_0015286","UMLS_C0282577"],"diseaseName":"Congenital Disorders of Glycosylation","diseaseType":"disease","diseaseUMLSCUI":"C0282577","diseaseClasses_MSH":["Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)","Nutritional and Metabolic Diseases (C18)"],"diseaseClasses_UMLS_ST":["Disease or Syndrome (T047)"],"diseaseClasses_DO":["physical disorder (0080015)","genetic disease (630)","disease of metabolism (0014667)"],"diseaseClasses_HPO":[],"chemicalsIncludedInEvidence":null,"numberPmidsWithChemsIncludedInEvidenceBySource":[],"score":0.6000000000000001,"yearInitial":2005,"yearFinal":2022,"el":"Limited","ei":1.0},{"assocID":20783606,"symbolOfGene":"ALG14","geneNcbiID":199857,"geneEnsemblIDs":["ENSG00000172339"],"geneNcbiType":"protein-coding","geneDSI":0.622,"geneDPI":0.542,"genepLI":6.9786E-5,"geneProteinStrIDs":["Q96F25"],"geneProteinClassIDs":[],"geneProteinClassNames":[],"diseaseVocabularies":["MESH_D020294","NCI_C84647","ORDO_590","DO_3635","MONDO_0018940","UMLS_C0751882"],"diseaseName":"Myasthenic Syndromes, Congenital","diseaseType":"disease","diseaseUMLSCUI":"C0751882","diseaseClasses_MSH":["Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)","Nervous System Diseases (C10)"],"diseaseClasses_UMLS_ST":["Disease or Syndrome (T047)"],"diseaseClasses_DO":["physical disorder (0080015)","disease of anatomical entity (7)"],"diseaseClasses_HPO":[],"chemicalsIncludedInEvidence":null,"numberPmidsWithChemsIncludedInEvidenceBySource":[],"score":0.6000000000000001,"yearInitial":null,"yearFinal":null,"el":null,"ei":1.0},{"assocID":25777482,"symbolOfGene":"ALG14","geneNcbiID":199857,"geneEnsemblIDs":["ENSG00000172339"],"geneNcbiType":"protein-coding","geneDSI":0.622,"geneDPI":0.542,"genepLI":6.9786E-5,"geneProteinStrIDs":["Q96F25"],"geneProteinClassIDs":[],"geneProteinClassNames":[],"diseaseVocabularies":["OMIM_612866","OMIM_616227","DO_0110658","MONDO_0014542","UMLS_C4015596"],"diseaseName":"MYASTHENIC SYNDROME, CONGENITAL, 15","diseaseType":"disease","diseaseUMLSCUI":"C4015596","diseaseClasses_MSH":[],"diseaseClasses_UMLS_ST":["Disease or Syndrome (T047)"],"diseaseClasses_DO":["genetic disease (630)","physical disorder (0080015)","disease of anatomical entity (7)"],"diseaseClasses_HPO":[],"chemicalsIncludedInEvidence":null,"numberPmidsWithChemsIncludedInEvidenceBySource":[],"score":0.6,"yearInitial":2013,"yearFinal":2021,"el":null,"ei":1.0},{"assocID":26279064,"symbolOfGene":"ALG14","geneNcbiID":199857,"geneEnsemblIDs":["ENSG00000172339"],"geneNcbiType":"protein-coding","geneDSI":0.622,"geneDPI":0.542,"genepLI":6.9786E-5,"geneProteinStrIDs":["Q96F25"],"geneProteinClassIDs":[],"geneProteinClassNames":[],"diseaseVocabularies":["OMIM_612866","OMIM_619036","MONDO_0033619","UMLS_C5436652"],"diseaseName":"MYOPATHY, EPILEPSY, AND PROGRESSIVE CEREBRAL ATROPHY","diseaseType":"disease","diseaseUMLSCUI":"C5436652","diseaseClasses_MSH":[],"diseaseClasses_UMLS_ST":["Disease or Syndrome (T047)"],"diseaseClasses_DO":[],"diseaseClasses_HPO":[],"chemicalsIncludedInEvidence":null,"numberPmidsWithChemsIncludedInEvidenceBySource":[],"score":0.5,"yearInitial":2017,"yearFinal":2017,"el":null,"ei":1.0},{"assocID":26354111,"symbolOfGene":"ALG14","geneNcbiID":199857,"geneEnsemblIDs":["ENSG00000172339"],"geneNcbiType":"protein-coding","geneDSI":0.622,"geneDPI":0.542,"genepLI":6.9786E-5,"geneProteinStrIDs":["Q96F25"],"geneProteinClassIDs":[],"geneProteinClassNames":[],"diseaseVocabularies":["ORDO_353327","EFO_0700079","UMLS_C5680989"],"diseaseName":"Congenital myasthenic syndromes with glycosylation defect","diseaseType":"disease","diseaseUMLSCUI":"C5680989","diseaseClasses_MSH":["Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)","Nervous System Diseases (C10)"],"diseaseClasses_UMLS_ST":["Disease or Syndrome (T047)"],"diseaseClasses_DO":[],"diseaseClasses_HPO":[],"chemicalsIncludedInEvidence":null,"numberPmidsWithChemsIncludedInEvidenceBySource":[],"score":0.4,"yearInitial":2013,"yearFinal":2013,"el":null,"ei":1.0},{"assocID":26279043,"symbolOfGene":"ALG14","geneNcbiID":199857,"geneEnsemblIDs":["ENSG00000172339"],"geneNcbiType":"protein-coding","geneDSI":0.622,"geneDPI":0.542,"genepLI":6.9786E-5,"geneProteinStrIDs":["Q96F25"],"geneProteinClassIDs":[],"geneProteinClassNames":[],"diseaseVocabularies":["OMIM_612866","OMIM_619031","MONDO_0033572","UMLS_C5436646"],"diseaseName":"INTELLECTUAL DEVELOPMENTAL DISORDER WITH EPILEPSY, BEHAVIORAL ABNORMALITIES, AND COARSE FACIES","diseaseType":"disease","diseaseUMLSCUI":"C5436646","diseaseClasses_MSH":[],"diseaseClasses_UMLS_ST":["Disease or Syndrome (T047)"],"diseaseClasses_DO":[],"diseaseClasses_HPO":[],"chemicalsIncludedInEvidence":null,"numberPmidsWithChemsIncludedInEvidenceBySource":[],"score":0.4,"yearInitial":2018,"yearFinal":2018,"el":null,"ei":1.0}]}'

        obtained_data, metadata = get_gene_disease(api_key="test", bridgedb_df=bridgedb_dataframe)

        expected_data = pd.Series(
            [
                [
                    {
                        "HPO": None,
                        "NCI": "NCI_C84615",
                        "OMIM": "",
                        "MONDO": "MONDO_0015286",
                        "ORDO": "ORDO_137",
                        "EFO": "",
                        "DO": "DO_5212",
                        "MESH": "MESH_D018981",
                        "UMLS": "UMLS_C0282577",
                        "disease_name": "Congenital Disorders of Glycosylation",
                        "disease_type": "disease",
                        "disease_umlscui": "C0282577",
                        "score": 0.6000000000000001,
                        "ei": 1.0,
                        "el": "Limited",
                    },
                    {
                        "HPO": None,
                        "NCI": "NCI_C84647",
                        "OMIM": "",
                        "MONDO": "MONDO_0018940",
                        "ORDO": "ORDO_590",
                        "EFO": "",
                        "DO": "DO_3635",
                        "MESH": "MESH_D020294",
                        "UMLS": "UMLS_C0751882",
                        "disease_name": "Myasthenic Syndromes, Congenital",
                        "disease_type": "disease",
                        "disease_umlscui": "C0751882",
                        "score": 0.6000000000000001,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "HPO": None,
                        "NCI": "",
                        "OMIM": "OMIM_612866, OMIM_616227",
                        "MONDO": "MONDO_0014542",
                        "ORDO": "",
                        "EFO": "",
                        "DO": "DO_0110658",
                        "MESH": "",
                        "UMLS": "UMLS_C4015596",
                        "disease_name": "MYASTHENIC SYNDROME, CONGENITAL, 15",
                        "disease_type": "disease",
                        "disease_umlscui": "C4015596",
                        "score": 0.6,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "HPO": None,
                        "NCI": "",
                        "OMIM": "OMIM_612866, OMIM_619036",
                        "MONDO": "MONDO_0033619",
                        "ORDO": "",
                        "EFO": "",
                        "DO": "",
                        "MESH": "",
                        "UMLS": "UMLS_C5436652",
                        "disease_name": "MYOPATHY, EPILEPSY, AND PROGRESSIVE CEREBRAL ATROPHY",
                        "disease_type": "disease",
                        "disease_umlscui": "C5436652",
                        "score": 0.5,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "HPO": None,
                        "NCI": "",
                        "OMIM": "",
                        "MONDO": "",
                        "ORDO": "ORDO_353327",
                        "EFO": "EFO_0700079",
                        "DO": "",
                        "MESH": "",
                        "UMLS": "UMLS_C5680989",
                        "disease_name": "Congenital myasthenic syndromes with glycosylation defect",
                        "disease_type": "disease",
                        "disease_umlscui": "C5680989",
                        "score": 0.4,
                        "ei": 1.0,
                        "el": None,
                    },
                    {
                        "HPO": None,
                        "NCI": "",
                        "OMIM": "OMIM_612866, OMIM_619031",
                        "MONDO": "MONDO_0033572",
                        "ORDO": "",
                        "EFO": "",
                        "DO": "",
                        "MESH": "",
                        "UMLS": "UMLS_C5436646",
                        "disease_name": "INTELLECTUAL DEVELOPMENTAL DISORDER WITH EPILEPSY, BEHAVIORAL ABNORMALITIES, AND COARSE FACIES",
                        "disease_type": "disease",
                        "disease_umlscui": "C5436646",
                        "score": 0.4,
                        "ei": 1.0,
                        "el": None,
                    },
                ]
            ]
        )
        expected_data.name = DISGENET_DISEASE_COL

        pd.testing.assert_series_equal(obtained_data[DISGENET_DISEASE_COL], expected_data)
        self.assertIsInstance(metadata, dict)
