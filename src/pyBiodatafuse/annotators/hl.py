#!/usr/bin/env python
# -*- coding: utf-8 -*-

import requests

ENSEMBL_REST_ENDPOINT = "https://rest.ensembl.org"

def get_human_homologs(mouse_genes):
    """Retrieve human homologs for mouse genes using ensembl API.
    
    :param mouse_genes: list of gene ids.
    :returns: dictionary mapping mouse genes to human homologs.
    """
    homologs = []

    for mouse_gene in mouse_genes:
        response = requests.get(
            f"{ENSEMBL_REST_ENDPOINT}/homology/id/mouse/{mouse_gene}",
            headers={"Content-Type": "application/json"},
            params={"target_species": "homo_sapiens"}
        )
        
        if response.status_code == 200:
            data = response.json()
            for homology in data.get("data", [])[0].get("homologies", []):
                if homology["target"]["species"] == "homo_sapiens":
                    homologs.append(homology["target"]["id"])
                    break

    return homologs
