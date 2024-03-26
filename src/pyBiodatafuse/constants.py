# coding: utf-8

"""Python constant file."""


from typing import Union

# Endpoints / API
BGEE_ENDPOINT = "https://www.bgee.org/sparql/"
DISGENET_ENDPOINT = "http://rdf.disgenet.org/sparql/"
MINERVA_ENDPOINT = "https://minerva-net.lcsb.uni.lu/api/"
MOLMEDB_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/molmedb"
OPENTARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
PUBCHEM_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/idsm"
STRING_ENDPOINT = "https://string-db.org/api"
WIKIDATA_ENDPOINT = "https://query.wikidata.org/sparql"
WIKIPATHWAYS_ENDPOINT = "https://sparql.wikipathways.org/sparql"

# Data sources
BGEE = "Bgee"
DISGENET = "DisGeNET"
MINERVA = "MINERVA"
MOLMEDB = "MolMeDB"
OPENTARGETS = "OpenTargets"
PUBCHEM = "PubChem"
STRING = "StringDB"
WIKIDATA = "Wikidata"
WIKIPATHWAYS = "WikiPathways"

# Input type for each data source
BGEE_INPUT_ID = "Ensembl"
DISGENET_INPUT_ID = "NCBI Gene"
MINERVA_INPUT_ID = "Ensembl"
MOLMEDB_GENE_INPUT_ID = "Uniprot-TrEMBL"
MOLMEDB_COMPOUND_INPUT_ID = "InChIKey"
OPENTARGETS_INPUT_ID = "Ensembl"
PUBCHEM_INPUT_ID = "Uniprot-TrEMBL"
STRING_INPUT_ID = "Ensembl"
WIKIDATA_INPUT_ID = ""  # TODO
WIKIPATHWAYS_INPUT_ID = "NCBI Gene"

# Output annotation for each data source
# Anatomical entity node
# Bgee
BGEE_OUTPUT_DICT = {
    "anatomical_entity_id": str,
    "anatomical_entity_name": str,
    "developmental_stage_id": str,
    "developmental_stage_name": str,
    "expression_level": float,
    "confidence_level_id": str,
    "confidence_level_name": str,
}
ANATOMICAL_ENTITY_ID = "UBERON"
DEVELOPMENTAL_STAGE_ID = "HsapDv|UBERON"
CONFIDENCE_LEVEL_ID = "CIO"

# Location node
# Open Targets - Location
OPENTARGETS_LOCATION_OUTPUT_DICT = {
    "location_id": str,
    "location": str,
    "subcellular_location": str,
}
LOCATION_ID = "SL"
OPENTARGETS_LOCATION_COL = f"{OPENTARGETS}_Location"

# Disease node
# DisGeNet
DISGENET_OUTPUT_DICT = {
    "disease_id": str,
    "disease_name": str,
    "score": float,
    "evidence_source": str,
}
# Open Targets - Disease
OPENTARGETS_DISEASE_OUTPUT_DICT = {"disease_id": str, "disease_name": str, "therapeutic_areas": str}
DISEASE_ID = "umls|EFO|MONDO"
OPENTARGETS_DISEASE_COL = f"{OPENTARGETS}_Diseases"

# Pathway node
# MINERVA
MINERVA_OUTPUT_DICT = {
    "pathway_id": int,
    "pathway_label": str,
    "pathway_gene_count": int,
}
# WikiPathways
WIKIPATHWAYS_OUTPUT_DICT = {"pathway_id": str, "pathway_label": str, "pathway_gene_count": int}
# Open Targets - Reactome
OPENTARGETS_REACTOME_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
}
PATHWAY_ID = "WP|R-"
OPENTARGETS_REACTOME_COL = f"{OPENTARGETS}_Reactome"

# GO
# Open Targets - GO processes
OPENTARGETS_GO_OUTPUT_DICT = {"go_id": str, "go_name": str, "go_type": str}
GO_ID = "GO"
OPENTARGETS_GO_COL = f"{OPENTARGETS}_GO"  # TODO: Cross-check if correct name

# Compound
# Open Targets - Compound
OPENTARGETS_COMPOUND_OUTPUT_DICT = {
    "chembl_id": str,
    "drugbank_id": Union[str, None, float],
    "compound_cid": Union[str, None, float],
    "compound_name": str,
    "is_approved": bool,
    "relation": str,
    "adverse_effect_count": Union[int, None, float],
    "adverse_effect": Union[str, None, float],
}
CHEMBL_ID = "CHEMBL"
RELATION = "inhibits|activates"
OPENTARGETS_COMPOUND_COL = f"{OPENTARGETS}_Compounds"  # TODO: Cross-check if correct name
# MolMeDB - Gene input
MOLMEDB_GENE_OUTPUT_DICT = {
    "compound_name": str,
    "InChIKey": str,
    "SMILES": str,
    "compound_cid": str,
    "molmedb_id": str,
    "source_doi": str,
    "source_pmid": str,
    "chebi_id": str,
    "pdb_ligand_id": str,
    "drugbank_id": str,
}
MOLMEDB_ID = "MM"
SOURCE_DOI = "doi"
DRUGBANK_ID = "DB"
MOLMEDB_INHIBITOR_COL = f"{MOLMEDB}_transporter_inhibitor"
# MolMeDB - Compound input
MOLMEDB_COMPOUND_OUTPUT_DICT = {
    "uniprot_trembl_id": str,
    "hgnc_symbol": str,
    "source_doi": str,
    "source_pmid": str,
}
UNIPROT_TREMBL_ID = "P"
MOLMEDB_INHIBITED_COL = f"{MOLMEDB}_transporter_inhibited"
# PubChem - Assays
PUBCHEM_OUTPUT_DICT = {
    "pubchem_assay_id": str,
    "assay_type": str,
    "outcome": str,
    "compound_cid": str,
    "compound_name": str,
    "SMILES": str,
    "InChI": str,
}
OUTCOME = "active|inactive"
INCHI = "InChI"
PUBCHEM_ASSAYS_COL = f"{PUBCHEM}_Assays"

# Gene Node
# STRING
# TODO: to be checked

# Wikidata
# TODO: to be checked

# Node and edge main lable and attributes for each data source
# Anatomical entity node
# Bgee
BGEE_NODE_LABELS = "Anatomical Entity"
BGEE_NODE_MAIN_LABEL = "anatomical_entity_id"
BGEE_NODE_ATTRS = {
    "source": BGEE,
    "name": None,
    "id": None,
    "developmental_stage_name": None,
    "developmental_stage_id": None,
    "labels": BGEE_NODE_LABELS,
}
BGEE_EDGE_LABEL = "expressed_in"
BGEE_EDGE_ATTRS = {
    "source": BGEE,
    "confidence_level_name": None,
    "confidence_level_id": None,
    "label": BGEE_EDGE_LABEL,
}

# Location node
# Open Targets - Location
OPENTARGETS_LOCATION_NODE_LABELS = "Location"
OPENTARGETS_LOCATION_NODE_MAIN_LABEL = "location_id"
OPENTARGETS_LOCATION_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "subcellular_location": None,
    "labels": OPENTARGETS_LOCATION_NODE_LABELS,
}
OPENTARGETS_LOCATION_EDGE_LABEL = "localized_in"
OPENTARGETS_LOCATION_EDGE_ATTRS = {"source": OPENTARGETS, "label": OPENTARGETS_LOCATION_EDGE_LABEL}

# Disease node
# DisGeNet
DISGENET_NODE_LABELS = "Disease"
DISGENET_NODE_MAIN_LABEL = "disease_id"
DISGENET_NODE_ATTRS = {
    "source": DISGENET,
    "name": None,
    "id": None,
    "evidence_source": None,
    "labels": DISGENET_NODE_LABELS,
}
DISGENET_EDGE_LABEL = "associated_with"
DISGENET_EDGE_ATTRS = {
    "source": DISGENET,
    "score": None,
    "label": DISGENET_EDGE_LABEL,
}

# Open Targets - Disease
OPENTARGETS_DISEASE_NODE_LABELS = "Disease"
OPENTARGETS_DISEASE_NODE_MAIN_LABEL = "disease_id"
OPENTARGETS_DISEASE_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "therapeutic_areas": None,
    "labels": OPENTARGETS_DISEASE_NODE_LABELS,
}
OPENTARGETS_DISEASE_EDGE_LABEL = "associated_with"
OPENTARGETS_DISEASE_EDGE_ATTRS = {"source": OPENTARGETS, "label": OPENTARGETS_DISEASE_EDGE_LABEL}

# Pathway node
# MINERVA
MINERVA_NODE_LABELS = "Pathway"
MINERVA_NODE_MAIN_LABEL = "pathway_id"
MINERVA_NODE_ATTRS = {
    "source": MINERVA,
    "name": None,
    "id": None,
    "gene_count": None,
    "labels": MINERVA_NODE_LABELS,
}
MOLMEDB_NODE_MAIN_LABEL = "Compound"
MINERVA_EDGE_LABEL = "part_of"
MINERVA_EDGE_ATTRS = {"source": MINERVA, "label": MINERVA_EDGE_LABEL}

# WikiPathways
WIKIPATHWAYS_NODE_LABELS = "Pathway"
WIKIPATHWAYS_NODE_MAIN_LABEL = "pathway_id"
WIKIPATHWAYS_NODE_ATTRS = {
    "source": WIKIPATHWAYS,
    "name": None,
    "id": None,
    "gene_count": None,
    "labels": WIKIPATHWAYS_NODE_LABELS,
}
WIKIPATHWAYS_EDGE_LABEL = "part_of"
WIKIPATHWAYS_EDGE_ATTRS = {"source": WIKIPATHWAYS, "label": WIKIPATHWAYS_EDGE_LABEL}

# Open Targets - Reactome
OPENTARGETS_REACTOME_NODE_LABELS = "Pathway"
OPENTARGETS_REACTOME_NODE_MAIN_LABEL = "pathway_id"
OPENTARGETS_REACTOME_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "labels": OPENTARGETS_REACTOME_NODE_LABELS,
}
OPENTARGETS_REACTOME_EDGE_LABEL = "part_of"
OPENTARGETS_REACTOME_EDGE_ATTRS = {"source": OPENTARGETS, "label": OPENTARGETS_REACTOME_EDGE_LABEL}

# GO
# Open Targets - GO processes
OPENTARGETS_GO_NODE_LABELS = "Gene Ontology"
OPENTARGETS_GO_NODE_MAIN_LABEL = "go_id"
OPENTARGETS_GO_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "type": None,
    "labels": OPENTARGETS_GO_NODE_LABELS,
}
OPENTARGETS_GO_EDGE_LABEL = "part_of"
OPENTARGETS_GO_EDGE_ATTRS = {"source": OPENTARGETS, "label": OPENTARGETS_GO_EDGE_LABEL}

# Compound
# Open Targets - Compound
OPENTARGETS_COMPOUND_NODE_LABELS = "Compound"
OPENTARGETS_COMPOUND_NODE_MAIN_LABEL = "compound_cid"
OPENTARGETS_COMPOUND_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "chembl_id": None,
    "DrugBank_id": None,
    "compound_cid": None,
    "is_approved": None,
    "adverse_effect_count": None,
    "labels": OPENTARGETS_COMPOUND_NODE_LABELS,
}
OPENTARGETS_COMPOUND_EDGE_ATTRS = {"source": OPENTARGETS, "label": None}
# Side effect
# Open Targets - Compound
OPENTARGETS_SIDE_EFFECT_NODE_LABELS = "Side Effect"
OPENTARGETS_SIDE_EFFECT_NODE_MAIN_LABEL = "adverse_effect"
# TODO: add the side effect id (adverse_effect id)
OPENTARGETS_SIDE_EFFECT_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "labels": OPENTARGETS_SIDE_EFFECT_NODE_MAIN_LABEL,
}
OPENTARGETS_SIDE_EFFECT_EDGE_LABEL = "has_side_effect"

# MolMeDB - Gene input
MOLMEDB_COMPOUND_NODE_LABELS = "Compound"
MOLMEDB_COMPOUND_NODE_MAIN_LABEL = "compound_cid"
MOLMEDB_COMPOUND_NODE_ATTRS = {
    "source": MOLMEDB,
    "name": None,
    "id": None,
    "MolMeDB_id": None,
    "ChEBI_id": None,
    "drugbank_id": None,
    "compound_cid": None,
    "pdb_ligand_id": None,
    "InChIKey": None,
    "SMILES": None,
    "source_doi": None,
    "source_pmid": None,
    "labels": MOLMEDB_COMPOUND_NODE_LABELS,
}
MOLMEDB_COMPOUND_EDGE_LABEL = "inhibits"
MOLMEDB_COMPOUND_EDGE_ATTRS = {
    "source": MOLMEDB,
    "label": MOLMEDB_COMPOUND_EDGE_LABEL,
}

# PubChem - Assays
PUBCHEM_NODE_LABELS = "Compound"
PUBCHEM_NODE_MAIN_LABEL = "compound_cid"
PUBCHEM_NODE_ATTRS = {
    "source": PUBCHEM,
    "name": None,
    "id": None,
    "InChI": None,
    "SMILES": None,
    "labels": PUBCHEM_NODE_LABELS,
}
PUBCHEM_EDGE_ATTRS = {
    "source": PUBCHEM,
    "assay_type": None,
    "pubchem_assay_id": None,
    "outcome": None,
    "label": None,
}

# MolMeDB - Compound input
# TODO: to be checked

# Gene Node
# STRING
# TODO: to be checked
STRING_EDGE_MAIN_LABEL = "stringdb_link_to"
STRING_EDGE_LABEL = "StringDB_ppi_interaction"

# Assay node

# Wikidata
# TODO: to be checked
