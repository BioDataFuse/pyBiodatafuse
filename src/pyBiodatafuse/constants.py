# coding: utf-8

"""Python constant file."""


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
MINERVA_INPUT_ID = "NCBI Gene"
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
DISGENET_OUTPUT_DICT = {"disease_id": str, "disease_name": str, "score": float, "evidence_source": str}
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
OPENTARGETS_GO_OUTPUT_DICT = {"go_id": str, "go_name": str}
GO_ID = "GO"
OPENTARGETS_GO_COL = f"{OPENTARGETS}_GO"  # TODO: Cross-check if correct name

# Compound
# Open Targets - Compound
OPENTARGETS_COMPOUND_OUTPUT_DICT = {
    "chembl_id": str,
    "compound_name": str,
    "is_approved": bool,
    "relation": str,
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

# Gene Node
# STRING
# TODO: to be checked

# Assay node
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
PUBCHEM_Assays_COL = f"{PUBCHEM}_Assays"

# Wikidata
# TODO: to be checked

# Node lables for each data source
BGEE_NODE_LABELS = "Anatomical Entity"
DISGENET_NODE_LABELS = "Disease"
MINERVA_NODE_LABELS = "Pathway"
MOLMEDB_NODE_LABELS = "Compound"
OPENTARGETS_REACTOME_NODE_LABELS = "Pathway"
OPENTARGETS_GO_NODE_LABELS = "Gene Ontology"
OPENTARGETS_LOCATION_NODE_LABELS = "Location"
OPENTARGETS_DISEASE_NODE_LABELS = "Disease"
OPENTARGETS_COMPOUND_NODE_LABELS = "Compound"
PUBCHEM_NODE_LABELS = "Compound"
WIKIPATHWAYS_NODE_LABELS = "Pathway"

# Edge label output for each data source
BGEE_EDGE_LABEL = "expressed_in"
DISGENET_EDGE_LABEL = "associated_with"
MINERVA_EDGE_LABEL = "part_of"
MOLMEDB_EDGE_LABEL = "inhibits"
OPENTARGETS_REACTOME_EDGE_LABEL = "part_of"
OPENTARGETS_GO_EDGE_LABEL = "part_of"
OPENTARGETS_LOCATION_EDGE_LABEL = "localized_in"
OPENTARGETS_DISEASE_EDGE_LABEL = "associated_with"
STRING_EDGE_LABEL = "StringDB_ppi_interaction"
WIKIPATHWAYS_EDGE_LABEL = "part_of"
