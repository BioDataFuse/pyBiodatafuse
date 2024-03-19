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
WIKIPATHWAY_ENDPOINT = "https://sparql.wikipathways.org/sparql"

# Data sources
BGEE = "Bgee"
DISGENET = "DisGeNET"
MINERVA = "MINERVA Net"
MOLMEDB = "MolMeDB"
OPENTARGETS = "OpenTargets"
PUBCHEM = "PubChem"
STRING = "StringDB"
WIKIDATA = "Wikidata"
WIKIPATHWAY = "WikiPathways"

# Input type for each data source
BGEE_INPUT_ID = "Ensembl"
DISGENET_INPUT_ID = "NCBI Gene"
MINERVA_INPUT_ID = "NCBI Gene"
MOLMEDB_GENE_INPUT_ID = "Uniprot-TrEMBL"
MOLMEDB_COMPOUND_INPUT_ID = "InChIKey"
OPENTARGETS_INPUT_ID = "Ensembl"
PUBCHEM_INPUT_ID = "Uniprot-TrEMBL"
STRING_INPUT_ID = "Ensembl"
WIKIDATA_INPUT_ID = ""
WIKIPATHWAY_INPUT_ID = "NCBI Gene"

# Output annotation for each data source
## Anatomical entity node
### Bgee
BGEE_OUTPUT_DICT = {
    "anatomical_entity_id": str,
    "anatomical_entity_name": str,
    "developmental_stage_id": str,
    "developmental_stage_name": str,
    "expression_level": float,
    "confidence_level_id": str,
}
ANATOMICAL_ENTITY_ID = "UBERON"
DEVELOPMENTAL_STAGE_ID = ("HsapDv", "UBERON")
CONFIDENCE_LEVEL_ID = "CIO"
## Location node
### Open Targets - Location
OPENTARGETS_LOCATION_OUTPUT_DICT = {
  "location_id": str,
  "location": str, 
  "subcellular_location": str
}
LOCATION_ID = "SL"
## Disease node
### DisGeNet
DISGENET_OUTPUT_DICT = {
    "disease_id": str,
    "disease_name": str,
    "score": float,
    "source": str
}
### Open Targets - Disease
OPENTARGETS_DISEASE_OUTPUT_DICT = {
  "disease_id": str,
  "disease_name": str,
  "therapeutic_areas": str
}
DISEASE_ID = "umls|EFO|MONDO"
## Pathway node
### MINERVA
MINERVA_OUTPUT_DICT = {
    "pathway_id": int,
    "pathway_label": str,
    "pathway_gene_count": int,
}
### WikiPathways
WIKIPATHWAY_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
    "pathway_gene_count": int}
### Open Targets - Reactome
OPENTARGETS_REACTOME_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
}
PATHWAY_ID = "WP|R-"
## GO
### Open Targets - GO processes
OPENTARGETS_GO_OUTPUT_DICT = {
    "go_id": str,
    "go_name": str
}
GO_ID = "GO"
## Compound
### Open Targets - Compound
OPENTARGETS_COMPOUND_OUTPUT_DICT = {
    "chembl_id": str,
    "drug_name": str,
    "is_approved": bool,
    "relation": str
}
CHEMBL_ID = "CHEMBL"
RELATION = "inhibits|activates"
### MolMeDB - Gene input
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
    "drugbank_id": str
}
MOLMEDB_ID = "MM"
SOURCE_DOI = "doi"
DRUGBANK_ID = "DB"
### MolMeDB - Compound input
MOLMEDB_COMPOUND_OUTPUT_DICT = {
    "uniprot_trembl_id": str,
    "hgnc_symbol": str,
    "source_doi": str,
    "source_pmid": str  
}
UNIPROT_TREMBL_ID = "P"
## Gene Node
### STRING
# TODO: to be checked 

## Assay node
PUBCHEM_OUTPUT_DICT = {
    "assay_type": str,
    "outcome": str,
    "compound_cid": str,
    "compound_name": str,
    "SMILES": str,
    "InChI": str
}
OUTCOME = "active|inactive"
INCHI = "InChI"

### Wikidata
# TODO: to be checked 

# Node type output for each data source
BGEE_NODE_TYPE = ""
DISGENET_NODE_TYPE = ""
MINERVA_NODE_TYPE = ""
MOLMEDB_NODE_TYPE = ""
OPENTARGETS_NODE_TYPE = ""
PUBCHEM_NODE_TYPE = ""
STRING_NODE_TYPE = ""
WIKIDATA_NODE_TYPE = ""
WIKIPATHWAY_NODE_TYPE = ""

# Edge type output for each data source
BGEE_EDGE_TYPE = ""
DISGENET_EDGE_TYPE = ""
MINERVA_EDGE_TYPE = ""
MOLMEDB_EDGE_TYPE = ""
OPENTARGETS_EDGE_TYPE = ""
PUBCHEM_EDGE_TYPE = ""
STRING_EDGE_TYPE = ""
WIKIDATA_EDGE_TYPE = ""
WIKIPATHWAY_EDGE_TYPE = ""
