# coding: utf-8

"""Python constant file."""


# Endpoints / API
BGEE_ENDPOINT = "https://www.bgee.org/sparql/"
DISGENET_ENDPOINT = "https://api.disgenet.com/api/v1/gda/summary"
MINERVA_ENDPOINT = "https://minerva-net.lcsb.uni.lu/api/"
MOLMEDB_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/molmedb"
OPENTARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
PUBCHEM_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/idsm"
STRING_ENDPOINT = "https://string-db.org/api"
WIKIDATA_ENDPOINT = "https://query.wikidata.org/sparql"
WIKIPATHWAYS_ENDPOINT = "https://sparql.wikipathways.org/sparql"

# Data sources
BGEE = "Bgee"
DISGENET = "DISGENET"
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
OPENTARGETS_GENE_INPUT_ID = "Ensembl"
OPENTARGETS_COMPOUND_INPUT_ID = "PubChem Compound"  # If using bridgedb mapping
OPENTARGETS_COMPOUND_QUERY_INPUT_ID = "chembl_id"
OPENTARGETS_DISEASE_INPUT_ID = "EFO"
PUBCHEM_INPUT_ID = "Uniprot-TrEMBL"
STRING_INPUT_ID = "Ensembl"
WIKIDATA_PUBLICATION_INPUT_ID = "NCBI Gene"
WIKIPATHWAYS_INPUT_ID = "NCBI Gene"

# Output annotation for each data source
# Bgee
BGEE_OUTPUT_DICT = {
    "anatomical_entity_id": str,
    "anatomical_entity_name": str,
    "expression_level": float,
    "confidence_level_id": str,
    "confidence_level_name": str,
    "developmental_stage_id": str,
    "developmental_stage_name": str,
}
ANATOMICAL_ENTITY_ID = "UBERON"
DEVELOPMENTAL_STAGE_ID = "HsapDv|UBERON"
CONFIDENCE_LEVEL_ID = "CIO"
ANATOMICAL_ENTITIES_LIST = [
    "blood",
    "bone marrow",
    "brain",
    "breast",
    "cardiovascular system",
    "digestive system",
    "heart",
    "immune organ",
    "kidney",
    "liver",
    "lung",
    "nervous system",
    "pancreas",
    "placenta",
    "reproductive system",
    "respiratory system",
    "skeletal system",
]
BGEE_EXPRESSION_LEVELS_COL = f"{BGEE}_expression_levels"

# DISGENET
DISGENET_OUTPUT_DICT = {
    "disease_name": str,
    "HPO": str,  # "HPO_HP:0100013"
    "NCI": str,  # "NCI_C2910"
    "OMIM": str,  # "OMIM_607906"
    "MONDO": str,  # "MONDO_0021100"
    "ORDO": str,  # "ORDO_137"
    "EFO": str,  # "EFO_0003756"
    "DO": str,  # "DO_0060041"
    "MESH": str,  # "MESH_D000067877"
    "UMLS": str,  # "UMLS_C1510586"
    "disease_type": str,
    "disease_umlscui": str,  # "C1510586"
    "score": float,
    "ei": float,
    "el": float,
}
DISGENET_DISEASE_COL = f"{DISGENET}_diseases"

# Open Targets - Disease
OPENTARGETS_DISEASE_OUTPUT_DICT = {
    "disease_name": str,
    "therapeutic_areas": str,
    "HPO": str,  # "HPO_HP:0100013"
    "NCI": str,  # "NCI_C2910"
    "OMIM": str,  # "OMIM_607906"
    "MONDO": str,  # "MONDO_0021100"
    "ORDO": str,  # "ORDO_137"
    "EFO": str,  # "EFO_0003756"
    "DO": str,  # "DO_0060041"
    "MESH": str,  # "MESH_D000067877"
    "UMLS": str,  # "UMLS_C1510586"
}  # TODO: Tooba please check if you want to add compound annotations too here in the dict
OPENTARGETS_IGNORE_DISEASE_IDS = [
    "icd10cm",
    "icd9",
    "snomedct",
    "sctid",
    "meddra",
    "icd10",
    "wikipedia",
    "snomedct_us",
    "oncotree",
    "nifstd",
    "gard",
    "nord",
    "icdo",
    "hgnc",
    "cohd",
    "kegg",
    "decipher",
    "http",
    "omimps",
    "csp",
]
DRUG_ID = "CHEMBL"  # TODO: Yojana, please check as this is not in the opentarget output dict
OPENTARGETS_DISEASE_COL = f"{OPENTARGETS}_diseases"

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
PATHWAY_ID = "WP|R-"  # ID Start with WP or R-
OPENTARGETS_REACTOME_COL = f"{OPENTARGETS}_reactome"

# Open Targets - GO processes
OPENTARGETS_GO_OUTPUT_DICT = {"go_id": str, "go_name": str, "go_type": str}
GO_ID = "GO"
OPENTARGETS_GO_COL = f"{OPENTARGETS}_go"

# Open Targets - Compound
OPENTARGETS_COMPOUND_OUTPUT_DICT = {
    "chembl_id": str,
    "drugbank_id": str,
    "compound_cid": str,
    "compound_name": str,
    "clincal_trial_phase": int,
    "is_approved": bool,
    "relation": str,
    "adverse_effect_count": int,
    "adverse_effect": dict,
}
CHEMBL_ID = "CHEMBL"
DRUGBANK_ID = "DB"
RELATION = "inhibits|activates"
OPENTARGETS_COMPOUND_COL = f"{OPENTARGETS}_compounds"

# MolMeDB - Gene/Protein input
MOLMEDB_PROTEIN_COMPOUND_OUTPUT_DICT = {
    "compound_name": str,
    "InChIKey": str,
    "SMILES": str,
    "compound_cid": str,
    "molmedb_id": str,
    "source_pmid": str,
    "chebi_id": str,
    "drugbank_id": str,
    "uniprot_trembl_id": str,  # uniprot id of isoform
}
MOLMEDB_ID = "MM"
DRUGBANK_ID = "DB"
MOLMEDB_INHIBITOR_COL = f"{MOLMEDB}_transporter_inhibitor"

# MolMeDB - Compound input
MOLMEDB_COMPOUND_PROTEIN_OUTPUT_DICT = {
    "uniprot_trembl_id": str,
    "hgnc_symbol": str,
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
PUBCHEM_ASSAYS_COL = f"{PUBCHEM}_assays"

# STRING
STRING_OUTPUT_DICT = {"stringdb_link_to": str, STRING_INPUT_ID: str, "score": int}
STRING_PPI_COL = f"{STRING}_ppi"

# Wikidata
WIKIDATA_CC_COL = f"{WIKIDATA}_cellular_components"

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
BGEE_EDGE_LABEL = "expressed_by"
BGEE_EDGE_ATTRS = {
    "source": BGEE,
    "expression_level": None,
    "confidence_level_name": None,
    "confidence_level_id": None,
    "label": BGEE_EDGE_LABEL,
}

# Disease node
# DISGENET
DISGENET_NODE_LABELS = "Disease"
DISGENET_NODE_MAIN_LABEL = "disease_id"
DISGENET_NODE_ATTRS = {
    "source": DISGENET,
    "name": None,
    "id": None,
    "disease_type": None,
    "ei": None,
    "el": None,
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
}  # TODO: Yojana, would it be possible to add pathway size here (gene_count)
OPENTARGETS_REACTOME_EDGE_LABEL = "part_of"
OPENTARGETS_REACTOME_EDGE_ATTRS = {"source": OPENTARGETS, "label": OPENTARGETS_REACTOME_EDGE_LABEL}

# GO nodes
# Open Targets - GO processes
OPENTARGETS_GO_BP_NODE_LABELS = "Biological Process"
OPENTARGETS_GO_MF_NODE_LABELS = "Molecular Function"
OPENTARGETS_GO_CC_NODE_LABELS = "Cellular Component"
OPENTARGETS_GO_NODE_MAIN_LABEL = "go_id"
OPENTARGETS_GO_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "labels": None,
}
OPENTARGETS_GO_EDGE_LABEL = "part_of"
OPENTARGETS_GO_EDGE_ATTRS = {"source": OPENTARGETS, "label": OPENTARGETS_GO_EDGE_LABEL}

# Compound node
# Open Targets - Compound
OPENTARGETS_COMPOUND_NODE_LABELS = "Compound"
OPENTARGETS_COMPOUND_NODE_MAIN_LABEL = "compound_cid"
OPENTARGETS_COMPOUND_NODE_ATTRS = {
    "source": OPENTARGETS,
    "name": None,
    "id": None,
    "chembl_id": None,
    "drugbank_id": None,
    "compound_cid": None,
    "clincal_trial_phase": None,
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
    "labels": OPENTARGETS_SIDE_EFFECT_NODE_LABELS,
}
OPENTARGETS_SIDE_EFFECT_EDGE_LABEL = "has_side_effect"

# MolMeDB - Gene/Protein input
MOLMEDB_PROTEIN_COMPOUND_NODE_LABELS = "Compound"
MOLMEDB_PROTEIN_COMPOUND_NODE_MAIN_LABEL = "compound_cid"
MOLMEDB_PROTEIN_COMPOUND_NODE_ATTRS = {
    "source": MOLMEDB,
    "name": None,
    "id": None,
    "MolMeDB_id": None,
    "ChEBI_id": None,
    "drugbank_id": None,
    "compound_cid": None,
    "InChIKey": None,
    "SMILES": None,
    "source_pmid": None,
    "labels": MOLMEDB_PROTEIN_COMPOUND_NODE_LABELS,
}
MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL = "inhibits"
MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS = {
    "source": MOLMEDB,
    "label": MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
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

# Gene node
GENE_NODE_LABELS = "Gene"
# STRING
# TODO: to be checked
STRING_EDGE_MAIN_LABEL = "stringdb_link_to"
STRING_EDGE_LABEL = "StringDB_ppi"
STRING_EDGE_ATTRS = {
    "source": STRING,
    "score": None,
    "label": STRING_EDGE_LABEL,
}

# Wikidata

# TODO: to be checked


# Output from the explorers
