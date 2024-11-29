# coding: utf-8

"""Python constant file."""


# Endpoints / API
BRIDGEDB_ENDPOINT = "https://webservice.bridgedb.org"
BGEE_ENDPOINT = "https://www.bgee.org/sparql/"
DISGENET_ENDPOINT = "https://api.disgenet.com/api/v1/gda/summary"
MINERVA_ENDPOINT = "https://minerva-net.lcsb.uni.lu/api/"
MOLMEDB_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/molmedb"
NCBI_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov"
OPENTARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
PUBCHEM_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/idsm"
STRING_ENDPOINT = "https://string-db.org/api"
UNIPROT_ENDPOINT = "https://rest.uniprot.org"
WIKIDATA_ENDPOINT = "https://query-main.wikidata.org/sparql"
WIKIPATHWAYS_ENDPOINT = "https://sparql.wikipathways.org/sparql"

# Data sources
BRIDGEDB = "BridgeDB"
BGEE = "Bgee"
DISGENET = "DISGENET"
MINERVA = "MINERVA"
MOLMEDB = "MolMeDB"
OPENTARGETS = "OpenTargets"
PUBCHEM = "PubChem"
STRING = "StringDB"
WIKIDATA = "Wikidata"
WIKIPATHWAYS = "WikiPathways"
OPENTARGETS_REACTOME = "OpenTargets_reactome"

# DataFrame Columns

TARGET_SOURCE_COL = "target.source"
TARGET_COL = "target"
IDENTIFIER_COL = "identifier"
IDENTIFIER_SOURCE_COL = "identifier.source"
TARGET_SOURCE_COL = "target.source"
BGEE_GENE_EXPRESSION_LEVELS_COL = f"{BGEE}_gene_expression_levels"
DISGENET_DISEASE_COL = f"{DISGENET}_diseases"
OPENTARGETS_DISEASE_COL = f"{OPENTARGETS}_diseases"
LITERATURE_DISEASE_COL = "literature_based_info"
OPENTARGETS_REACTOME_COL = f"{OPENTARGETS}_reactome"
OPENTARGETS_GO_COL = f"{OPENTARGETS}_go"
OPENTARGETS_DISEASE_COMPOUND_COL = f"{OPENTARGETS}_disease_compounds"
OPENTARGETS_GENE_COMPOUND_COL = f"{OPENTARGETS}_gene_compounds"
MOLMEDB_PROTEIN_COMPOUND_COL = f"{MOLMEDB}_transporter_inhibitor"
MOLMEDB_COMPOUND_PROTEIN_COL = f"{MOLMEDB}_transporter_inhibited"
PUBCHEM_COMPOUND_ASSAYS_COL = f"{PUBCHEM}_assays"
STRING_PPI_COL = f"{STRING}_ppi"
WIKIDATA_CC_COL = f"{WIKIDATA}_cellular_components"

# Input type for each data source
BGEE_GENE_INPUT_ID = "Ensembl"
DISGENET_GENE_INPUT_ID = "NCBI Gene"
MINERVA_GENE_INPUT_ID = "Ensembl"
MOLMEDB_PROTEIN_INPUT_ID = "Uniprot-TrEMBL"
MOLMEDB_COMPOUND_INPUT_ID = "InChIKey"
OPENTARGETS_GENE_INPUT_ID = "Ensembl"
OPENTARGETS_COMPOUND_INPUT_ID = "PubChem Compound"  # If using bridgedb mapping
OPENTARGETS_COMPOUND_QUERY_INPUT_ID = "chembl_id"
OPENTARGETS_DISEASE_INPUT_ID_1 = "EFO"
OPENTARGETS_DISEASE_INPUT_ID_2 = "MONOD"
PUBCHEM_COMPOUND_INPUT_ID = "Uniprot-TrEMBL"
STRING_GENE_INPUT_ID = "Ensembl"
WIKIDATA_GENE_INPUT_ID = "NCBI Gene"
WIKIPATHWAYS_GENE_INPUT_ID = "NCBI Gene"
PATENT_INPUT_ID = "PubChem Compound"

# Output annotation for each data source
# Bgee
BGEE_GENE_EXPRESSION_OUTPUT_DICT = {
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

# DISGENET
DISGENET_DISEASE_OUTPUT_DICT = {
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
    "el": str,
}
HPO = "HPO"
NCI = "NCI"
OMIM = "OMIM|MIM"
MONDO = "MONDO"
ORDO = "ORDO"
EFO = "EFO"
DO = "DO"
MESH = "MESH"
UMLS = "UMLS"
LITERATURE_DISEASE_OUTPUT_DICT = {
    "disease_name": str,
    "UMLS": str,
    "MONDO": str,
    "source": str,
}

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

# MINERVA
MINERVA_PATHWAY_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
    "pathway_gene_count": int,
}

# WikiPathways
WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
    "pathway_gene_count": int,
}

# Open Targets - Reactome
OPENTARGETS_REACTOME_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
}
PATHWAY_ID = "MINERVA|WP|R-"  # ID Start with WP or R-

# Open Targets - GO processes
OPENTARGETS_GO_OUTPUT_DICT = {"go_id": str, "go_name": str, "go_type": str}
GO_ID = "GO"

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
    "adverse_effect": list,
    # "mechanisms_of_action": list,
}
CHEMBL_ID = "CHEMBL"
DRUGBANK_ID = "DrugBank"
RELATION = "inhibits|activates"
OPENTARGETS_COMPOUND_DISEASE_RELATION = "treats"


# MolMeDB - Gene/Protein input
MOLMEDB_PROTEIN_COMPOUND_OUTPUT_DICT = {
    "compound_name": str,
    "inchikey": str,
    "smiles": str,
    "compound_cid": str,
    "molmedb_id": str,
    "source_pmid": str,
    "chebi_id": str,
    "drugbank_id": str,
    "uniprot_trembl_id": str,  # uniprot id of isoform
}
MOLMEDB_ID = "MM"
DRUGBANK_ID = "DrugBank"

# MolMeDB - Compound input
MOLMEDB_COMPOUND_PROTEIN_OUTPUT_DICT = {
    "uniprot_trembl_id": str,
    "hgnc_symbol": str,
    "source_pmid": str,
}
UNIPROT_TREMBL_ID = "P"

# PubChem - Assays
PUBCHEM_COMPOUND_OUTPUT_DICT = {
    "pubchem_assay_id": str,
    "assay_type": str,
    "outcome": str,
    "compound_cid": str,
    "compound_name": str,
    "smiles": str,
    "inchi": str,
}
OUTCOME = "active|inactive"
INCHI = "InChI"

# STRING
STRING_OUTPUT_DICT = {"stringdb_link_to": str, STRING_GENE_INPUT_ID: str, "score": float}


""" Node and edge main lable and attributes for each data source """
# Anatomical entity node
# Bgee
ANATOMICAL_NODE_LABELS = "Anatomical Entity"
BGEE_ANATOMICAL_NODE_MAIN_LABEL = "anatomical_entity_id"
BGEE_ANATOMICAL_NODE_ATTRS = {
    "datasource": BGEE,
    "name": None,
    "id": None,
    "labels": ANATOMICAL_NODE_LABELS,
}
BGEE_GENE_ANATOMICAL_EDGE_LABEL = "expressed_by"
BGEE_EDGE_ATTRS = {
    "datasource": BGEE,
    "expression_level": None,
    "developmental_stage_name": None,
    "developmental_stage_id": None,
    "confidence_level_name": None,
    "confidence_level_id": None,
    "label": BGEE_GENE_ANATOMICAL_EDGE_LABEL,
}

# Disease node
# DISGENET
DISEASE_NODE_LABELS = "Disease"
DISEASE_NODE_MAIN_LABEL = "UMLS"
DISGENET_DISEASE_NODE_ATTRS = {
    "datasource": DISGENET,
    "name": None,
    "id": None,
    "HPO": None,
    "NCI": None,
    "OMIM": None,
    "MONDO": None,
    "ORDO": None,
    "EFO": None,
    "DO": None,
    "MESH": None,
    "UMLS": None,
    "disease_type": None,
    "disease_umlscui": None,
    "labels": DISEASE_NODE_LABELS,
}
GENE_DISEASE_EDGE_LABEL = "associated_with"
DISGENET_EDGE_ATTRS = {
    "datasource": DISGENET,
    "score": None,
    "ei": None,
    "el": None,
    "label": GENE_DISEASE_EDGE_LABEL,
}

# Literature
LITERATURE_NODE_MAIN_LABEL = "UMLS"
LITERATURE_DISEASE_NODE_ATTRS = {
    "datasource": None,
    "name": None,
    "id": None,
    "MONDO": None,
    "UMLS": None,
    "labels": DISEASE_NODE_LABELS,
}
LITERATURE_DISEASE_EDGE_ATTRS = {
    "datasource": None,
    "label": GENE_DISEASE_EDGE_LABEL,
}

# TODO: The disease annotations are not curated and will be used again when the OpenTarget annotation improves.
# Open Targets - Disease
# OPENTARGETS_DISEASE_NODE_ATTRS = {
#     "datasource": OPENTARGETS,
#     "name": None,
#     "id": None,
#     "therapeutic_areas": None,
#     "labels": DISEASE_NODE_LABELS,
# }
# OPENTARGETS_DISEASE_EDGE_ATTRS = {
#     "datasource": OPENTARGETS,
#     "label": GENE_DISEASE_EDGE_LABEL,
# }


# Pathway node
# MINERVA, WikiPathways, Open Targets - Reactome
PATHWAY_NODE_LABELS = "Pathway"
PATHWAY_NODE_MAIN_LABEL = "pathway_id"
PATHWAY_NODE_ATTRS = {
    "datasource": None,
    "name": None,
    "id": None,
    "gene_count": None,
    "labels": PATHWAY_NODE_LABELS,
}  # TODO: Yojana, would it be possible to add pathway size here (gene_count)
GENE_PATHWAY_EDGE_LABEL = "part_of"
GENE_PATHWAY_EDGE_ATTRS = {"datasource": None, "label": GENE_PATHWAY_EDGE_LABEL}

# GO nodes
# Open Targets - GO processes
GO_BP_NODE_LABELS = "Biological Process"
GO_MF_NODE_LABELS = "Molecular Function"
GO_CC_NODE_LABELS = "Cellular Component"
GO_NODE_MAIN_LABEL = "go_id"
GO_NODE_ATTRS = {
    "datasource": OPENTARGETS,
    "name": None,
    "id": None,
    "labels": None,
}
GENE_GO_EDGE_LABEL = "part_of"
GENE_GO_EDGE_ATTRS = {"datasource": OPENTARGETS, "label": GENE_GO_EDGE_LABEL}

# Compound node
# Open Targets - Compound
COMPOUND_NODE_LABELS = "Compound"
COMPOUND_NODE_MAIN_LABEL = "compound_cid"
OPENTARGETS_COMPOUND_NODE_ATTRS = {
    "datasource": OPENTARGETS,
    "name": None,
    "id": None,
    "chembl_id": None,
    "drugbank_id": None,
    "compound_cid": None,
    "clincal_trial_phase": None,
    "is_approved": None,
    "adverse_effect_count": None,
    "labels": COMPOUND_NODE_LABELS,
}
OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS = {"datasource": OPENTARGETS, "label": None}
# Side effect
# Open Targets - Compound
SIDE_EFFECT_NODE_LABELS = "Side Effect"
SIDE_EFFECT_NODE_MAIN_LABEL = "adverse_effect"
# TODO: add the side effect id (adverse_effect id)
SIDE_EFFECT_NODE_ATTRS = {
    "datasource": OPENTARGETS,
    "name": None,
    "labels": SIDE_EFFECT_NODE_LABELS,
}
COMPOUND_SIDE_EFFECT_EDGE_LABEL = "has_side_effect"
COMPOUND_SIDE_EFFECT_EDGE_ATTRS = {
    "datasource": OPENTARGETS,
    "label": COMPOUND_SIDE_EFFECT_EDGE_LABEL,
}

# MolMeDB - Gene/Protein input
MOLMEDB_COMPOUND_NODE_ATTRS = {
    "datasource": MOLMEDB,
    "name": None,
    "id": None,
    "molmedb_id": None,
    "chebi_id": None,
    "drugbank_id": None,
    "compound_cid": None,
    "inchikey": None,
    "smiles": None,
    "source_pmid": None,
    "labels": COMPOUND_NODE_LABELS,
}
MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL = "inhibits"
MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS = {
    "datasource": MOLMEDB,
    "label": MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
}

# PubChem - Assays
PUBCHEM_COMPOUND_NODE_ATTRS = {
    "datasource": PUBCHEM,
    "name": None,
    "id": None,
    "inchi": None,
    "smiles": None,
    "labels": COMPOUND_NODE_LABELS,
}
PUBCHEM_GENE_COMPOUND_EDGE_ATTRS = {
    "datasource": PUBCHEM,
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
STRING_PPI_EDGE_MAIN_LABEL = "stringdb_link_to"
STRING_PPI_EDGE_LABEL = "interacts_with"
STRING_PPI_EDGE_ATTRS = {
    "datasource": STRING,
    "score": None,
    "label": STRING_PPI_EDGE_LABEL,
}

# Disease - Compound edge
OPENTARGETS_DISEASE_COMPOUND_EDGE_ATTRS = {
    "datasource": OPENTARGETS,
    "label": None,
}

# Wikidata

# TODO: to be checked

# Mapper from namespace to BridgeDB datasource
COMPOUND_NAMESPACE_MAPPER = {"pubchem.compound": "PubChem Compound", "CHEMBL": "ChEMBL compound"}

# RDF (rdflib constants and namespaces)


DATA_TYPES_SOURCES = {
    "identifier": IDENTIFIER_COL,
    "identifier_source": IDENTIFIER_SOURCE_COL,
    "target": TARGET_COL,
    "target_source": TARGET_SOURCE_COL,
    "gene_expression": BGEE_GENE_EXPRESSION_LEVELS_COL,
    "assays": PUBCHEM_COMPOUND_ASSAYS_COL,
    "go": OPENTARGETS_GO_COL,
    "gene_compounds": OPENTARGETS_GENE_COMPOUND_COL,
    "literature_disease": LITERATURE_DISEASE_COL,
    "transporter_inhibitor": MOLMEDB_PROTEIN_COMPOUND_COL,
    "ppi": STRING_PPI_COL,
    "disgenet": DISGENET_DISEASE_COL,
    "opentargets_disease": OPENTARGETS_DISEASE_COL,
}

# Dictionary to store namespace strings
NAMESPACE_BINDINGS = {
    "sio": "http://semanticscience.org/resource/",
    "hgnc": "http://bio2rdf.org/hgnc:",
    "obo": "http://purl.obolibrary.org/obo/",
    "ensembl": "https://identifiers.org/ensembl:",
    "dcat": "http://www.w3.org/ns/dcat#",
    "bdf": "https://biodatafuse.org/",
    "foaf": "http://xmlns.com/foaf/0.1/",
    "skos": "http://www.w3.org/2004/02/skos/core#",
    "owl": "http://www.w3.org/2002/07/owl#",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    "xsd": "http://www.w3.org/2001/XMLSchema#",
    "bdfo": "https://biodatafuse.org/onto/bdf#",
    "mondo": "https://monarchinitiative.org/disease/",
    "umls": "https://www.ncbi.nlm.nih.gov/medgen/",
    "so": "http://purl.obolibrary.org/obo/so#",
}

# Patterns URIs for nodes (one for each node in the schema)
URIS = {
    "gene_disease_association": "gene_disease_association",
    "gene_base_node": "gene",
    "gene_symbol_base_node": "gene_symbol",
    "source_base_node": "source",
    "data_source_base_node": "datasource",
    "score_base_node": "score",
    "experimental_process_node": "experimental_process",
    "anatomical_entity_base_node": "anatomical_entity",
    "life_cycle_base_node": "life_cycle",
    "gene_expression_value_base_node": "gene_expression_value",
}

# NODE TYPES
NODE_TYPES = {
    "gene_node": f"{NAMESPACE_BINDINGS['obo']}NCIT_C16612",
    "disease_node": f"{NAMESPACE_BINDINGS['obo']}NCIT_C7057",
    "gene_disease_association": f"{NAMESPACE_BINDINGS['sio']}SIO_000983",
    "score_node": f"{NAMESPACE_BINDINGS['obo']}NCIT_C25338",
    "data_source_node": "http://semanticscience.org/resource/SIO_000750",
    "gene_expression_value_node": f"{NAMESPACE_BINDINGS['sio']}SIO_001077",
    "anatomical_entity_node": "http://semanticscience.org/resource/SIO_001262",
    "tested_compound_node": "http://semanticscience.org/resource/SIO_010038",
    "source_database": f"{NAMESPACE_BINDINGS['sio']}SIO_000750",
    "experimental_process_node": "http://www.ebi.ac.uk/efo/EFO_0002694",
    "pathway_node": f"{NAMESPACE_BINDINGS['obo']}PW_0000001",
    "adverse_event_node": f"{NAMESPACE_BINDINGS['obo']}OAE_0000001",
    "ensemble": "http://identifiers.org/ensembl/",
    "ncbi_disease": "https://www.ncbi.nlm.nih.gov/medgen/",
    "article": f"{NAMESPACE_BINDINGS['obo']}IAO_0000013",
    "protein_node": "http://purl.obolibrary.org/obo/NCIT_C17021",
    "ppi_node": "http://purl.obolibrary.org/obo/NCIT_C18469",
    "compound_node": "http://semanticscience.org/resource/SIO_010038",
    "approved_compound": "http://purl.obolibrary.org/obo/NCIT_C172573",
    "aid": "http://purl.obolibrary.org/obo/CLO_0037244",
    "developmental_stage_node": "http://purl.obolibrary.org/obo/NCIT_C43531",
    "el_node": "https://biodatafuse.org/onto/bdf#DisGeNET_Evidence_Level",
    "ei_node": "https://biodatafuse.org/onto/bdf#DisGeNET_Evidence_Index",
}

# PREDICATES
PREDICATES = {
    "sio_refers_to": f"{NAMESPACE_BINDINGS['sio']}SIO_000628",
    "sio_has_measurement_value": f"{NAMESPACE_BINDINGS['sio']}SIO_000216",
    "sio_has_source": f"{NAMESPACE_BINDINGS['sio']}SIO_000253",
    "sio_is_associated_with": f"{NAMESPACE_BINDINGS['sio']}SIO_001403",
    "sio_has_value": f"{NAMESPACE_BINDINGS['sio']}SIO_000300",
    "sio_has_input": f"{NAMESPACE_BINDINGS['sio']}SIO_000230",
    "sio_has_output": f"{NAMESPACE_BINDINGS['sio']}SIO_000229",
    "chebi_inchi": "http://purl.obolibrary.org/obo/chebi/inchi",
    "chebi_smiles": "http://purl.obolibrary.org/obo/chebi/smiles",
    "cheminf_compound_id": "http://semanticscience.org/resource/CHEMINF_000140",
    "sio_is_part_of": f"{NAMESPACE_BINDINGS['sio']}SIO_000068",
    "has_gene_count": "",  # TODO gene count
    "precedes": "http://semanticscience.org/resource/SIO_000248",
    "is_preceded_by": "http://semanticscience.org/resource/SIO_000248",
    "sio_has_part": f"{NAMESPACE_BINDINGS['sio']}SIO_000028",
    "negatively_regulates": f"{NAMESPACE_BINDINGS['obo']}RO_0002449",
    "phase": "http://purl.obolibrary.org/obo/PATO_0000083",
    "translation_of": "http://purl.obolibrary.org/obo/so#translation_of",
    "translates_to": "http://purl.obolibrary.org/obo/so#translates_to",
    "variant_of": "http://purl.obolibrary.org/obo/so#variant_of",
}

# Classes for clinical phases

CLINICAL_PHASES = {
    "1.0": "http://purl.obolibrary.org/obo/OPMI_0000368",
    "2.0": "http://purl.obolibrary.org/obo/OPMI_0000369",
    "3.0": "http://purl.obolibrary.org/obo/OPMI_0000370",
    "4.0": "http://purl.obolibrary.org/obo/OPMI_0000371",
}

# GO Types
GO_TYPES = {
    "C": "http://purl.obolibrary.org/obo/GO_0005575",
    "P": "http://purl.obolibrary.org/obo/GO_0008150",
    "F": "http://purl.obolibrary.org/obo/GO_0003674",
}

# Compound MoA

MOAS = {
    "ANTAGONIST": "http://purl.obolibrary.org/obo/RO_0018029",
    "AGONIST": "http://purl.obolibrary.org/obo/RO_0018027",
    "BLOCKER": "http://purl.obolibrary.org/obo/RO_0003002",
    "INHIBITOR": "http://purl.obolibrary.org/obo/RO_0012006",  # TODO this predicate has the wrong range
    "MODULATOR": "http://purl.obolibrary.org/obo/RO_0011002",  # TODO needs more granularity
    "PARTIAL AGONIST": "http://purl.obolibrary.org/obo/RO_0018027",  # TODO needs more granularity
    "INVERSE AGONIST": "http://purl.obolibrary.org/obo/RO_0018028",
}

# Data sources
DATA_SOURCES = {
    DISGENET: "https://disgenet.com/",
    WIKIPATHWAYS: "https://wikipathways.org",
    MINERVA: "https://minerva.pages.uni.lu/doc/",
    BRIDGEDB: "https://www.bridgedb.org/",
    STRING: "https://string-db.org/",
    OPENTARGETS: "https://www.opentargets.org/",
    BGEE: "https://www.bgee.org/",
    MOLMEDB: "https://molmedb.upol.cz",
    PUBCHEM: "https://pubchem.ncbi.nlm.nih.gov/",
    WIKIDATA: "https://wikidata.org",
    OPENTARGETS_REACTOME: "https://www.opentargets.org/",
}

DISEASE_IDENTIFIER_TYPES = [
    "HPO",
    "NCI",
    "OMIM",
    "MONDO",
    "ORDO",
    "EFO",
    "DO",
    "MESH",
    "UMLS",
]

BASE_URLS_DBS = {
    "uniprot": "https://www.uniprot.org/uniprotkb/",
    "ensembl": "http://identifiers.org/ensembl:",
    "stringdb": "https://string-db.org/network/",
}

NAMESPACE_SHAPES = {
    "http://www.w3.org/1999/02/22-rdf-syntax-ns#": "rdf",
    "http://example.org/": "ex",
    "http://weso.es/shapes/": ":",
    "http://www.w3.org/2001/XMLSchema#": "xsd",
    "http://www.w3.org/2002/07/owl#": "owl",
    "http://purl.obolibrary.org/obo/": "obo",
    "http://purl.obolibrary.org/obo/so#": "so",
    "https://biodatafuse.org/onto/bdf#": "bdfo",
    "https://minerva-net.lcsb.uni.lu/api/": "minerva",
    "https://reactome.org/content/detail/": "reactome",
    "https://www.uniprot.org/uniprotkb/": "uniprot",
    "http://identifiers.org/ensembl#": "ensembl",
}
