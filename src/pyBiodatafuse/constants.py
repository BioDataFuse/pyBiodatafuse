# coding: utf-8

"""Python constant file."""


"""
API endpoints for each data source
"""
BRIDGEDB_ENDPOINT = "https://webservice.bridgedb.org"
BGEE_ENDPOINT = "https://www.bgee.org/sparql/"
DISGENET_ENDPOINT = "https://api.disgenet.com/api/v1/gda/summary"
ENSEMBL_ENDPOINT = "https://rest.ensembl.org"
INTACT_ENDPOINT = "https://www.ebi.ac.uk/intact"
KEGG_ENDPOINT = "https://rest.kegg.jp"
MINERVA_ENDPOINT = "https://minerva-net.lcsb.uni.lu/api/"
MOLMEDB_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/molmedb"
OPENTARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
PUBCHEM_ENDPOINT = "https://idsm.elixir-czech.cz/sparql/endpoint/idsm"
STRING_ENDPOINT = "https://string-db.org/api"
WIKIDATA_ENDPOINT = "https://query-main.wikidata.org/sparql"
WIKIPATHWAYS_ENDPOINT = "https://sparql.wikipathways.org/sparql"
AOPWIKI_ENDPOINT = "https://aopwiki.rdf.bigcat-bioinformatics.org/sparql/"
NCBI_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov"  # part of StringDB
UNIPROT_ID_MAPPER_ENDPOINT = "https://rest.uniprot.org/idmapping/"

"""
All data sources
"""
BRIDGEDB = "BridgeDB"
BGEE = "Bgee"
DISGENET = "DISGENET"
ENSEMBL = "Ensembl"
INTACT = "IntAct"
KEGG = "KEGG"
MINERVA = "MINERVA"
MOLMEDB = "MolMeDB"
OPENTARGETS = "OpenTargets"
PUBCHEM = "PubChem"
STRING = "StringDB"
WIKIDATA = "Wikidata"
WIKIPATHWAYS = "WikiPathways"
AOPWIKIRDF = "AOP Wiki RDF"

"""
All weblinks for the datasources
"""
DATA_SOURCES = {
    BRIDGEDB: "https://www.bridgedb.org/",
    BGEE: "https://www.bgee.org/",
    DISGENET: "https://disgenet.com/",
    ENSEMBL: "https://www.ensembl.org/",
    INTACT: "https://www.ebi.ac.uk/intact/",
    KEGG: "https://www.genome.jp/kegg/pathway.html",
    MINERVA: "https://minerva.pages.uni.lu/doc/",
    MOLMEDB: "https://molmedb.upol.cz",
    OPENTARGETS: "https://www.opentargets.org/",
    PUBCHEM: "https://pubchem.ncbi.nlm.nih.gov/",
    STRING: "https://string-db.org/",
    WIKIDATA: "https://wikidata.org",
    WIKIPATHWAYS: "https://wikipathways.org",
    AOPWIKIRDF: "https://aopwiki.rdf.bigcat-bioinformatics.org",
}

"""
Annotator statistics constants
"""
QUERY = "query"
NUM_NODES = "number_of_added_nodes"
NUM_EDGES = "number_of_added_edges"

# DataFrame Columns

TARGET_COL = "target"
TARGET_SOURCE_COL = "target.source"
IDENTIFIER_COL = "identifier"
IDENTIFIER_SOURCE_COL = "identifier.source"

BGEE_GENE_EXPRESSION_LEVELS_COL = f"{BGEE}_gene_expression_levels"
DISGENET_DISEASE_COL = f"{DISGENET}_diseases"
ENSEMBL_HOMOLOG_COL = f"{ENSEMBL}_homologs"
INTACT_INTERACT_COL = f"{INTACT}_gene_interactions"
INTACT_COMPOUND_INTERACT_COL = f"{INTACT}_compound_interactions"
KEGG_PATHWAY_COL = f"{KEGG}_pathways"
LITERATURE_DISEASE_COL = "literature_based_info"
OPENTARGETS_REACTOME_COL = f"{OPENTARGETS}_reactome"
OPENTARGETS_GO_COL = f"{OPENTARGETS}_go"
OPENTARGETS_DISEASE_COL = f"{OPENTARGETS}_diseases"
OPENTARGETS_DISEASE_COMPOUND_COL = f"{OPENTARGETS}_disease_compounds"
OPENTARGETS_GENE_COMPOUND_COL = f"{OPENTARGETS}_gene_compounds"
MINERVA_PATHWAY_COL = f"{MINERVA}_pathways"
MOLMEDB_PROTEIN_COMPOUND_COL = f"{MOLMEDB}_transporter_inhibitor"
MOLMEDB_COMPOUND_PROTEIN_COL = f"{MOLMEDB}_transporter_inhibited"
PUBCHEM_COMPOUND_ASSAYS_COL = f"{PUBCHEM}_assays"
STRING_INTERACT_COL = f"{STRING}_ppi"
WIKIDATA_CC_COL = f"{WIKIDATA}_cellular_components"
WIKIPATHWAYS_MOLECULAR_COL = f"{WIKIPATHWAYS}_molecular"
AOPWIKI_GENE_COL = "aop_gene"  # todo fix this
AOPWIKI_COMPOUND_COL = "pubchem_compound"  # todo fix this

# Ontologies and vocabularies namespaces
PUBCHEM_COMPOUND = "PubChem Compound"
HPO = "HPO"
NCI = "NCI"
OMIM_IDS = "OMIM|MIM"
OMIM = "OMIM"
MONDO = "MONDO"
ORDO = "ORDO"
EFO = "EFO"
DO = "DO"
MESH = "MESH"
UMLS = "UMLS"
NCBI_GENE = "NCBI Gene"
CHEBI = "ChEBI"
UNIPROT = "Uniprot-TrEMBL"
CHEMBL = "chembl_id"
INCHIKEY = "InChIKey"
SMILES = "SMILES"
DRUGBANK = "DrugBank"

# Input type for each data source
BGEE_GENE_INPUT_ID = ENSEMBL
DISGENET_GENE_INPUT_ID = NCBI_GENE
ENSEMBL_GENE_INPUT_ID = ENSEMBL
INTACT_GENE_INPUT_ID = ENSEMBL
INTACT_COMPOUND_INPUT_ID = CHEBI
KEGG_GENE_INPUT_ID = NCBI_GENE
OPENTARGETS_GENE_INPUT_ID = ENSEMBL
OPENTARGETS_COMPOUND_INPUT_ID = PUBCHEM_COMPOUND
OPENTARGETS_COMPOUND_QUERY_INPUT_ID = CHEMBL
OPENTARGETS_DISEASE_INPUT_ID_1 = EFO
OPENTARGETS_DISEASE_INPUT_ID_2 = MONDO
MINERVA_GENE_INPUT_ID = ENSEMBL
MOLMEDB_PROTEIN_INPUT_ID = UNIPROT
MOLMEDB_COMPOUND_INPUT_ID = INCHIKEY
PUBCHEM_COMPOUND_INPUT_ID = UNIPROT
STRING_GENE_INPUT_ID = ENSEMBL
STRING_GENE_LINK_ID = f"{ENSEMBL}_link"
WIKIDATA_GENE_INPUT_ID = NCBI_GENE
WIKIPATHWAYS_GENE_INPUT_ID = NCBI_GENE
PATENT_INPUT_ID = PUBCHEM_COMPOUND
AOPWIKI_GENE_INPUT_ID = ENSEMBL
AOPWIKI_COMPOUND_INPUT_ID = PUBCHEM_COMPOUND

PATHWAY_ID = "pathway_id"
PATHWAY_LABEL = "pathway_label"
PATHWAY_GENE_COUNTS = "pathway_gene_counts"
PATHWAY_COMPOUNDS = "pathway_compounds"
PATHWAYS = "pathways"

SOURCE_PMID = "source_pmid"

ENTITY_SYMBOL = "symbol"
ENTITY_REFS = "references"
ENTITY_TYPE = "type"
ENTITY_NAME = "name"

# Output annotation for each data source
"""
BGEE variables
"""
ANATOMICAL_ID = "anatomical_entity_id"
ANATOMICAL_ENTITY_ID = "UBERON"
ANATOMICAL_NAME = "anatomical_entity_name"
EXPRESSION_LEVEL = "expression_level"
CONFIDENCE_ID = "confidence_level_id"
CONFIDENCE_LEVEL_ID = "CIO"
CONFIDENCE_LEVEL_NAME = "confidence_level_name"
DEVELOPMENTAL_ID = "developmental_stage_id"
DEVELOPMENTAL_STAGE_ID = "HsapDv|UBERON"
DEVELOPMENTAL_STAGE_NAME = "developmental_stage_name"

BGEE_GENE_EXPRESSION_OUTPUT_DICT = {
    ANATOMICAL_ID: str,
    ANATOMICAL_NAME: str,
    EXPRESSION_LEVEL: float,
    CONFIDENCE_ID: str,
    CONFIDENCE_LEVEL_NAME: str,
    DEVELOPMENTAL_ID: str,
    DEVELOPMENTAL_STAGE_NAME: str,
}

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

"""
DISGENET variables
"""
DISEASE_NAME = "disease_name"
DISEASE_TYPE = "disease_type"
DISEASE_UMLSCUI = "disease_umlscui"
DISGENET_SCORE = "score"
DISGENET_EI = "ei"
DISGENET_EL = "el"

DISGENET_DISEASE_OUTPUT_DICT = {
    DISEASE_NAME: str,
    HPO: str,  # "HPO_HP:0100013"
    NCI: str,  # "NCI_C2910"
    OMIM: str,  # "OMIM_607906"
    MONDO: str,  # "MONDO_0021100"
    ORDO: str,  # "ORDO_137"
    EFO: str,  # "EFO_0003756"
    DO: str,  # "DO_0060041"
    MESH: str,  # "MESH_D000067877"
    UMLS: str,  # "UMLS_C1510586"
    DISEASE_TYPE: str,
    DISGENET_SCORE: float,
    DISGENET_EI: float,
    DISGENET_EL: str,
}

VALUE_CHECK_LIST = [
    HPO,
    NCI,
    OMIM,
    MONDO,
    ORDO,
    EFO,
    DO,
    MESH,
    UMLS,
]

"""
INTACT variables
"""
INTACT_INTERACTION_ID = "interaction_id"
INTACT_INTERACTOR_ID_A = "interactor_id_A"
INTACT_INTERACTOR_ID_B = "interactor_id_B"
INTACT_SCORE = "score"
INTACT_BIOLOGICAL_ROLE_A = "biological_role_A"
INTACT_BIOLOGICAL_ROLE_B = "biological_role_B"
INTACT_TYPE = "type"
INTACT_DETECTION_METHOD = "detection_method"
INTACT_HOST_ORGANISM = "host_organism"
INTACT_INTERACTOR_A_NAME = "interactor_A_name"
INTACT_INTERACTOR_B_NAME = "interactor_B_name"
INTACT_INTERACTOR_A_SPECIES = "interactor_A_species"
INTACT_INTERACTOR_B_SPECIES = "interactor_B_species"
INTACT_MOLECULE_A = "molecule_A"
INTACT_MOLECULE_B = "molecule_B"
INTACT_ID_A = "id_A"
INTACT_ID_B = "id_B"
INTACT_PUBMED_PUBLICATION_ID = "pubmed_publication_id"

INTACT_COMPOUND_INTERACTION_TYPES = ["compound_compound", "compound_gene", "both_compounds"]
INTACT_GENE_INTERACTION_TYPES = ["gene_gene", "gene_compound", "both"]

INTACT_OUTPUT_DICT = {
    INTACT_INTERACTION_ID: str,
    INTACT_INTERACTOR_ID_A: str,
    INTACT_INTERACTOR_ID_B: str,
    INTACT_SCORE: float,
    INTACT_BIOLOGICAL_ROLE_A: str,
    INTACT_BIOLOGICAL_ROLE_B: str,
    INTACT_TYPE: str,
    INTACT_DETECTION_METHOD: str,
    INTACT_HOST_ORGANISM: str,
    INTACT_INTERACTOR_A_NAME: str,
    INTACT_INTERACTOR_B_NAME: str,
    INTACT_INTERACTOR_A_SPECIES: str,
    INTACT_INTERACTOR_B_SPECIES: str,
    INTACT_MOLECULE_A: str,
    INTACT_MOLECULE_B: str,
    INTACT_ID_A: str,
    INTACT_ID_B: str,
    INTACT_PUBMED_PUBLICATION_ID: str,
}

"""
KEGG variables
"""
KEGG_IDENTIFIER = f"{KEGG}_id"
KEGG_COMPOUND_NAME = f"{KEGG}_compound_name"
KEGG_COMPOUND_OUTPUT_DICT = {
    KEGG_IDENTIFIER: str,
    KEGG_COMPOUND_NAME: str,
}

KEGG_PATHWAY_OUTPUT_DICT = {
    KEGG_IDENTIFIER: str,
    PATHWAYS: {
        PATHWAY_ID: str,
        PATHWAY_LABEL: str,
        PATHWAY_GENE_COUNTS: int,
        PATHWAY_COMPOUNDS: list,
    },
}

"""
MINERVA variables
"""
VALID_MINERVA_ENTITIES = [
    "Compartment",
    "Complex",
    "Drug",
    "Gene",
    "Ion",
    "Phenotype",
    "Protein",
    "RNA",
    "Simple molecule",
    "Pathway",
]
INTERESTED_INFO = [ENTITY_TYPE, ENTITY_REFS, ENTITY_SYMBOL, ENTITY_NAME, ENSEMBL]

MINERVA_PATHWAY_OUTPUT_DICT = {
    PATHWAY_ID: str,
    PATHWAY_LABEL: str,
    PATHWAY_GENE_COUNTS: int,
}

MINERVA_PATHWAY_DEFAULT_ID = MINERVA

"""
MolMeDB variables
"""
MOLMEDB_COMPOUND_NAME = f"{MOLMEDB}_compound_name"
MOLMEDB_INCHIKEY = f"{MOLMEDB}_inchikey"
MOLMEDB_SMILES = f"{MOLMEDB}_smiles"
MOLMEDB_ID = f"{MOLMEDB}_id"
MOLMEDB_HGNC_SYMBOL = f"{MOLMEDB}_hgnc_symbol"

MAP_COMPOUND_COL_NAMES = {
    "transporterID": TARGET_COL,
    "label": MOLMEDB_COMPOUND_NAME,
    "InChIKey": MOLMEDB_INCHIKEY,
    "SMILES": MOLMEDB_SMILES,
    "molmedb_id": MOLMEDB_ID,
    "pubchem_compound_id": "compound_cid",
}

MOLMEDB_PROTEIN_COMPOUND_OUTPUT_DICT = {
    MOLMEDB_COMPOUND_NAME: str,
    MOLMEDB_INCHIKEY: str,
    MOLMEDB_SMILES: str,
    MOLMEDB_ID: str,
    SOURCE_PMID: str,
}
MOLMEDB_COMPOUND_DEFAULT_ID = "MM"

MOLMEDB_INHIBITOR_INCHIKEY = "inhibitorInChIKey"
MOLMEDB_UNIPROT_TREMBL_ID = "uniprot_trembl_id"
MOLMEDB_HGNC_ID = "hgcn_id"

MAP_GENE_COL_NAMES = {
    MOLMEDB_INHIBITOR_INCHIKEY: TARGET_COL,
    MOLMEDB_HGNC_ID: MOLMEDB_HGNC_SYMBOL,
}

# MolMeDB - Compound input
MOLMEDB_COMPOUND_PROTEIN_OUTPUT_DICT = {
    # MOLMEDB_UNIPROT_TREMBL_ID: str,
    MOLMEDB_HGNC_SYMBOL: str,
    SOURCE_PMID: str,
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


# WikiPathways
WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
    "pathway_gene_count": int,
}

WIKIPATHWAYS_MOLECULAR_GENE_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
    "targetGene": str,
    "targetProtein": str,
    "targetMetabolite": str,
    "mimtype": str,
    "rhea_id": str,
}

# Open Targets - Reactome
OPENTARGETS_REACTOME_OUTPUT_DICT = {
    "pathway_id": str,
    "pathway_label": str,
}
# PATHWAY_ID = "MINERVA|WP|R-"  # ID Start with WP or R-

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
STRING_OUTPUT_DICT = {
    "stringdb_link_to": str,
    STRING_GENE_INPUT_ID: str,
    "score": float,
    "Uniprot-TrEMBL": str,
    "Uniprot-TrEMBL_link": str,
}

# AOPWIKI
AOPWIKI_GENE_OUTPUT_DICT = {
    "aop": str,
    "aop_title": str,
    "MIEtitle": str,
    "MIE": str,
    "KE_downstream": str,
    "KE_downstream_title": str,
    "KER": str,
    "ao": str,
    "ao_title": str,
    "KE_upstream": str,
    "KE_upstream_title": str,
    "KE_upstream_organ": str,
    "KE_downstream_organ": str,
    "pubchem_compound": str,
}

AOPWIKI_COMPOUND_OUTPUT_DICT = {
    "aop": str,
    "aop_title": str,
    "MIEtitle": str,
    "MIE": str,
    "KE_downstream": str,
    "KE_downstream_title": str,
    "KER": str,
    "ao": str,
    "ao_title": str,
    "KE_upstream": str,
    "KE_upstream_title": str,
    "KE_upstream_organ": str,
    "KE_downstream_organ": str,
}


"""
Node and edge main lable and attributes
"""
# Common attributes
DATASOURCE = "datasource"
NAME = ENTITY_NAME
ID = "id"
LABEL = "label"
EDGE_HASH = "edge_hash"

# Node types
GENE_NODE_LABEL = "Gene"
DISEASE_NODE_LABEL = "Disease"
COMPOUND_NODE_LABEL = "Compound"
ANATOMICAL_NODE_LABEL = "Anatomical Entity"
PATHWAY_NODE_LABEL = "Pathway"
GO_BP_NODE_LABEL = "Biological Process"
GO_MF_NODE_LABEL = "Molecular Function"
GO_CC_NODE_LABEL = "Cellular Component"
SIDE_EFFECT_NODE_LABEL = "Side Effect"
HOMOLOG_NODE_LABEL = "Homolog"

"""
Anatomical entity nodes
"""
BGEE_ANATOMICAL_NODE_MAIN_LABEL = ANATOMICAL_ID
BGEE_ANATOMICAL_NODE_ATTRS = {
    DATASOURCE: BGEE,
    NAME: None,
    ID: None,
    LABEL: ANATOMICAL_NODE_LABEL,
}
BGEE_GENE_ANATOMICAL_EDGE_LABEL = "expressed_by"
BGEE_EDGE_ATTRS = {
    DATASOURCE: BGEE,
    EXPRESSION_LEVEL: None,
    DEVELOPMENTAL_STAGE_NAME: None,
    DEVELOPMENTAL_ID: None,
    CONFIDENCE_LEVEL_NAME: None,
    CONFIDENCE_ID: None,
    LABEL: BGEE_GENE_ANATOMICAL_EDGE_LABEL,
}

"""
Disease nodes
"""
DISEASE_NODE_MAIN_LABEL = UMLS
DISGENET_DISEASE_NODE_ATTRS = {
    DATASOURCE: DISGENET,
    NAME: None,
    ID: None,
    HPO: None,
    NCI: None,
    OMIM: None,
    MONDO: None,
    ORDO: None,
    EFO: None,
    DO: None,
    MESH: None,
    UMLS: None,
    DISEASE_TYPE: None,
    LABEL: DISEASE_NODE_LABEL,
}
GENE_DISEASE_EDGE_LABEL = "associated_with"
DISGENET_EDGE_ATTRS = {
    DATASOURCE: DISGENET,
    DISGENET_SCORE: None,
    DISGENET_EI: None,
    DISGENET_EL: None,
    LABEL: GENE_DISEASE_EDGE_LABEL,
}

"""
Pathway nodes
"""
GENE_COUNTS = "gene_counts"
PATHWAY_NODE_ATTRS = {
    DATASOURCE: None,
    NAME: None,
    ID: None,
    LABEL: PATHWAY_NODE_LABEL,
    # GENE_COUNTS: None,
}

GENE_PATHWAY_EDGE_LABEL = "part_of"
GENE_PATHWAY_EDGE_ATTRS = {
    DATASOURCE: None,
    LABEL: GENE_PATHWAY_EDGE_LABEL,
}

# IntAct interactions
SPECIES = "species"
INTACT_INTERACTION_TYPE = "interaction_type"
INTACT_PPI_EDGE_MAIN_LABEL = "intact_link_to"
INTACT_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
INTACT_NODE_ATTRS.update(
    {
        DATASOURCE: INTACT,
        ID: None,
        LABEL: None,
    }
)
INTACT_PPI_EDGE_ATTRS = {
    DATASOURCE: INTACT,
    INTACT_DETECTION_METHOD: None,
    INTACT_TYPE: None,
    LABEL: INTACT_PPI_EDGE_MAIN_LABEL,
}

KEGG_PATHWAY_NODE_MAIN_LABEL = PATHWAY_ID
KEGG_PATHWAY_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
KEGG_PATHWAY_NODE_ATTRS.update(
    {
        DATASOURCE: KEGG,
        ID: None,
        LABEL: PATHWAY_NODE_LABEL,
        GENE_COUNTS: None,
    }
)

MINERVA_PATHWAY_NODE_MAIN_LABEL = PATHWAY_ID
MINERVA_PATHWAY_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
MINERVA_PATHWAY_NODE_ATTRS.update(
    {
        DATASOURCE: MINERVA,
        ID: None,
        LABEL: PATHWAY_NODE_LABEL,
        GENE_COUNTS: None,
    }
)

# molecular pathway node
MOLECULAR_PATHWAY_NODE_MAIN_LABEL = "pathway_id"
MOLECULAR_GENE_NODE_ATTRS = {"datasource": WIKIPATHWAYS, "label": ""}
MOLECULAR_PATHWAY_NODE_ATTRS = {
    "pathway_id": "str",
    "pathway_label": "str",
    "id": "str",
    "labels": PATHWAY_NODE_LABEL,
}
MOLECULAR_GENE_PATHWAY_EDGE_LABEL = "part_of"
MOLECULAR_INTERACTION_EDGE_ATTRS = {"interaction_type": "str", "rhea_id": str}

GO_NODE_MAIN_LABEL = "go_id"
GO_NODE_ATTRS = {
    "datasource": OPENTARGETS,
    "name": None,
    "id": None,
    "labels": None,
}
GENE_GO_EDGE_LABEL = "part_of"
GENE_GO_EDGE_ATTRS = {"datasource": OPENTARGETS, "label": GENE_GO_EDGE_LABEL}

"""
Compound nodes
"""
COMPOUND_NODE_MAIN_LABEL = "compound_cid"

# IntAct compounds
INTACT_COMPOUND_NODE_MAIN_LABEL = "compounds"
MOLECULE = "molecule"
INTACT_COMPOUND_NODE_ATTRS = {
    DATASOURCE: INTACT,
    ID: None,
    NAME: None,
    LABEL: COMPOUND_NODE_LABEL,
}
INTACT_COMPOUND_EDGE_ATTRS = {
    DATASOURCE: INTACT,
    LABEL: None,
}

# KEGG compounds
KEGG_COMPOUND_NODE_MAIN_LABEL = "compounds"
KEGG_COMPOUND_NODE_ATTRS = {
    DATASOURCE: KEGG,
    ID: None,
    LABEL: COMPOUND_NODE_LABEL,
}

KEGG_COMPOUND_EDGE_LABEL = "contains"
KEGG_COMPOUND_EDGE_ATTRS = {
    DATASOURCE: KEGG,
    LABEL: None,
}

# MolMeDB
MOLMEDB_COMPOUND_NODE_ATTRS = {
    DATASOURCE: MOLMEDB,
    NAME: None,
    ID: None,
    MOLMEDB_INCHIKEY: None,
    MOLMEDB_SMILES: None,
    SOURCE_PMID: None,
    LABEL: COMPOUND_NODE_LABEL,
}
MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL = "inhibits"
MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS = {
    "datasource": MOLMEDB,
    "label": MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
}


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
    "labels": COMPOUND_NODE_LABEL,
}
OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS = {"datasource": OPENTARGETS, "label": None}
# Side effect
# Open Targets - Compound

SIDE_EFFECT_NODE_MAIN_LABEL = "adverse_effect"
# TODO: add the side effect id (adverse_effect id)
SIDE_EFFECT_NODE_ATTRS = {
    "datasource": OPENTARGETS,
    "name": None,
    "labels": SIDE_EFFECT_NODE_LABEL,
}
COMPOUND_SIDE_EFFECT_EDGE_LABEL = "has_side_effect"
COMPOUND_SIDE_EFFECT_EDGE_ATTRS = {
    "datasource": OPENTARGETS,
    "label": COMPOUND_SIDE_EFFECT_EDGE_LABEL,
}


# PubChem - Assays
PUBCHEM_COMPOUND_NODE_ATTRS = {
    "datasource": PUBCHEM,
    "name": None,
    "id": None,
    "inchi": None,
    "smiles": None,
    "labels": COMPOUND_NODE_LABEL,
}
PUBCHEM_GENE_COMPOUND_EDGE_ATTRS = {
    "datasource": PUBCHEM,
    "assay_type": None,
    "pubchem_assay_id": None,
    "outcome": None,
    "label": None,
}

"""
Gene node
"""


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

# Ensembl Homologs
ENSEMBL_HOMOLOG_NODE_ATTRS = {
    "datasource": ENSEMBL,
    "id": None,
    "labels": HOMOLOG_NODE_LABEL,
}
ENSEMBL_HOMOLOG_MAIN_LABEL = "homolog"
ENSEMBL_HOMOLOG_EDGE_LABEL = "is_homolog_of"
ENSEMBL_HOMOLOG_EDGE_ATTRS = {
    "datasource": ENSEMBL,
    "label": ENSEMBL_HOMOLOG_EDGE_LABEL,
}


# Literature
LITERATURE_NODE_MAIN_LABEL = "UMLS"
LITERATURE_DISEASE_NODE_ATTRS = {
    "datasource": None,
    "name": None,
    "id": None,
    "MONDO": None,
    "UMLS": None,
    "labels": DISEASE_NODE_LABEL,
}
LITERATURE_DISEASE_EDGE_ATTRS = {
    "datasource": None,
    "label": GENE_DISEASE_EDGE_LABEL,
}

"""
RDF Specific constants
"""

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
    "ppi": STRING_INTERACT_COL,
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
    "protein_node": f"{NAMESPACE_BINDINGS['obo']}NCIT_C17021",
    "ppi_node": f"{NAMESPACE_BINDINGS['obo']}NCIT_C18469",
    "compound_node": "http://semanticscience.org/resource/SIO_010038",
    "approved_compound": f"{NAMESPACE_BINDINGS['obo']}NCIT_C172573",
    "aid": f"{NAMESPACE_BINDINGS['obo']}CLO_0037244",
    "developmental_stage_node": f"{NAMESPACE_BINDINGS['obo']}NCIT_C43531",
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
    "chebi_inchi": f"{NAMESPACE_BINDINGS['obo']}chebi/inchi",
    "chebi_smiles": f"{NAMESPACE_BINDINGS['obo']}chebi/smiles",
    "cheminf_compound_id": "http://semanticscience.org/resource/CHEMINF_000140",
    "sio_is_part_of": f"{NAMESPACE_BINDINGS['sio']}SIO_000068",
    "has_gene_count": "",  # TODO gene count
    "precedes": "http://semanticscience.org/resource/SIO_000248",
    "is_preceded_by": "http://semanticscience.org/resource/SIO_000248",
    "sio_has_part": f"{NAMESPACE_BINDINGS['sio']}SIO_000028",
    "negatively_regulates": f"{NAMESPACE_BINDINGS['obo']}RO_0002449",
    "phase": f"{NAMESPACE_BINDINGS['obo']}PATO_0000083",
    "translation_of": f"{NAMESPACE_BINDINGS['obo']}so#translation_of",
    "translates_to": f"{NAMESPACE_BINDINGS['obo']}so#translates_to",
    "variant_of": f"{NAMESPACE_BINDINGS['obo']}so#variant_of",
}

# Classes for clinical phases

CLINICAL_PHASES = {
    "1.0": f"{NAMESPACE_BINDINGS['obo']}OPMI_0000368",
    "2.0": f"{NAMESPACE_BINDINGS['obo']}OPMI_0000369",
    "3.0": f"{NAMESPACE_BINDINGS['obo']}OPMI_0000370",
    "4.0": f"{NAMESPACE_BINDINGS['obo']}OPMI_0000371",
}

# GO Types
GO_TYPES = {
    "C": f"{NAMESPACE_BINDINGS['obo']}GO_0005575",
    "P": f"{NAMESPACE_BINDINGS['obo']}GO_0008150",
    "F": f"{NAMESPACE_BINDINGS['obo']}GO_0003674",
}

# Compound MoA

MOAS = {
    "ANTAGONIST": f"{NAMESPACE_BINDINGS['obo']}RO_0018029",
    "AGONIST": f"{NAMESPACE_BINDINGS['obo']}RO_0018027",
    "BLOCKER": f"{NAMESPACE_BINDINGS['obo']}RO_0003002",
    "INHIBITOR": f"{NAMESPACE_BINDINGS['obo']}RO_0012006",  # TODO this predicate has the wrong range
    "MODULATOR": f"{NAMESPACE_BINDINGS['obo']}RO_0011002",  # TODO needs more granularity
    "PARTIAL AGONIST": f"{NAMESPACE_BINDINGS['obo']}RO_0018027",  # TODO needs more granularity
    "INVERSE AGONIST": f"{NAMESPACE_BINDINGS['obo']}RO_0018028",
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
    f"{NAMESPACE_BINDINGS['obo']}": "obo",
    f"{NAMESPACE_BINDINGS['obo']}so#": "so",
    "https://biodatafuse.org/onto/bdf#": "bdfo",
    "https://minerva-net.lcsb.uni.lu/api/": "minerva",
    "https://reactome.org/content/detail/": "reactome",
    "https://www.uniprot.org/uniprotkb/": "uniprot",
    "http://identifiers.org/ensembl#": "ensembl",
}

AOPWIKI_LIFESTAGE_MAPPINGS = {
    "All life stages": "",
    "Adult": "",
    "Juvenile": "",
    "Development": "",
    "Birth to < 1 month": "",
    "Not Otherwise Specified": "",
    "Adults": "",
    "Adult, reproductively mature": "",
    "During development and at adulthood": "",
    "Old Age": "",
    "Embryo": "",
    "During brain development": "",
    "During brain development, adulthood and aging": "",
    "Foetal": "",
    "Larvae": "",
}

SOURCE_NAMESPACES = {
    "CAS": "https://commonchemistry.cas.org/detail?cas_rn=",
    "ChEBI": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:",
    "ChemSpider": "http://www.chemspider.com/Chemical-Structure.",
    "ChEMBL compound": "https://www.ebi.ac.uk/chembl/compound_report_card/",
    "DrugBank": "https://go.drugbank.com/drugs/",
    "HMDB": "https://hmdb.ca/metabolites/",
    "InChIKey": "https://www.chemspider.com/Search.aspx?q=",
    "KEGG Compound": "https://www.genome.jp/dbget-bin/www_bget?cpd:",
    "KEGG Drug": "https://www.genome.jp/dbget-bin/www_bget?dr:",
    "KEGG Glycan": "https://www.genome.jp/dbget-bin/www_bget?gl:",
    "LIPID MAPS": "https://www.lipidmaps.org/data/LMSDRecord.php?LMID=",
    "LipidBank": "http://lipidbank.jp/cgi-bin/detail.cgi?id=",
    "PharmGKB Drug": "https://www.pharmgkb.org/chemical/",
    "PubChem Compound": "https://pubchem.ncbi.nlm.nih.gov/compound/",
    "PubChem Substance": "https://pubchem.ncbi.nlm.nih.gov/substance/",
    "SwissLipids": "https://www.swisslipids.org/#/entity/",
    "TTD Drug": "http://db.idrblab.net/ttd/data/drug/details/",
    "Wikidata": "https://www.wikidata.org/wiki/",
    "Wikipedia": "https://en.wikipedia.org/wiki/",
    # TODO ADD ALL
}
