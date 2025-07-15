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
TFLINK_DOWNLOAD_URL = "https://cdn.netbiol.org/tflink/download_files"
GPROFILER_VERSION_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/util/data_versions"
MITOCARTA_DOWNLOAD_URL = "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/"
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
OPENTARGETS_REACTOME = "OpenTargets_reactome"
TFLINK = "TFLink"
GPROFILER = "g:Profiler"
MITOCARTA = "MitoCarta"
AOPWIKIRDF = "AOP_Wiki_RDF"

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
METADATA = "metadata"

# DataFrame Columns

TARGET_COL = "target"
TARGET_SOURCE_COL = "target.source"
IDENTIFIER_COL = "identifier"
IDENTIFIER_SOURCE_COL = "identifier.source"
SOURCE_COL = "source"

BGEE_GENE_EXPRESSION_LEVELS_COL = f"{BGEE}_gene_expression_levels"
DISGENET_DISEASE_COL = f"{DISGENET}_diseases"
ENSEMBL_HOMOLOG_COL = f"{ENSEMBL}_homologs"
INTACT_INTERACT_COL = f"{INTACT}_gene_interactions"
INTACT_COMPOUND_INTERACT_COL = f"{INTACT}_compound_interactions"
KEGG_PATHWAY_COL = f"{KEGG}_pathways"
KEGG_COMPOUND_COL = f"{KEGG}_compounds"
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
STRING_INTERACT_COL = f"{STRING}_interactions"
WIKIDATA_CC_COL = f"{WIKIDATA}_cellular_components"
WIKIPATHWAYS_MOLECULAR_COL = f"{WIKIPATHWAYS}_molecular"
WIKIPATHWAYS_PATHWAY_COL = f"{WIKIPATHWAYS}_pathway"
MITOCART_PATHWAY_COL = f"{MITOCARTA}_pathways"
AOPWIKI_GENE_COL = f"{AOPWIKIRDF}_genes"
AOPWIKI_COMPOUND_COL = f"{AOPWIKIRDF}_compounds"

# Ontologies and vocabularies namespaces
PUBCHEM_COMPOUND = "PubChem Compound"
PUBCHEM_COMPOUND_CID = "CID"
HPO = "HPO"
NCI = "NCI"
OMIM = "OMIM"
MONDO = "MONDO"
ORDO = "ORDO"
EFO = "EFO"
DO = "DO"
MESH = "MESH"
UMLS = "UMLS"
NCBI_GENE = "NCBI Gene"
NCBI_GENE_ID = "NCBIGene"
CHEBI = "ChEBI"
UNIPROT_TREMBL = "Uniprot-TrEMBL"
CHEMBL = "CHEMBL"
CHEMBL_ID = "chembl_id"
KEGG_COMPOUND = "KEGG Compound"
INCHI = "InChI"
INCHIKEY = "InChIKey"
SMILES = "SMILES"
DRUGBANK = "DrugBank"
DRUGBANK_ID = "drugbank_id"
GO = "GO"
REACTOME = "Reactome"
WP = "WikiPathways"
BIODATAFUSE = "Biodatafuse"
UBERON = "UBERON"
CIO = "CIO"
WIKIPATHWAY = "WP"
AOP_PATHWAY = "AOP"
MOL_INITIATING_EVENT = "MIE"
KEY_EVENT = "KE"
ADVERSE_OUTCOME = "AO"

# Input type for each data source
BGEE_GENE_INPUT_ID = ENSEMBL
DISGENET_GENE_INPUT_ID = NCBI_GENE
ENSEMBL_GENE_INPUT_ID = ENSEMBL
INTACT_GENE_INPUT_ID = ENSEMBL
INTACT_COMPOUND_INPUT_ID = CHEBI
KEGG_GENE_INPUT_ID = NCBI_GENE
KEGG_COMPOUND_INPUT_ID = KEGG_COMPOUND
OPENTARGETS_GENE_INPUT_ID = ENSEMBL
OPENTARGETS_COMPOUND_INPUT_ID = PUBCHEM_COMPOUND
OPENTARGETS_COMPOUND_QUERY_INPUT_ID = CHEMBL_ID
OPENTARGETS_DISEASE_INPUT_ID = EFO
MINERVA_GENE_INPUT_ID = ENSEMBL
MOLMEDB_PROTEIN_INPUT_ID = UNIPROT_TREMBL
MOLMEDB_COMPOUND_INPUT_ID = INCHIKEY
PUBCHEM_COMPOUND_INPUT_ID = UNIPROT_TREMBL
STRING_GENE_INPUT_ID = ENSEMBL
STRING_GENE_LINK_ID = f"{ENSEMBL}_link"
WIKIDATA_GENE_INPUT_ID = NCBI_GENE
WIKIPATHWAYS_GENE_INPUT_ID = NCBI_GENE
PATENT_INPUT_ID = PUBCHEM_COMPOUND
TFLINK_GENE_INPUT_ID = NCBI_GENE
GPROFILER_GENE_INPUT_ID = ENSEMBL
MITOCARTA_GENE_INPUT_ID = ENSEMBL
AOPWIKI_GENE_INPUT_ID = ENSEMBL
AOPWIKI_COMPOUND_INPUT_ID = PUBCHEM_COMPOUND

PATHWAY_ID = "pathway_id"
PATHWAY_LABEL = "pathway_label"
PATHWAY_GENE_COUNTS = "pathway_gene_counts"
PATHWAY_GENES = "pathway_genes"
PATHWAY_COMPOUND_COUNTS = "pathway_compound_counts"
PATHWAY_COMPOUNDS = "pathway_compounds"
PATHWAYS = "pathways"

GO_ID = "go_id"
GO_NAME = "go_name"
GO_TYPE = "go_type"

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
ANATOMICAL_NAME = "anatomical_entity_name"
EXPRESSION_LEVEL = "expression_level"
CONFIDENCE_ID = "confidence_level_id"
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

BGEE_VALUE_CHECK_LIST = [UBERON, CIO, DEVELOPMENTAL_STAGE_ID]

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
DISGENET_SCORE = "score"

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
    f"{OMIM}|MIM",
    MONDO,
    ORDO,
    EFO,
    DO,
    MESH,
    UMLS,
]

"""
GProfiler variables
"""
ORGANISM = "organism"
GPROFILER_ID = "id"
GPROFILER_NAME = "name"
GPROFILER_INTERSECTIONS = "intersections"
GPROFILER_RESULT_COL = "gprofiler"
GPROFILER_COLS_TO_KEEP = ["source", GPROFILER_ID, GPROFILER_INTERSECTIONS, GPROFILER_RESULT_COL]

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
Mitocart variables
"""
MITO_SELECTED_COLUMNS = {
    "human": [
        "EnsemblGeneID_mapping_version_20200130",
        "Description",
        "MitoCarta3.0_Evidence",
        "MitoCarta3.0_SubMitoLocalization",
        "MitoCarta3.0_MitoPathways",
        "HPA_Main_Location_2020 (Reliability)",
        "Tissues",
    ],
    "mouse": [
        "EnsemblGeneID",
        "Description",
        "MitoCarta3.0_Evidence",
        "MitoCarta3.0_SubMitoLocalization",
        "MitoCarta3.0_MitoPathways",
        "HPA_Main_Location_2020 (Reliability)",
        "Tissues",
    ],
}

MITO_PATHWAYS = "mito_pathways"
MITO_ENSEMBL_ID = "ensembl_id"

MITOCART_COL_MAPPER = {
    "EnsemblGeneID_mapping_version_20200130": TARGET_COL,
    "EnsemblGeneID": TARGET_COL,
    "Description": "gene_description",
    "MitoCarta3.0_Evidence": "evidence",
    "MitoCarta3.0_SubMitoLocalization": "sub_mito_localization",
    "MitoCarta3.0_MitoPathways": MITO_PATHWAYS,
    "HPA_Main_Location_2020 (Reliability)": "hpa_location",
    "Tissues": "tissue_expression",
}
MITOCART_OUTPUT = [
    "gene_description",
    "evidence",
    "sub_mito_localization",
    "hpa_location",
    "tissue_expression",
    "mito_pathways",
]


"""
MolMeDB variables
"""
MOLMEDB_COMPOUND_NAME = f"{MOLMEDB}_compound_name"
MOLMEDB_INCHIKEY = f"{MOLMEDB}_inchikey"
MOLMEDB_SMILES = f"{MOLMEDB}_smiles"
MOLMEDB_ID = f"{MOLMEDB}_id"

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
MOLMEDB_UNIPROT_ID = f"{MOLMEDB}_uniprot_trembl_id"
MOLMEDB_HGNC_ID = "hgcn_id"
MOLMEDB_HGNC_SYMBOL = f"{MOLMEDB}_hgnc_symbol"

MAP_GENE_COL_NAMES = {
    MOLMEDB_INHIBITOR_INCHIKEY: TARGET_COL,
    MOLMEDB_HGNC_ID: MOLMEDB_HGNC_SYMBOL,
    MOLMEDB_UNIPROT_TREMBL_ID: MOLMEDB_UNIPROT_ID,
}

MOLMEDB_COMPOUND_PROTEIN_OUTPUT_DICT = {
    MOLMEDB_UNIPROT_ID: str,
    MOLMEDB_HGNC_SYMBOL: str,
    SOURCE_PMID: str,
}
MOLMEDB_UNIPROT_TREMBL_DEFAULT_ID = UNIPROT_TREMBL

"""
Open Targets variables
"""
OPENTARGETS_GO_ID = GO_ID
OPENTARGETS_GO_NAME = GO_NAME
OPENTARGETS_GO_TYPE = GO_TYPE

OPENTARGETS_GO_OUTPUT_DICT = {
    OPENTARGETS_GO_ID: str,
    OPENTARGETS_GO_NAME: str,
    OPENTARGETS_GO_TYPE: str,
}

OPENTARGETS_POSSIBLE_PATHWAY_IDS = f"{REACTOME}|{WP}"
OPENTARGETS_REACTOME_OUTPUT_DICT = {
    PATHWAY_ID: str,
    PATHWAY_LABEL: str,
}

OPENTARGETS_COMPOUND_RELATION = "relation"
OPENTARGETS_ADVERSE_EFFECT_COUNT = "adverse_effect_count"
OPENTARGETS_ADVERSE_EFFECT = "adverse_effect"
OPENTARET_COMPOUND_COLS = [
    "chembl_id",
    "compound_name",
    "is_approved",
    "clincal_trial_phase",
    "cross_references",
    "adverse_events",
]
OPENTARGETS_COMPOUND_OUTPUT_DICT = {
    CHEMBL_ID: str,
    DRUGBANK_ID: str,
    "compound_cid": str,
    "compound_name": str,
    "clincal_trial_phase": int,
    "is_approved": bool,
    "relation": str,
    OPENTARGETS_ADVERSE_EFFECT_COUNT: int,
    OPENTARGETS_ADVERSE_EFFECT: list,
}
OPENTARGETS_COMPOUND_VALUE_CHECK_LIST = [CHEMBL, DRUGBANK, PUBCHEM_COMPOUND_CID]
OPENTARGETS_COMPOUND_DISEASE_RELATION = "treats"

"""
PubChem variables
"""
ASSAY_ENDPOINT_TYPES = {
    "http://www.bioassayontology.org/bao#BAO_0000188": "EC50",
    "http://www.bioassayontology.org/bao#BAO_0000190": "IC50",
    "http://www.bioassayontology.org/bao#BAO_0002146": "MIC",
}
PUBCHEM_ASSAY_ID = "pubchem_assay_id"
PUBCHEM_COMPOUND_NAME = "compound_name"
PUBCHEM_SMILES = "smiles"
PUBCHEM_INCHI = "inchi"
PUBCHEM_UNIPROT_IRI = "http://purl.uniprot.org/uniprot/"
PUBCHEM_BIOASSAY_IRI = "http://rdf.ncbi.nlm.nih.gov/pubchem/bioassay/"
PUBCHEM_OUTCOME_IRI = "http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#"
PUBCHEM_COMPOUND_IRI = "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/"

PUBCHEM_COMPOUND_OUTPUT_DICT = {
    PUBCHEM_ASSAY_ID: str,
    "assay_type": str,
    "outcome": str,
    "compound_cid": str,
    PUBCHEM_COMPOUND_NAME: str,
    PUBCHEM_SMILES: str,
    PUBCHEM_INCHI: str,
}
PUBCHEM_POSSIBLE_OUTCOMES = "active|inactive"

"""
StringDB variables
"""
STRING_PREFERRED_NAME_A = "preferredName_A"
STRING_PREFERRED_NAME_B = "preferredName_B"
STRING_PPI_LINK_TO = "stringdb_link_to"
STRING_PPI_SCORE = "score"
UNIPROT_TREMBL_A = "Uniprot-TrEMBL_A"
UNIPROT_TREMBL_B = "Uniprot-TrEMBL_B"
STRING_PPI_INTERACTS_WITH = "interacts_with"
STRING_OUTPUT_DICT = {
    STRING_PPI_LINK_TO: str,
    STRING_GENE_INPUT_ID: str,
    STRING_PPI_SCORE: float,
    UNIPROT_TREMBL_A: str,
    UNIPROT_TREMBL_B: str,
}

"""
TFLink variables
"""
ITS_TARGET_COL = "its_target"
ITS_TF_COL = "its_tf"
IS_TF_COL = "is_tf"
IS_TARGET_COL = "is_target"

TFLINK_GENE_ID_TF = "NCBI.GeneID.TF"
TFLINK_GENE_ID_TARGET = "NCBI.GeneID.Target"
ENSEMBL_GENE_ID_TF = "Ensembl.GeneID.TF"
ENSEMBL_GENE_ID_TARGET = "Ensembl.GeneID.Target"

TF_COLS_TO_KEEP = [
    TFLINK_GENE_ID_TARGET,
    "Ensembl.GeneID.Target",
    "Name.Target",
    "UniprotID.Target",
    "Target.TFLink.ortho",
    "Target.nonTFLink.ortho",
    "Detection.method",
    "PubmedID",
    "Source.database",
    "Small-scale.evidence",
]

TARGET_COLS_TO_KEEP = [
    "NCBI.GeneID.TF",
    "Ensembl.GeneID.TF",
    "Name.TF",
    "UniprotID.TF",
    "TF.TFLink.ortho",
    "TF.nonTFLink.ortho",
    "Detection.method",
    "PubmedID",
    "Source.database",
    "Small-scale.evidence",
]

"""
Wikidata variables
"""
WIKIDATA_ID_COL = "wikidata_id"
WIKIDATA_LABEL_COL = "wikidata_label"
WIKIDATA_GO_COL = "go_id"
WIKIDATA_OUTPUT_DICT = {
    WIKIDATA_ID_COL: str,
    WIKIDATA_LABEL_COL: str,
    WIKIDATA_GO_COL: str,
}

"""
WikiPathways variables
"""
WIKIPATHWAYS_GENE_ID = "gene_id"

WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT = {
    PATHWAY_ID: str,
    PATHWAY_LABEL: str,
    PATHWAY_GENE_COUNTS: int,
}

WIKIPATHWAYS_TARGET_GENE = "targetGene"
WIKIPATHWAYS_TARGET_PROTEIN = "targetProtein"
WIKIPATHWAYS_TARGET_METABOLITE = "targetMetabolite"
WIKIPATHWAYS_MIM_TYPE = "mimtype"
WIKIPATHWAYS_RHEA_ID = "rhea_id"

WIKIPATHWAYS_MOLECULAR_GENE_OUTPUT_DICT = {
    PATHWAY_ID: str,
    PATHWAY_LABEL: str,
    WIKIPATHWAYS_TARGET_GENE: str,
    WIKIPATHWAYS_TARGET_PROTEIN: str,
    WIKIPATHWAYS_TARGET_METABOLITE: str,
    WIKIPATHWAYS_MIM_TYPE: str,
    WIKIPATHWAYS_RHEA_ID: str,
}

WIKIPATHWAY_ID_CLEANER_DICT = {
    WIKIPATHWAYS_GENE_ID: "https://identifiers.org/ncbigene/",
    WIKIPATHWAYS_TARGET_GENE: "https://identifiers.org/ncbigene/",
    WIKIPATHWAYS_TARGET_PROTEIN: "https://identifiers.org/uniprot/",
    WIKIPATHWAYS_TARGET_METABOLITE: "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/",
    WIKIPATHWAYS_MIM_TYPE: "http://vocabularies.wikipathways.org/wp#",
}

WIKIPATHWAY_NAMESPACE_DICT = {
    WIKIPATHWAYS_GENE_ID: NCBI_GENE_ID,
    WIKIPATHWAYS_TARGET_GENE: NCBI_GENE_ID,
    WIKIPATHWAYS_TARGET_PROTEIN: UNIPROT_TREMBL,
    WIKIPATHWAYS_TARGET_METABOLITE: PUBCHEM_COMPOUND_CID,
}

"""
AOPWIKI variables
"""
AOP = "aop"
AOP_TITLE = "aop_title"
MIE_TITLE = "MIE_title"
MIE = "MIE"
KE_DOWNSTREAM = "KE_downstream"
KE_DOWNSTREAM_TITLE = "KE_downstream_title"
KER = "KER"
AO = "ao"
AO_TITLE = "ao_title"
KE_UPSTREAM = "KE_upstream"
KE_UPSTREAM_TITLE = "KE_upstream_title"
KE_UPSTREAM_ORGAN = "KE_upstream_organ"
KE_DOWNSTREAM_ORGAN = "KE_downstream_organ"
PUBCHEM_COMPOUND_AOPWIKI = "pubchem_compound"

AOPWIKI_GENE_OUTPUT_DICT = {
    AOP: str,
    AOP_TITLE: str,
    MIE_TITLE: str,
    MIE: str,
    KE_DOWNSTREAM: str,
    KE_DOWNSTREAM_TITLE: str,
    KER: str,
    AO: str,
    AO_TITLE: str,
    KE_UPSTREAM: str,
    KE_UPSTREAM_TITLE: str,
    KE_UPSTREAM_ORGAN: str,
    KE_DOWNSTREAM_ORGAN: str,
    PUBCHEM_COMPOUND_AOPWIKI: str,
}

AOPWIKI_COMPOUND_OUTPUT_DICT = {
    AOP: str,
    AOP_TITLE: str,
    MIE_TITLE: str,
    MIE: str,
    KE_DOWNSTREAM: str,
    KE_DOWNSTREAM_TITLE: str,
    KER: str,
    AO: str,
    AO_TITLE: str,
    KE_UPSTREAM: str,
    KE_UPSTREAM_TITLE: str,
    KE_UPSTREAM_ORGAN: str,
    KE_DOWNSTREAM_ORGAN: str,
}


# QUERY_PROCESS = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-get-by-biological-process.rq.rq")

AOPWIKI_GENE_PARAMS = {
    "type_of_input": "genes",
    "uri": "<https://identifiers.org/ensembl/",
    "column": AOPWIKI_GENE_INPUT_ID,
    "input_col": AOPWIKI_GENE_INPUT_ID,
    "output_dict": AOPWIKI_GENE_OUTPUT_DICT,
    "query_file_name": "aopwiki-gene.rq",
    "output_col": AOPWIKI_GENE_COL,
}
AOPWIKI_COMPOUND_PARAMS = {
    "type_of_input": "compounds",
    "uri": "<https://identifiers.org/pubchem.compound/",
    "column": AOPWIKI_COMPOUND_INPUT_ID,
    "input_col": "pubchem_compound",
    "output_dict": AOPWIKI_COMPOUND_OUTPUT_DICT,
    "query_file_name": "aopwiki-compound.rq",
    "output_col": AOPWIKI_COMPOUND_COL,
}

# TODO: Look into this
ENSEMBL_HOMOLOGS = "Ensembl_homologs"

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
PHENOTYPE_NODE_LABEL = "Phenotype"
MIRNA_NODE_LABEL = "miRNA"
TRANS_FACTOR_NODE_LABEL = "Transcription Factor"
MITOCHONDRIAL_PATHWAY_NODE_LABEL = "Mitochondrial Pathway"
KEY_EVENT_NODE_LABEL = "Key Event"
MIE_NODE_LABEL = "Molecular Initiating Event"
AOP_NODE_LABEL = "Adverse Outcome Pathway"
AO_NODE_LABEL = "Adverse Outcome"

# Edge types (Neo4j)
INTERACTS_WITH = "INTERACTS_WITH"
PART_OF = "PART_OF"
ASSOCIATED_WITH = "ASSOCIATED_WITH"
ACTIVATES = "ACTIVATES"
HAS_SIDE_EFFECT = "HAS_SIDE_EFFECT"
TREATS = "TREATS"
INHIBITS = "INHIBITS"
EXPRESSED_BY = "EXPRESSED_BY"
UPSTREAM_OF = "UPSTREAM_OF"
DOWNSTREAM_OF = "DOWNSTREAM_OF"

"""
Anatomical entity nodes and edges
"""
BGEE_ANATOMICAL_NODE_MAIN_LABEL = ANATOMICAL_ID
BGEE_ANATOMICAL_NODE_ATTRS = {
    DATASOURCE: BGEE,
    NAME: None,
    ID: None,
    UBERON: None,
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
Disease nodes and edges
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

OPENTARGET_DISEASE_NODE_ATTRS = {
    NAME: None,
    ID: None,
    DATASOURCE: OPENTARGETS,
    LABEL: DISEASE_NODE_LABEL,
}

OPENTARGETS_DISEASE_COMPOUND_EDGE_ATTRS = {
    DATASOURCE: OPENTARGETS,
    LABEL: None,
}

"""
Pathway nodes and edges
"""
GENE_COUNTS = "gene_counts"
PATHWAY_NODE_ATTRS = {
    DATASOURCE: None,
    NAME: None,
    ID: None,
    LABEL: PATHWAY_NODE_LABEL,
    GENE_COUNTS: None,
}

GENE_PATHWAY_EDGE_LABEL = "part_of"
GENE_PATHWAY_EDGE_ATTRS = {
    DATASOURCE: None,
    LABEL: GENE_PATHWAY_EDGE_LABEL,
}

# GO nodes
# Open Targets - GO processes
GO_BP_NODE_LABEL = "Biological Process"
GO_MF_NODE_LABEL = "Molecular Function"
GO_CC_NODE_LABEL = "Cellular Component"
GO_NODE_MAIN_LABEL = "go_id"
GO_NODE_ATTRS = {
    "datasource": OPENTARGETS,
    "name": None,
    "id": None,
    "labels": None,
}
GENE_GO_EDGE_LABEL = "part_of"
GENE_GO_EDGE_ATTRS = {"datasource": OPENTARGETS, "label": GENE_GO_EDGE_LABEL}

# IntAct interactions
SPECIES = "species"
# INTACT_INTERACTION_TYPE = "interaction_type"
# INTACT_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
# INTACT_NODE_ATTRS.update(
#     {
#         DATASOURCE: INTACT,
#         ID: None,
#         LABEL: None,
#     }
# )

KEGG_PATHWAY_NODE_MAIN_LABEL = PATHWAY_ID
KEGG_PATHWAY_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
KEGG_PATHWAY_NODE_ATTRS.update(
    {
        DATASOURCE: KEGG,
        GENE_COUNTS: None,
    }
)

MINERVA_PATHWAY_NODE_MAIN_LABEL = PATHWAY_ID
MINERVA_PATHWAY_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
MINERVA_PATHWAY_NODE_ATTRS.update(
    {
        DATASOURCE: MINERVA,
        GENE_COUNTS: None,
    }
)

OPENTARGETS_GO_NODE_MAIN_LABEL = GO_ID
OPENTARGETS_GO_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
OPENTARGETS_GO_NODE_ATTRS.update(
    {
        DATASOURCE: OPENTARGETS,
    }
)

OPENTARGETS_GENE_GO_EDGE_ATTRS = GENE_PATHWAY_EDGE_ATTRS.copy()
OPENTARGETS_GENE_GO_EDGE_ATTRS.update({DATASOURCE: OPENTARGETS})

OPENTARGETS_REACTOME_NODE_MAIN_LABEL = PATHWAY_ID
OPENTARGETS_REACTOME_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
OPENTARGETS_REACTOME_NODE_ATTRS.update({DATASOURCE: OPENTARGETS})
OPENTARGETS_GENE_REACTOME_EDGE_ATTRS = GENE_PATHWAY_EDGE_ATTRS.copy()
OPENTARGETS_GENE_REACTOME_EDGE_ATTRS.update({DATASOURCE: OPENTARGETS})

WIKIPATHWAYS_NODE_MAIN_LABEL = PATHWAY_ID
WIKIPATHWAYS_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
WIKIPATHWAYS_NODE_ATTRS.update({DATASOURCE: WIKIPATHWAYS, GENE_COUNTS: None})

MOLECULAR_PATHWAY_NODE_MAIN_LABEL = PATHWAY_ID
MOLECULAR_PATHWAY_NODE_ATTRS = PATHWAY_NODE_ATTRS.copy()
MOLECULAR_PATHWAY_NODE_ATTRS.update({DATASOURCE: WIKIPATHWAYS})

WIKIPATHWAYS_INTERACTION_TYPE = "interaction_type"
MOLECULAR_INTERACTION_EDGE_ATTRS = GENE_PATHWAY_EDGE_ATTRS.copy()  # type: ignore
MOLECULAR_INTERACTION_EDGE_ATTRS.update(
    {
        WIKIPATHWAYS_INTERACTION_TYPE: None,
        WIKIPATHWAYS_RHEA_ID: None,
        DATASOURCE: WIKIPATHWAYS,
    }
)

# GProfiler
GPROFILER_GENE_NODE_MAIN_LABEL = "name"
GPROFILER_TERM_SIZE = "term_size"
P_VALUE = "p_value"
SIGNIFICANT = "significant"
GPROFILER_GENE_NODE_ATTRS = {
    NAME: None,
    ID: None,
    GENE_COUNTS: None,
    P_VALUE: None,
    SIGNIFICANT: None,
    LABEL: None,
    DATASOURCE: GPROFILER,
}

GPROFILER_EDGE_ATTRS = {
    DATASOURCE: GPROFILER,
    LABEL: None,
}

GPROFILER_PATHWAY_TYPE = "pathway_type"
GPROFILE_GENE_HP_EDGE_LABEL = "linked_to"
GPROFILE_GENE_HPA_EDGE_LABEL = "expressed_in"
GPROFILER_GENE_MIRNA_EDGE_LABEL = "regulates"
GPROFILER_GENE_TF_EDGE_LABEL = "regulated_by"

# Mitocart
MITOCART_NODE_MAIN_LABEL = "mito_pathways"
MITOCART_HPA_LOCATION = "hpa_location"
MITOCART_SUB_MITO_LOCALIZATION = "sub_mito_localization"
EVIDENCE = "evidence"

MITOCART_PATHWAY_LABEL = "Mitochondrial_Pathway"
MITOCART_EDGE_LABEL = "encodes_mitochondrial_protein"

MITOCART_NODE_ATTRS = {
    DATASOURCE: MITOCARTA,
    NAME: None,
    ID: None,
    LABEL: MITOCHONDRIAL_PATHWAY_NODE_LABEL,
    MITOCART_HPA_LOCATION: None,
    MITOCART_SUB_MITO_LOCALIZATION: None,
    EVIDENCE: None,
}

MITOCART_GENE_PATHWAY_EDGE_LABEL = "encodes_mitochondrial_protein"

"""
Compound nodes and edges
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
    DATASOURCE: MOLMEDB,
    LABEL: MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
}

# Opentargets
OPENTARGETS_COMPOUND_NODE_MAIN_LABEL = "compound_cid"
OPENTARGETS_COMPOUND_CID = "compound_cid"
OPENTARGETS_COMPOUND_CLINICAL_TRIAL_PHASE = "clincal_trial_phase"
OPENTARGETS_COMPOUND_IS_APPROVED = "is_approved"

OPENTARGETS_COMPOUND_NODE_ATTRS = {
    DATASOURCE: OPENTARGETS,
    NAME: None,
    ID: None,
    DRUGBANK_ID: None,
    OPENTARGETS_COMPOUND_CID: None,
    OPENTARGETS_COMPOUND_CLINICAL_TRIAL_PHASE: None,
    OPENTARGETS_COMPOUND_IS_APPROVED: None,
    OPENTARGETS_ADVERSE_EFFECT_COUNT: None,
    LABEL: COMPOUND_NODE_LABEL,
}
OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS = {DATASOURCE: OPENTARGETS, LABEL: None, ID: None}

"""
Side effect nodes and edges
"""
# Opentarget
COMPOUND_SIDE_EFFECT_NODE_MAIN_LABEL = "adverse_effect"
COMPOUND_SIDE_EFFECT_NODE_LABEL = "name"
COMPOUND_SIDE_EFFECT_NODE_ATTRS = {
    DATASOURCE: OPENTARGETS,
    NAME: None,
    ID: None,
    LABEL: SIDE_EFFECT_NODE_LABEL,
}

# Opentarget
COMPOUND_SIDE_EFFECT_EDGE_LABEL = "has_side_effect"
COMPOUND_SIDE_EFFECT_EDGE_ATTRS = {
    DATASOURCE: OPENTARGETS,
    LABEL: COMPOUND_SIDE_EFFECT_EDGE_LABEL,
}

# PubChem - Assays (TODO: check this)
PUBCHEM_COMPOUND_NODE_ATTRS = {
    DATASOURCE: PUBCHEM,
    NAME: None,
    ID: None,
    INCHI: None,
    SMILES: None,
    LABEL: COMPOUND_NODE_LABEL,
}

PUBCHEM_ASSAY_TYPE = "assay_type"
PUBCHEM_ASSAY_ID = "pubchem_assay_id"
PUBCHEM_EDGE_LABEL_MAPPER = {
    "active": "activates",
    "inactive": "inhibits",
}

PUBCHEM_GENE_COMPOUND_EDGE_ATTRS = {
    DATASOURCE: PUBCHEM,
    PUBCHEM_ASSAY_TYPE: None,
    PUBCHEM_ASSAY_ID: None,
    LABEL: None,
}

"""
PPI interactions
"""
INTACT_PPI_EDGE_MAIN_LABEL = "intact_link_to"
INTACT_PPI_EDGE_ATTRS = {
    DATASOURCE: INTACT,
    INTACT_DETECTION_METHOD: None,
    INTACT_TYPE: None,
    LABEL: INTACT_PPI_EDGE_MAIN_LABEL,
}

# STRING
STRING_PPI_EDGE_MAIN_LABEL = "interacts_with"
STRING_PPI_EDGE_ATTRS = {
    DATASOURCE: STRING,
    STRING_PPI_SCORE: None,
    LABEL: STRING_PPI_EDGE_MAIN_LABEL,
}

"""
Other nodes and edges
"""

# Ensembl Homologs TODO: check this
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

# TFLink
NAME_TARGET = "name_target"
UNIPROTID_TARGET = "uniprotid_target"
DETECTION_METHOD = "detection_method"
PUBMEDID = "pubmedid"
SOURCE_DATABASE = "source_database"
SMALL_SCALE_EVIDENCE = "small_scale_evidence"

TFLINK_EDGE_LABEL = "tf_regulates"
TFLINK_GENE_TF_EDGE_ATTRS = {
    DATASOURCE: TFLINK,
    LABEL: TFLINK_EDGE_LABEL,
    NAME_TARGET: None,
    UNIPROTID_TARGET: None,
    DETECTION_METHOD: None,
    PUBMEDID: None,
    SOURCE_DATABASE: None,
    SMALL_SCALE_EVIDENCE: None,
}


# Literature
LITERATURE_NODE_MAIN_LABEL = UMLS
LITERATURE_DISEASE_NODE_ATTRS = {
    DATASOURCE: None,
    NAME: None,
    ID: None,
    MONDO: None,
    UMLS: None,
    LABEL: DISEASE_NODE_LABEL,
}
LITERATURE_DISEASE_EDGE_ATTRS = {
    DATASOURCE: None,
    LABEL: GENE_DISEASE_EDGE_LABEL,
}

# Wikidata


# AOPWIKI
AOP_NODE_MAIN_LABEL = "aop"
AOP_EDGE_LABEL = "associated_with"  # Gene to AOP edge
MIE_NODE_MAIN_LABEL = "MIE"
MIE_AOP_EDGE_LABEL = "upstream_of"  # MIE to AOP edge
KEY_EVENT_UPSTREAM_NODE_MAIN_LABEL = "KE_upstream"
KE_UPSTREAM_MIE_EDGE_LABEL = "upstream_of"  # KE upstream to MIE edge
KEY_EVENT_DOWNSTREAM_NODE_MAIN_LABEL = "KE_downstream"
KE_DOWNSTREAM_KE_EDGE_LABEL = "downstream_of"  # KE downstream to KE upstream edge
AO_NODE_MAIN_LABEL = "ao"
AO_KE_EDGE_LABEL = "associated_with"  # KE downstream to AO edge

AOPWIKI_EDGE_ATTRS = {
    DATASOURCE: AOPWIKIRDF,
    LABEL: None,
}

AOPWIKI_NODE_ATTRS = {
    DATASOURCE: AOPWIKIRDF,
    NAME: None,
    ID: None,
    LABEL: None,
    "organ": None,
}

"""
Cytoscape variables
"""
NODE_TYPE = "node_type"
SOURCE = "source"
TARGET = "target"
INTERACTION = "interaction"

ALL_NODE_LABELS = {
    GENE_NODE_LABEL: "ELLIPSE",
    TRANS_FACTOR_NODE_LABEL: "ELLIPSE",
    HOMOLOG_NODE_LABEL: "ELLIPSE",
    DISEASE_NODE_LABEL: "HEXAGON",
    COMPOUND_NODE_LABEL: "DIAMOND",
    ANATOMICAL_NODE_LABEL: "RECTANGLE",
    PATHWAY_NODE_LABEL: "ROUND_RECTANGLE",
    GO_BP_NODE_LABEL: "ROUND_RECTANGLE",
    GO_MF_NODE_LABEL: "ROUND_RECTANGLE",
    GO_CC_NODE_LABEL: "ROUND_RECTANGLE",
    SIDE_EFFECT_NODE_LABEL: "OCTAGON",
    PHENOTYPE_NODE_LABEL: "HEXAGON",
    MIRNA_NODE_LABEL: "HEXAGON",
    MITOCHONDRIAL_PATHWAY_NODE_LABEL: "ROUND_RECTANGLE",
    KEY_EVENT_NODE_LABEL: "OCTAGON",
    MIE_NODE_LABEL: "ROUND_RECTANGLE",
    AOP_NODE_LABEL: "ROUND_RECTANGLE",
    AO_NODE_LABEL: "ROUND_RECTANGLE",
}

COLOR_MAPPER = {
    GENE_NODE_LABEL: "#42d4f4",
    DISEASE_NODE_LABEL: "#4363d8",
    COMPOUND_NODE_LABEL: "#e6194B",
    ANATOMICAL_NODE_LABEL: "#ff7b00",
    PATHWAY_NODE_LABEL: "#ffa652",
    GO_BP_NODE_LABEL: "#ffcd90",
    GO_MF_NODE_LABEL: "#3cb44b",
    GO_CC_NODE_LABEL: "#ffd700",
    SIDE_EFFECT_NODE_LABEL: "#aaffc3",
    HOMOLOG_NODE_LABEL: "#9b59b6",
    PHENOTYPE_NODE_LABEL: "#000075",
    MIRNA_NODE_LABEL: "#e91e63",
    TRANS_FACTOR_NODE_LABEL: "#8e24aa",
    MITOCHONDRIAL_PATHWAY_NODE_LABEL: "#00bcd4",
    KEY_EVENT_NODE_LABEL: "#ff5722",
    MIE_NODE_LABEL: "#795548",
    AOP_NODE_LABEL: "#607d8b",
    AO_NODE_LABEL: "#f44336",
}

"""
RDF variables
"""
# Mapper from namespace to BridgeDB datasource
COMPOUND_NAMESPACE_MAPPER = {"pubchem.compound": "PubChem Compound", "CHEMBL": "ChEMBL compound"}


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
    "Opentargets_reactome": OPENTARGETS_REACTOME_COL,
}

NODE_URI_PREFIXES = {
    ENSEMBL: "https://identifiers.org/ensembl#",
    "medgen": "https://www.ncbi.nlm.nih.gov/medgen/",
    "pubchem_assay": "https://pubchem.ncbi.nlm.nih.gov/bioassay/",
    PUBCHEM: "https://pubchem.ncbi.nlm.nih.gov/compound/",
    MOLMEDB: "https://molmedb.upol.cz/mol/",
    "uniprot": "https://www.uniprot.org/uniprotkb/",
    "pubmed": "https://pubmed.ncbi.nlm.nih.gov/",
    WIKIPATHWAYS: "https://www.wikipathways.org/pathways/",
    REACTOME: "https://reactome.org/content/detail/",
    MINERVA: "https://minerva-net.lcsb.uni.lu/api/",
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
    "interaction": "interaction",
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
    "ao": "http://aopkb.org/aop_ontology#AdverseOutcome",
    "ke": "http://aopkb.org/aop_ontology#KeyEvent",
    "mie": "http://aopkb.org/aop_ontology#MolecularInitiatingEvent",
    "ker": "http://aopkb.org/aop_ontology#KeyEventRelationship",
    "interaction": "https://vocabularies.wikipathways.org/wp#Interaction",
}

# PREDICATES
PREDICATES = {
    "sameAs": f"{NAMESPACE_BINDINGS['owl']}sameAs",
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
    "has_upstreamkey_event": "http://aopkb.org/aop_ontology#has_upstream_key_event",
    "has_downstreamkey_event": "http://aopkb.org/aop_ontology#has_downstream_key_event",
    "occurs_in": f"{NAMESPACE_BINDINGS['obo']}BFO_0000066",
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
    OPENTARGETS: "https://www.opentargets.org/",
    AOPWIKIRDF: "https://aopwiki.rdf.bigcat-bioinformatics.org",
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
    "aopwiki": "https://identifiers.org/",
    "WikiPathways": "https://www.wikipathways.org/",
    "Reactome": "https://reactome.org/",
    "Minerva": "https://minerva-net.lcsb.uni.lu/",
    # TODO ADD ALL
}
