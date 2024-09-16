# coding: utf-8

"""Python module to produce an RDF graph from the property table."""

import logging
from importlib import resources
import pandas as pd
from bioregistry import get_iri, normalize_curie
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, SKOS, XSD, DC, DCTERMS
import numpy as np
from datetime import datetime
from pyBiodatafuse.constants import (
    BGEE_GENE_EXPRESSION_LEVELS_COL,
    DISGENET_DISEASE_COL,
    NAMESPACE_BINDINGS,
    NODE_TYPES,
    PREDICATES,
    URIS,
    CLINICAL_PHASES,
    GO_TYPES,
    MOAS,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def replace_na_none(item):
    """
    Recursively replaces 'na', 'nan', 'none' strings, Python None, and NaN values with None.

    Handles:
    - Strings: Converts 'na', 'nan', 'none' to None.
    - None and NaN: Converts to None.
    - Lists, dictionaries, and NumPy arrays: Recursively processes elements.

    Args:
        item: The value to process (string, None, NaN, list, dict, or np.ndarray).

    Returns:
        The processed item with 'na', 'nan', 'none', None, and NaN replaced by None.
    """
    # Replace 'na' strings
    if isinstance(item, str) and item.lower() in ["na", "nan", "none"]:
        return None
    # Replace Python None (and handle NaNs)
    elif item is None or (isinstance(item, float) and pd.isna(item)):
        return None
    # Recursively handle lists
    elif isinstance(item, list):
        return [replace_na_none(sub_item) for sub_item in item]
    # Recursively handle dictionaries
    elif isinstance(item, dict):
        return {key: replace_na_none(value) for key, value in item.items()}
    # Recursively handle NumPy arrays
    elif isinstance(item, np.ndarray):
        return np.array([replace_na_none(sub_item) for sub_item in item], dtype=object)
    # Return the item as-is if no replacement needed
    else:
        return item


def add_gene_node(g: Graph, row) -> URIRef:
    """Create and add a gene node to the RDF graph.

    :param g: RDF graph to which the gene node will be added.
    :param row: DataFrame row containing gene information.
    :return: URIRef for the created gene node.
    """
    target = row["target"]
    if target:
        gene_node = URIRef(f"http://identifiers.org/ensembl/{target}")
        g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((gene_node, RDFS.label, Literal(row["identifier"], datatype=XSD.string)))
        return gene_node
    else:
        return None


def add_disease_node(g: Graph, disease_data: dict) -> URIRef:
    """Create and add a disease node to the RDF graph.

    :param g: RDF graph to which the disease node will be added.
    :param disease_data: Dictionary containing disease information.
    :return: URIRef for the created disease node.
    """
    # UMLS IRIs not in Bioregistry
    disease_curie = disease_data.get("disease_umlscui")
    if disease_curie is None:
        disease_curie = disease_data.get("UMLS")
    disease_iri = f"https://www.ncbi.nlm.nih.gov/medgen/{disease_curie}"
    disease_node = URIRef(disease_iri)
    g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
    g.add(
        (disease_node, RDFS.label, Literal(disease_data.get("disease_name"), datatype=XSD.string))
    )
    for identifier_type in [
        "HPO",
        "NCI",
        "OMIM",
        "MONDO",
        "ORDO",
        "EFO",
        "DO",
        "MESH",
        "UMLS",
    ]:
        if disease_data.get(identifier_type):
            curie_field = disease_data.get(identifier_type)
            if "," in curie_field:
                curies = [i for i in curie_field.split(", ")]
            else:
                curies = [curie_field]
            for curie in curies:
                disease_source_iri = get_iri(curie)
                if disease_source_iri is None:
                    if ":" in curie:
                        curie = curie.split(":")[1]
                    disease_source_iri = get_iri("obo:" + curie)
                    g.add(
                        (disease_node, SKOS.closeMatch, URIRef(disease_source_iri))
                    )  # Some of the data does not look like a skos:exactMatch
                else:
                    g.add((disease_node, SKOS.closeMatch, URIRef(disease_source_iri)))
        else:
            pass
    return disease_node


def add_gene_disease_association_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, new_uris: dict
) -> URIRef:
    """Create and add a gene-disease association node to the RDF graph.

    :param g: RDF graph to which the association node will be added.
    :param id_number: Unique identifier for the association node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the gene.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :return: URIRef for the created association node.
    """
    gene_disease_assoc_node = URIRef(
        f"{new_uris['gene_disease_association']}/{id_number}/{source_idx}_assoc_{disease_id}"
    )
    g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene_disease_association"])))
    return gene_disease_assoc_node


def add_score_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, score: float, new_uris: dict, i
) -> URIRef:
    """Create and add a score node to the RDF graph.

    :param g: RDF graph to which the score node will be added.
    :param id_number: Unique identifier for the score node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the score.
    :param score: Score value for the gene-disease association.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: Index or iteration number used for node URI construction.
    :return: URIRef for the created score node.
    """
    score_node = URIRef(f"{new_uris['score_base_node']}/{id_number}/{i}/{source_idx}_{disease_id}")
    g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
    g.add(
        (
            score_node,
            URIRef(NAMESPACE_BINDINGS["sio"] + "has_value"),
            Literal(score, datatype=XSD.double),
        )
    )
    return score_node


def add_evidence_idx_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, ei: float, new_uris: dict, i
) -> URIRef:
    """Create and add an evidence index (EI) node to the RDF graph.

    :param g: RDF graph to which the EI node will be added.
    :param id_number: Unique identifier for the EI node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the EI.
    :param ei: Evidence index value for the gene-disease association.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: Index or iteration number used for node URI construction.
    :return: URIRef for the created EI node.
    """
    evidence_idx_node = URIRef(
        f"{new_uris['score_base_node']}/{id_number}/{i}/{source_idx}_{disease_id}"
    )
    g.add((evidence_idx_node, RDF.type, URIRef(NODE_TYPES["evidence_idx_node"])))
    g.add(
        (
            evidence_idx_node,
            URIRef(NAMESPACE_BINDINGS["sio"] + "has_value"),
            Literal(ei, datatype=XSD.double),
        )
    )
    return evidence_idx_node


def add_data_source_node(g: Graph) -> URIRef:
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :return: URIRef for the created data source node.
    """
    # for source in sources (eg disgenet)
    data_source_name = Literal("DisGeNET", datatype=XSD.string)
    data_source_url = URIRef("https://disgenet.com/")
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url


def add_gene_disease_associations(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    disease_data: list,
    new_uris: dict,
    i: int,
) -> None:
    """Process and add gene-disease association data to the RDF graph.

    :param g: RDF graph to which the associations will be added.
    :param id_number: Unique identifier for the associations.
    :param source_idx: Source index for the associations.
    :param gene_node: URIRef of the gene node associated with the disease data.
    :param disease_data: List of dictionaries containing disease association information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    """

    data = disease_data
    disease_node = add_disease_node(g, data)
    if disease_node:
        gene_disease_assoc_node = URIRef(
            f"{new_uris['gene_disease_association']}/{id_number}/{i}/{source_idx}"
        )
        g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene_disease_association"])))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), gene_node))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), disease_node))
        if data.get("score"):
            score_node = add_score_node(
                g=g,
                id_number=id_number,
                source_idx=source_idx,
                score=data["score"],
                new_uris=new_uris,
                i=i,
                disease_id=disease_data.get("disease_umlscui"),
            )  # SEE
            g.add(
                (
                    gene_disease_assoc_node,
                    URIRef(PREDICATES["sio_has_measurement_value"]),
                    score_node,
                )
            )
        data_source_node = add_data_source_node(g)
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_has_source"]), data_source_node))


def add_gene_expression_data(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    expression_data: list,
    experimental_process_data: list,
    new_uris: dict,
) -> None:
    """Process and add gene expression data to the RDF graph.

    :param g: RDF graph to which the expression data will be added.
    :param id_number: Unique identifier for the expression data.
    :param source_idx: Source index for the expression data.
    :param gene_node: URIRef of the gene node associated with the expression data.
    :param expression_data: List of dictionaries containing gene expression information.
    :param experimental_process_data: List of experimental process data.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    """
    for data in expression_data:
        if pd.isna(data["anatomical_entity_id"]):
            continue

        gene_expression_value_node = URIRef(
            f"{new_uris['gene_expression_value_base_node']}/{id_number}/{source_idx}"
        )
        developmental_stage_node = URIRef(get_iri(data["developmental_stage_id"].replace("_", ":")))
        anatomical_entity_node = URIRef(get_iri(data["anatomical_entity_id"].replace("_", ":")))
        exp_uri = new_uris['gene_expression_value_base_node']
        anatomical_entity = data['anatomical_entity_id']
        gene_expression_value_node = URIRef(
            f"{exp_uri}/{id_number}/{source_idx}_{anatomical_entity}"
        )

        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), gene_expression_value_node))
        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), anatomical_entity_node))
        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), developmental_stage_node))

        # g.add((gene_node, URIRef(PREDICATES['sio_is_part_of']), cellular_component))
        # g.add((gene_node, URIRef(PREDICATES['sio_is_represented_by']), gene_symbol))

        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_has_value"]),
                Literal(data["expression_level"], datatype=XSD.double),
            )
        )

        g.add((anatomical_entity_node, RDF.type, URIRef(NODE_TYPES["anatomical_entity_node"])))
        g.add(
            (
                anatomical_entity_node,
                RDFS.label,
                Literal(data["anatomical_entity_name"], datatype=XSD.string),
            )
        )
        g.add(
            (
                developmental_stage_node,
                RDFS.label,
                Literal(data["developmental_stage_name"], datatype=XSD.string),
            )
        )
        g.add((anatomical_entity_node, SKOS.exactMatch, anatomical_entity_node))

        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add((gene_expression_value_node, URIRef(PREDICATES["sio_has_input"]), gene_node))
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_has_output"]),
                gene_expression_value_node,
            )
        )

        # Add experimental node
        for exp in experimental_process_data:
            experimental_process_node = add_experimental_process_node(
                g=g,
                id_number=id_number,
                source_idx=source_idx,
                data=exp,
                new_uris=new_uris,
            )
            if experimental_process_node:
                g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), experimental_process_node))
                g.add((experimental_process_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
        ## Input gene, anatomical entity
        if experimental_process_node:
            g.add((experimental_process_node, URIRef(PREDICATES["sio_has_input"]), gene_node))
            g.add(
                (
                    experimental_process_node,
                    URIRef(PREDICATES["sio_has_input"]),
                    anatomical_entity_node,
                )
            )


def add_tested_substance_node(
    g: Graph,
    inchi: str,
    smiles: str,
    compound_id: str,
    compound_name: str,
):
    """Create and add a tested substance node to the RDF graph.

    :param g: RDF graph to which the tested substance node will be added.
    :param inchi: InChI identifier of the compound.
    :param smiles: SMILES identifier of the compound.
    :param compound_id: Compound ID.
    :param compound_name: Name of the compound.
    :return: URIRef for the created tested substance node.
    """
    uri_cas = "https://pubchem.ncbi.nlm.nih.gov/compound/" + str(compound_id).strip("CID")
    tested_substance_node = URIRef((uri_cas))
    g.add((tested_substance_node, RDF.type, URIRef(NODE_TYPES["tested_substance_node"])))
    g.add((tested_substance_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))
    g.add(
        (
            tested_substance_node,
            URIRef(PREDICATES["chebi_inchi"]),
            Literal(inchi, datatype=XSD.string),
        )
    )
    g.add(
        (
            tested_substance_node,
            URIRef(PREDICATES["chebi_smiles"]),
            Literal(smiles, datatype=XSD.string),
        )
    )
    g.add(
        (
            tested_substance_node,
            URIRef(PREDICATES["cheminf_compound_id"]),
            Literal(compound_id, datatype=XSD.string),
        )
    )
    return tested_substance_node


def add_experimental_process_node(
    g: Graph,
    id_number: str,
    source_idx: str,
    data: list,
    new_uris: dict,
):
    """Create and add an experimental process node to the RDF graph.

    :param g: RDF graph to which the experimental process node will be added.
    :param id_number: Unique identifier for the experimental process node.
    :param source_idx: Source index for the experimental process.
    :param data: List of dictionaries containing experimental process information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :return: URIRef for the created experimental process node.
    """
    pubchem_assay_id = data["pubchem_assay_id"]
    if pubchem_assay_id is not None:
        pubchem_assay_iri = "https://pubchem.ncbi.nlm.nih.gov/bioassay/" + str(
            pubchem_assay_id
        ).strip("AID")
        assay_type = data["assay_type"]
        inchi = data["InChI"]
        outcome = data["outcome"]
        compound_cid = data["compound_cid"]
        compound_name = data["compound_name"]
        smiles = data["SMILES"]
        experimental_process_node = URIRef(
            (f"{new_uris['experimental_process_node']}/{id_number}/{source_idx}_{pubchem_assay_id}")
        )
        g.add(
            (experimental_process_node, RDF.type, URIRef(NODE_TYPES["experimental_process_node"]))
        )
        g.add((experimental_process_node, RDF.type, URIRef(pubchem_assay_iri)))
        g.add((experimental_process_node, RDFS.label, Literal(assay_type, datatype=XSD.string)))
        g.add(
            (
                experimental_process_node,
                DC.identifier,
                Literal(pubchem_assay_id, datatype=XSD.string),
            )
        )
        # has tested substance
        tested_substance_node = add_tested_substance_node(
            g=g, inchi=inchi, smiles=smiles, compound_id=compound_cid, compound_name=compound_name
        )
        g.add(
            (experimental_process_node, URIRef(PREDICATES["sio_has_input"]), tested_substance_node)
        )
        # has outcome
        g.add((experimental_process_node, URIRef(PREDICATES["sio_has_output"]), Literal(outcome)))
        return experimental_process_node


def add_pathway_node(g: Graph, data: list, source: str):
    """Create and add a pathway node to the RDF graph.

    :param g: RDF graph to which the pathway node will be added.
    :param data: Dictionary containing pathway information.
    :param source: Source of the pathway information (e.g., WikiPathways, Reactome).
    :return: URIRef for the created pathway node.
    """
    pathway_label = data["pathway_label"]
    pathway_id = data["pathway_id"]
    if pathway_id is not None:
        if source == "WikiPathways":
            pathway_iri = "https://www.wikipathways.org/pathways/" + str(pathway_id)
            iri_source = "https://www.wikipathways.org/"

        if source == "OpenTargets_reactome":
            pathway_iri = "https://reactome.org/content/detail/" + str(pathway_id)
            iri_source = "https://reactome.org/"
            source = "Reactome"
        if source == "MINERVA":
            pathway_iri = "https://minerva-net.lcsb.uni.lu/api/" + str(pathway_id)
            iri_source = "https://minerva-net.lcsb.uni.lu/"
        pathway_node = URIRef(pathway_iri)
        g.add((pathway_node, RDF.type, URIRef(NODE_TYPES["pathway_node"])))
        g.add((pathway_node, RDFS.label, Literal(pathway_label, datatype=XSD.string)))
        g.add((pathway_node, URIRef(PREDICATES["sio_has_source"]), URIRef(iri_source)))
        g.add((URIRef(iri_source), RDFS.label, Literal(source, datatype=XSD.string)))
        return pathway_node
    else:
        return None


def add_compound_node(g: Graph, compound: dict, gene_node: URIRef) -> URIRef:
    """Create and add a compound node to the RDF graph.

    :param g: RDF graph to which the compound node will be added.
    :param compound: Dictionary containing compound information.
    :param gene_node: URIRef of the gene node associated with the compound.
    :return: URIRef for the created compound node.
    """
    chembl_id = compound["chembl_id"]
    drugbank_id = compound["drugbank_id"]
    compound_name = compound["compound_name"]
    compound_cid = compound["compound_cid"]
    is_approved = compound["is_approved"]
    clincal_trial_phase = compound["clincal_trial_phase"]
    relation = compound["relation"]
    compound_node = URIRef(f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}")
    g.add((compound_node, RDFS.label, Literal(str(compound_name), datatype=XSD.string)))
    g.add((compound_node, RDF.type, URIRef(NODE_TYPES["tested_substance_node"])))
    if is_approved:
        g.add(
            (
                compound_node,
                URIRef(PREDICATES["sio_has_value"]),
                URIRef("http://purl.obolibrary.org/obo/NCIT_C172573"),
            )
        )
    if relation in MOAS:
        relation_iri = MOAS[relation]
        g.add((compound_node, URIRef(relation_iri), gene_node))
    if clincal_trial_phase:
        clinical_phase_iri = CLINICAL_PHASES[str(clincal_trial_phase)]
        g.add(
            (
                compound_node,
                URIRef("http://purl.obolibrary.org/obo/PATO_0000083"),
                URIRef(clinical_phase_iri),
            )
        )
    for source_id in [drugbank_id, compound_cid]:
        if source_id:
            iri = None
            if source_id == drugbank_id:
                iri = f"https://go.drugbank.com/drugs/{source_id}"
            if source_id == compound_cid:
                iri = "https://pubchem.ncbi.nlm.nih.gov/compound/" + str(compound_cid)
            if iri:
                id_node = URIRef(iri)
                g.add((compound_node, SKOS.exactMatch, id_node))
                g.add((id_node, SKOS.exactMatch, compound_node))
                g.add((id_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))
                g.add((id_node, RDF.type, URIRef(NODE_TYPES["tested_substance_node"])))
            if isinstance(compound["adverse_effect"], list):
                for i in compound["adverse_effect"]:
                    ae = i["name"]
                    ae_node = add_ae_node(g, ae)
                    if ae_node:
                        g.add(((compound_node, URIRef(PREDICATES["has_adverse_event"]), ae_node)))
    return compound_node


def add_ae_node(g: Graph, ae: str) -> URIRef:
    """Create and add an adverse event node to the RDF graph.

    :param g: RDF graph to which the adverse event node will be added.
    :param ae: Name of the adverse event.
    :return: URIRef for the created adverse event node.
    """
    if ae:
        ae_f = str(ae).replace(" ", "_")
        ae_node = URIRef(f"https://biodatafuse.org/rdf/ae/{ae_f}")
        g.add((ae_node, RDFS.label, Literal(ae, datatype=XSD.string)))
        g.add((ae_node, RDF.type, URIRef(NODE_TYPES["adverse_event_node"])))
        return ae_node
    else:
        return None


def add_go_cpf(g: Graph, process_data: dict) -> URIRef:
    """Create and add a Gene Ontology (GO) node to the RDF graph.

    :param g: RDF graph to which the GO node will be added.
    :param process_data: Dictionary containing Gene Ontology information.
    :return: URIRef for the created GO node.
    """
    curie = process_data["go_id"]
    if curie:
        iri = get_iri(curie)
        label = process_data["go_name"]
        go_type = (
            URIRef(GO_TYPES[process_data["go_type"]])
            if process_data["go_type"] in GO_TYPES
            else None
        )
        go_cpf = URIRef(iri)
        g.add((go_cpf, RDFS.label, Literal(label, datatype=XSD.string)))
        if go_type:
            g.add((go_cpf, RDF.type, go_type))
        return go_cpf
    else:
        return None


def add_metadata(graph, graph_uri, version_iri=None, author=None, orcid=None):
    """
    Add metadata to the RDF graph, including creation date, version, author, and ORCID.

    Args:
        graph (rdflib.Graph): The RDF graph to which metadata is added.
        graph_uri (str): URI identifying the RDF graph.
        version_iri (str, optional): Version IRI to add.
        author (str, optional): Author's name.
        orcid (str, optional): Author's ORCID.

    Returns:
        None
    """
    # Automatically get the current date and time in ISO 8601 format
    creation_date = datetime.utcnow().isoformat() + "Z"

    # Define the URI for the graph resource
    graph_resource = URIRef(graph_uri)

    # Add the creation date to the graph resource
    graph.add((graph_resource, DCTERMS.created, Literal(creation_date, datatype=XSD.dateTime)))

    if version_iri:
        graph.add((graph_resource, DCTERMS.identifier, URIRef(version_iri)))

    if author:
        graph.add((graph_resource, DCTERMS.creator, Literal(author)))
        if orcid:
            graph.add((graph_resource, URIRef("https://w3id.org/reproduceme#ORCID"), URIRef(orcid)))


def generate_rdf(
    df: pd.DataFrame, base_uri: str, version_iri: str, author: str, orcid: str
) -> Graph:
    """Generate an RDF graph from the provided DataFrame.

    :param df: DataFrame containing the data to be converted to RDF.
    :param base_uri: Base URI to use for the nodes in the graph.
    :param version_iri: Version IRI to use as the graph URI and add metadata.
    :param author: Author's name.
    :param orcid: Author's ORCID.
    :return: RDF graph constructed from the DataFrame.
    """
    g = Graph()

    # Load ontology or base RDF
    with resources.path("pyBiodatafuse.resources", "base.owl") as owl_path:
        g.parse(owl_path)

    # Define the base URI and namespace bindings
    new_uris = URIS
    for key, value in new_uris.items():
        new_value = base_uri + value
        new_uris[key] = new_value
        g.bind(key, new_value)

    for key, value in NAMESPACE_BINDINGS.items():
        g.bind(key, value)

    df = df.applymap(replace_na_none)

    for i, row in df.iterrows():
        # Unpack the relevant columns
        source_idx = row["identifier"]
        source_namespace = row["identifier.source"]
        target_idx = row["target"]
        target_namespace = row["target.source"]
        expression_data = row[BGEE_GENE_EXPRESSION_LEVELS_COL]
        experimental_process_data = row["PubChem_assays"]
        disease_data = []
        for source in [DISGENET_DISEASE_COL, "OpenTargets_diseases"]:
            source_el = row[source]
            if isinstance(source_el, list):
                disease_data += source_el

        if (
            pd.isna(source_idx)
            or pd.isna(source_namespace)
            or pd.isna(target_idx)
            or pd.isna(target_namespace)
        ):
            continue

        source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
        target_curie = normalize_curie(f"{target_namespace}:{target_idx}")

        if not source_curie or not target_curie:
            continue

        id_number = f"{i:06d}"
        gene_node = add_gene_node(g, row)

        if len(disease_data) > 0:
            for j, disease in enumerate(disease_data):
                add_gene_disease_associations(
                    g, id_number, source_idx, gene_node, disease, new_uris, j
                )

        add_gene_expression_data(
            g,
            id_number,
            source_idx,
            gene_node,
            expression_data,
            experimental_process_data,
            new_uris,
        )

        for source in ["WikiPathways", "MINERVA", "OpenTargets_reactome"]:
            for pathway_data in row[source]:
                if pathway_data["pathway_id"]:
                    pathway_node = add_pathway_node(g=g, data=pathway_data, source=source)
                    g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node))
                    g.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))

        processes_data = row["OpenTargets_go"]
        for process_data in processes_data:
            go_cpf = add_go_cpf(g, process_data)
            if go_cpf:
                g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), go_cpf))
                g.add((go_cpf, URIRef(PREDICATES["sio_has_part"]), gene_node))

        compound_data = row["OpenTargets_compounds"]
        for compound in compound_data:
            add_compound_node(g, compound, gene_node)

    # Add metadata to the RDF graph
    add_metadata(g, version_iri, version_iri, author, orcid)

    return g
