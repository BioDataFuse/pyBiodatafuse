# coding: utf-8

"""Python module to produce an RDF graph from the property table."""

import logging
from datetime import datetime
from importlib import resources

import numpy as np
import pandas as pd
from bioregistry import get_iri, normalize_curie
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import DC, DCTERMS, OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import (
    BGEE_GENE_EXPRESSION_LEVELS_COL,
    CLINICAL_PHASES,
    DATA_SOURCES,
    DISGENET_DISEASE_COL,
    GO_TYPES,
    DISEASE_IDENTIFIER_TYPES,
    MOAS,
    NAMESPACE_BINDINGS,
    NODE_TYPES,
    OPENTARGETS_DISEASE_COL,
    PREDICATES,
    URIS,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def replace_na_none(item):
    """
    Recursively replaces 'na', 'nan', 'none' strings, Python None, and NaN values with None.

    :param item: The value to process (string, None, NaN, list, dict, or np.ndarray).
    :return: The processed item with 'na', 'nan', 'none', None, and NaN replaced by None.
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


def add_gene_protein_nodes(g: Graph, row) -> URIRef:
    """Create and add a gene node to the RDF graph.

    :param g: RDF graph to which the gene node will be added.
    :param row: DataFrame row containing gene information.
    :return: List of URIRefs for the created gene node and protein node
    """
    target = row["target"]
    if target:
        gene_node = URIRef(f"http://identifiers.org/ensembl#{target}")
        protein_node = URIRef(f"http://identifiers.org/ensembl#{target}xProtein")
        # TODO replace by UniProt ID when available
        g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((gene_node, RDFS.label, Literal(row["identifier"], datatype=XSD.string)))
        g.add((protein_node, URIRef(PREDICATES["has_gene_template"]), gene_node))
        g.add((protein_node, RDF.type, URIRef(NODE_TYPES['protein_node'])))
        g.add((protein_node, RDFS.label, Literal(f"{target}xProtein", datatype=XSD.string)))
        return [gene_node, protein_node]
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
    disease_data.get()
    disease_node = URIRef(disease_iri)
    g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
    g.add(
        (disease_node, RDFS.label, Literal(disease_data.get("disease_name"), datatype=XSD.string))
    )
    for identifier_type in DISEASE_IDENTIFIER_TYPES:
        curie_field = disease_data.get(identifier_type, None)
        if curie_field and 'ncit' not in curie_field:
            curies = [curie_field]
            if "," in (curie_field):
                curies = [i for i in curie_field.split(", ")]
            for curie in curies:
                disease_source_iri = get_iri(curie)
                if disease_source_iri is None:
                    if ":" in curie:
                        curie = curie.split(":")[1]
                    disease_source_iri = get_iri("obo:" + curie)
                    g.add(
                        (disease_node, OWL.SameAs, URIRef(disease_source_iri))
                    )  # Some of the data does not look like a skos:exactMatch
                else:
                    g.add((disease_node, OWL.sameAs, URIRef(disease_source_iri)))
    return disease_node


def add_score_node(
    g: Graph, id_number: str, source_idx: str, 
    disease_id: str, score: float, new_uris: dict, i: int, gene_id: str
) -> URIRef:
    """Create and add a score node for gene-disease associations.

    :param g: RDF graph to which the score node will be added.
    :param id_number: Unique identifier for the score node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the score.
    :param score: Score value for the gene-disease association.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: Index or iteration number used for node URI construction.
    :param gene_id: String value of the gene ID
    :return: URIRef for the created score node.
    """
    score_node = URIRef(f"{new_uris['score_base_node']}/{id_number}{i}{source_idx}_{disease_id}{gene_id}")
    g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
    g.add(
        (
            score_node,
            URIRef(PREDICATES["sio_has_value"]),
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
        f"{new_uris['score_base_node']}/{id_number}{i}{source_idx}_{disease_id}"
    )
    g.add((evidence_idx_node, RDF.type, URIRef(NODE_TYPES["evidence_idx_node"])))
    g.add(
        (
            evidence_idx_node,
            URIRef(PREDICATES["sio_has_value"]),
            Literal(ei, datatype=XSD.double),
        )
    )
    return evidence_idx_node


def add_data_source_node(g: Graph, source: str) -> URIRef:
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :param source: String containing the name of the source of the data
    :return: URIRef for the created data source node.
    """
    data_source_name = Literal(source, datatype=XSD.string)
    data_source_url = URIRef(DATA_SOURCES[source])
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url


def add_gene_disease_associations(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    disease_data: dict,
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
    :param i: the index of the row
    """
    data = disease_data
    gene_id = g.value(subject=gene_node, predicate=RDFS.label)
    disease_node = add_disease_node(g, data)
    if disease_node:
        gene_disease_assoc_node = URIRef(
            f"{new_uris['gene_disease_association']}/{id_number}{i}{source_idx}"
        )
        g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene_disease_association"])))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), gene_node))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), disease_node))
        if data.get("score"):
            disease_umlscui = disease_data.get("disease_umlscui", None)
            if disease_umlscui:
                score_node = add_score_node(
                    g=g,
                    id_number=id_number,
                    source_idx=source_idx,
                    score=data["score"],
                    new_uris=new_uris,
                    i=i,
                    disease_id=disease_umlscui,
                    gene_id=gene_id
                )
                g.add(
                    (
                        gene_disease_assoc_node,
                        URIRef(PREDICATES["sio_has_measurement_value"]),
                        score_node,
                    )
                )
        data_source_node = add_data_source_node(g, "DISGENET")
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
        if not data.get("anatomical_entity_id", None):
            continue

        gene_expression_value_node = URIRef(
            f"{new_uris['gene_expression_value_base_node']}/{id_number}/{source_idx}"
        )
        developmental_stage_node = URIRef(get_iri(data["developmental_stage_id"].replace("_", ":")))
        anatomical_entity_node = URIRef(get_iri(data["anatomical_entity_id"].replace("_", ":")))
        exp_uri = new_uris["gene_expression_value_base_node"]
        anatomical_entity = data["anatomical_entity_id"]
        gene_expression_value_node = URIRef(
            f"{exp_uri}/{id_number}/{source_idx}_{anatomical_entity}"
        )

        g.add(
            (gene_node, URIRef(PREDICATES["sio_has_measurement_value"]), gene_expression_value_node)
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_is_associated_with"]),
                anatomical_entity_node,
            )
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_is_associated_with"]),
                developmental_stage_node,
            )
        )

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
        g.add(
            (
                developmental_stage_node,
                RDF.type,
                URIRef(NODE_TYPES["developmental_stage_node"]),
            )
        )
        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add((gene_expression_value_node, URIRef(PREDICATES["sio_has_input"]), gene_node))
        g.add(
            (
                gene_node,
                URIRef(PREDICATES["sio_has_output"]),
                gene_expression_value_node,
            )
        )

        # Add experimental node
        if experimental_process_data:
            for exp in experimental_process_data:
                experimental_process_node = add_experimental_process_node(
                    g=g,
                    data=exp,
                )
                if experimental_process_node:
                    g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), experimental_process_node))
                    g.add(
                        (experimental_process_node, URIRef(PREDICATES["sio_has_part"]), gene_node)
                        )
                    g.add(
                        (
                            experimental_process_node,
                            URIRef(PREDICATES["sio_has_part"]),
                            anatomical_entity_node,
                        )
                    )
                    data_source_node = add_data_source_node(g, "Bgee")
                    g.add(
                        (
                            gene_expression_value_node,
                            URIRef(PREDICATES["sio_has_source"]),
                            data_source_node,
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
    :param inchi: inchi identifier of the compound.
    :param smiles: smiles identifier of the compound.
    :param compound_id: Compound ID.
    :param compound_name: Name of the compound.
    :return: URIRef for the created tested substance node.
    """
    uri_cas = "https://pubchem.ncbi.nlm.nih.gov/compound/" + str(compound_id).strip("CID")
    tested_substance_node = URIRef((uri_cas))
    g.add((tested_substance_node, RDF.type, URIRef(NODE_TYPES["tested_substance_node"])))
    g.add((tested_substance_node, RDF.type, URIRef(NODE_TYPES["drug_node"])))
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
    data: dict,
):
    """Create and add an experimental process node to the RDF graph.

    :param g: RDF graph to which the experimental process node will be added.
    :param data: List of dictionaries containing experimental process information.
    :return: URIRef for the created experimental process node.
    """
    pubchem_assay_id = data.get("pubchem_assay_id", None)
    if pubchem_assay_id:
        pubchem_assay_iri = "https://pubchem.ncbi.nlm.nih.gov/bioassay/" + str(
            pubchem_assay_id
        ).strip("AID")
        assay_type = data["assay_type"]
        inchi = data["inchi"]
        outcome = data["outcome"]
        compound_cid = data["compound_cid"]
        compound_name = data["compound_name"]
        smiles = data["smiles"]
        experimental_process_node = URIRef(pubchem_assay_iri)
        g.add(
            (experimental_process_node, RDF.type, URIRef(NODE_TYPES["experimental_process_node"]))
        )
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
        data_source_node = add_data_source_node(g, "PubChem")
        g.add((experimental_process_node, URIRef(PREDICATES["sio_has_source"]), data_source_node))
        return experimental_process_node


def add_pathway_node(g: Graph, data: dict, source: str):
    """Create and add a pathway node to the RDF graph.

    :param g: RDF graph to which the pathway node will be added.
    :param data: Dictionary containing pathway information.
    :param source: Source of the pathway information (e.g., WikiPathways, Reactome).
    :return: URIRef for the created pathway node.
    """
    pathway_label = data.get("pathway_label", None)
    pathway_id = data.get("pathway_id", None)

    if pathway_id:
        pathway_iri = ""
        iri_source = None
        if source == "WikiPathways":
            pathway_iri = "https://www.wikipathways.org/pathways/" + str(pathway_id.split(":")[1])
            iri_source = "https://www.wikipathways.org/"

        if source == "OpenTargets_reactome":
            pathway_id = pathway_id.split(":")[1]
            pathway_iri = "https://reactome.org/content/detail/" + str(pathway_id)
            iri_source = "https://reactome.org/"
            source = "Reactome"
        if source == "MINERVA":
            pathway_iri = "https://minerva-net.lcsb.uni.lu/api/" + str(pathway_id)
            iri_source = "https://minerva-net.lcsb.uni.lu/"
        pathway_node = URIRef(pathway_iri)
        g.add((pathway_node, RDF.type, URIRef(NODE_TYPES["pathway_node"])))
        g.add((pathway_node, RDFS.label, Literal(pathway_label, datatype=XSD.string)))
        if iri_source:
            g.add((pathway_node, URIRef(PREDICATES["sio_has_source"]), URIRef(iri_source)))
            g.add((URIRef(iri_source), RDFS.label, Literal(source, datatype=XSD.string)))
            # data_source_node = add_data_source_node(g, source)
            g.add((URIRef(iri_source), RDF.type, URIRef(NODE_TYPES["data_source_node"])))
        return pathway_node


def add_drug_node(g: Graph, compound: dict, protein_node: URIRef) -> URIRef:
    """Create and add a compound node to the RDF graph.

    :param g: RDF graph to which the compound node will be added.
    :param compound: Dictionary containing compound information.
    :param protein_node: URIRef of the protein node associated with the compound.
    :return: URIRef for the created compound node.
    """
    chembl_id = compound.get("chembl_id", None)
    drugbank_id = compound.get("drugbank_id", None)
    compound_name = compound.get("compound_name", None)
    compound_cid = compound.get("compound_cid", None)
    is_approved = compound.get("is_approved", None)
    clincal_trial_phase = compound.get("clincal_trial_phase", None)
    relation = compound.get("relation", None)
    if chembl_id:
        chembl_id = chembl_id.split("CHEMBL:")[1]
        drug_node = URIRef(f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}")
        g.add((drug_node, RDFS.label, Literal(str(compound_name), datatype=XSD.string)))
        data_source_node = add_data_source_node(g, "OpenTargets_reactome")
        g.add((drug_node, URIRef(PREDICATES["sio_has_source"]), data_source_node))
        if is_approved:
            g.add(
                (
                    drug_node,
                    RDF.type,  # TODO check
                    URIRef(NODE_TYPES["approved_drug"]),
                )
            )
        else:
            g.add((drug_node, RDF.type, URIRef(NODE_TYPES["tested_substance_node"])))
        if relation in MOAS:
            relation_iri = MOAS[relation]
            g.add((drug_node, URIRef(relation_iri), protein_node))
        if clincal_trial_phase:
            clinical_phase_iri = CLINICAL_PHASES[str(clincal_trial_phase)]
            g.add(
                (
                    drug_node,
                    URIRef(PREDICATES["sio_has_value"]),
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
                    g.add((drug_node, OWL.sameAs, id_node))
                    g.add((id_node, OWL.sameAs, drug_node))
                    g.add((id_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))
                    g.add((id_node, RDF.type, URIRef(NODE_TYPES["tested_substance_node"])))
                ae_list = compound.get("adverse_effect", None)
                if ae_list:
                    for i in compound["adverse_effect"]:
                        ae = i.get("name", None)
                        if ae:
                            ae_node = add_ae_node(g, ae)
                            if ae_node:
                                g.add((ae_node, URIRef(PREDICATES["is_preceded_by"]), drug_node))
                                g.add((drug_node, URIRef(PREDICATES["precedes"]), ae_node))
        return drug_node


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
            g.add((go_cpf, RDFS.subClassOf, go_type))
        data_source_node = add_data_source_node(g, "OpenTargets_reactome")
        g.add((go_cpf, URIRef(PREDICATES["sio_has_source"]), data_source_node))
        return go_cpf
    else:
        return None


def add_metadata(
    g: Graph, graph_uri: str, metadata: dict, version_iri=None, author=None, orcid=None
):
    """
    Add metadata to the RDF graph, including creation date, version, author, and ORCID.

    :param g: The RDF graph to which metadata is added.
    :param graph_uri: URI identifying the RDF graph.
    :param metadata: Combined metadata for a BioDatafuse query.
    :param version_iri: Version IRI to add.
    :param author: Author's name.
    :param orcid: Author's ORCID.
    """
    # Automatically get the current date and time in ISO 8601 format
    creation_date = datetime.utcnow().isoformat() + "Z"

    # Define the URI for the graph resource
    graph_resource = URIRef(graph_uri)

    # Add the creation date to the graph resource
    g.add((graph_resource, DCTERMS.created, Literal(creation_date, datatype=XSD.dateTime)))

    if version_iri:
        g.add((graph_resource, DCTERMS.identifier, URIRef(version_iri)))

    if orcid:
        orcid_node = URIRef(orcid)
        g.add((orcid_node, OWL.sameAs, Literal(author, datatype=XSD.string)))
        g.add((graph_resource, DCTERMS.creator, orcid_node))
    # Metadata file: add version, api version, and dates to sources
    api_version = None
    version = None
    for source in [
        "DISGENET",
        "Open Targets GraphQL & REST API Beta",
        "MINERVA",
        "WikiPathways",
        "StringDB",
        "BridgeDb",
    ]:
        entries = [
            item for item in metadata if isinstance(item, dict) and item.get("datasource") == source
        ]
        for entry in entries:
            data_version = None
            year = None
            month = None
            version = None
            api_version = None
            # date = None
            url_service = None
            source_node = None
            # input_type = entry.get("query").get("input_type")
            # date = entry.get("query").get("date")
            query = entry.get("query", None)
            if query:
                url_service = query.get("url", None)
            source = entry.get("source", None)
            if source == "Open Targets GraphQL & REST API Beta":
                source_node = URIRef(DATA_SOURCES["OpenTargets_reactome"])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("The Open Targets Platform", datatype=XSD.string),
                    )
                )
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    data_version = metadata_f.get("data_version", None)
                    if data_version:
                        year = data_version.get("year")
                        month = data_version.get("month")
                        version = f"{year}-{month}-01"
                source_version = metadata_f.get("source_version") if metadata_f else None
                api_version_data = source_version.get("apiVersion") if source_version else None
                if api_version_data:
                    api_version = ".".join([str(i) for i in api_version_data.values()])
                else:
                    api_version = None
            elif source == "DISGENET":
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("DisGeNET", datatype=XSD.string),
                    )
                )
                if metadata_f:
                    data_version = metadata_f.get("version", None)
                    version = metadata_f.get("version", None)

            elif source == "MINERVA":
                source_node = URIRef(DATA_SOURCES[source])
                # project = entry.get('query').get('MINERVA project')
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    version = metadata_f.get("source_version")
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("The MINERVA Platform", datatype=XSD.string),
                    )
                )

            elif source == "WikiPathways":
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    version = metadata_f.get("source_version", None)
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("WikiPathways", datatype=XSD.string),
                    )
                )

            elif source == "BridgeDb":  # Several metadata fields left unRDFied
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("BridgeDb", datatype=XSD.string),
                    )
                )
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    version = metadata_f.get("bridgedb.version", None)
                    api_version = metadata_f.get("webservice.version", None)

            elif source == "StringDB":
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("STRING", datatype=XSD.string),
                    )
                )

            # Add api node for source
            if url_service:
                api_node = URIRef(url_service)
                g.add((api_node, RDF.type, URIRef("https://schema.org/WebAPI")))
            if source_node:
                g.add((api_node, URIRef("https://schema.org/provider"), source_node))
                g.add((source_node, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
                if api_version:
                    g.add(
                        (
                            api_node,
                            URIRef("http://www.geneontology.org/formats/oboInOwl#version"),
                            Literal(api_version),
                        )
                    )
            if version:
                g.add(
                    (
                        source_node,
                        URIRef("http://www.geneontology.org/formats/oboInOwl#version"),
                        Literal(version),
                    )
                )


def add_transporter_inhibitor_node(g: Graph, transporter_inhibitor_data: dict, base_uri: str):
    """Add a transporter inhibitor node.

    :param g: RDFLib graph
    :param transporter_inhibitor_data: dictionary with the membrane-compound interaction data
    :param base_uri: The project base uri
    """
    data = transporter_inhibitor_data
    compound_name = data.get("compound_name", None)
    inchikey = data.get("inchikey", None)
    smiles = data.get("smiles", None)
    compound_cid = data.get("compound_cid", None)
    molmedb_id = data.get("molmedb_id", None)
    # source_pmid = data.get("source_pmid", None)
    chebi_id = data.get("chebi_id", None)
    drugbank_id = data.get("drugbank_id", None)
    uniprot_trembl_id = data.get("uniprot_trembl_id", None)
    if compound_cid:
        drug_node = URIRef(f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_cid}")
        g.add((drug_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))
        g.add(
            (
                drug_node,
                URIRef(PREDICATES["chebi_inchi"]),
                Literal(inchikey, datatype=XSD.string),
            )
        )
        g.add((drug_node, URIRef(PREDICATES["chebi_smiles"]), Literal(smiles)))
        g.add(
            (
                drug_node,
                OWL.sameAs,
                URIRef(f"https://molmedb.upol.cz/mol/{molmedb_id}"),
            )
        )
        g.add(
            (
               drug_node,
               OWL.sameAs,
               URIRef(f"https://identifiers.org/CHEBI#{chebi_id}")
            )
        )
        g.add(
            (
                drug_node,
                OWL.sameAs,
                URIRef(f"https://www.drugbank.ca/drugs/{drugbank_id}"),
            )
        )
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{uniprot_trembl_id}"),
                URIRef(PREDICATES["sio_has_part"]),
                drug_node,
            )
        )
        g.add(
            (
                URIRef(f"https://identifiers.org/CHEBI#{chebi_id}"),
                RDF.type,
                URIRef("https://purl.obolibrary.org/CHEMINF_000407"),
            )
        )
        g.add(
            (
                drug_node,
                OWL.sameAs,
                URIRef(f"https://identifiers.org/CHEBI#{chebi_id}"),
            )
        )
        g.add(
            (
                URIRef(f"https://www.drugbank.ca/drugs/{drugbank_id}"),
                RDF.type,
                URIRef("https://purl.obolibrary.org/CHEMINF_000406"),
            )
        )
        g.add(
            (
                drug_node,
                OWL.sameAs,
                URIRef(f"https://www.drugbank.ca/drugs/{drugbank_id}"),
            )
        )
        inhibition_node = URIRef(base_uri + f"inhibition/{uniprot_trembl_id}_{compound_cid}")
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{uniprot_trembl_id}"),
                URIRef(PREDICATES["sio_is_part_of"]),
                inhibition_node,
            )
        )
        g.add(
            (
                drug_node,
                URIRef(PREDICATES["sio_is_part_of"]),
                inhibition_node,
            )
        )
        g.add(
            (
                inhibition_node,
                RDF.type,
                URIRef("https://purl.obolibrary.org/GO_0032410"),
            )
        )
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{uniprot_trembl_id}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )


def add_ppi_data(
    g: Graph, protein_node: URIRef, protein_name: str, entry: dict, base_uri: str, new_uris: dict
) -> URIRef:
    """Add a protein protein interaction node.

    :param g: RDFLib graph
    :param protein_node: URIRef for the target protein
    :param protein_name: the name of the target protein
    :param entry: the ppi dictionary
    :param base_uri: the base URI for the project
    :param new_uris: dictionary with project node URIs
    :return: a ppi URIRef node
    """
    stringdb_link_to = entry.get("stringdb_link_to", None)
    ensembl = entry.get("Ensembl", None).split(":")[1]
    score = entry.get("score", None)
    if score:
        score = float(score)
        # Nodes
        ppi_node = URIRef(base_uri + f"ppi/{protein_name}_{ensembl}")
        g.add(
            (
                ppi_node,
                RDF.type,
                URIRef(NODE_TYPES["ppi_node"]),
            )
        )
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )
        g.add(
            (
                ppi_node,
                URIRef(PREDICATES["sio_has_part"]),
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
            )
        )
        g.add(
            (
                ppi_node,
                URIRef(PREDICATES["sio_has_part"]),
                protein_node,
            )
        )
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
                URIRef(PREDICATES["sio_is_part_of"]),
                ppi_node,
            )
        )
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )
        g.add(
            (
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
                RDFS.label,
                Literal(stringdb_link_to, datatype=XSD.string),
            )
        )
        g.add(
            (
                URIRef(f"http://identifiers.org/ensembl#{ensembl}"),
                OWL.sameAs,
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
            )
        )
        g.add(
            (
                URIRef(f"http://identifiers.org/ensembl#{ensembl}"),
                OWL.sameAs,
                URIRef(f"https://www.uniprot.org/uniprotkb/{stringdb_link_to}"),
            )
        )
        g.add(
            (
                URIRef(f"http://identifiers.org/ensembl#{ensembl}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )
        score_node = URIRef(f"{new_uris['score_base_node']}/{protein_name}_{ensembl}")
        g.add((ppi_node, URIRef(PREDICATES["sio_has_measurement_value"]), score_node))
        g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
        g.add(
            (
                score_node,
                URIRef(PREDICATES["sio_has_value"]),
                Literal(score, datatype=XSD.double),
            )
        )
        return ppi_node


def add_literature_based_data(
    g: Graph,
    entry: dict,
    gene_node: URIRef
):
    """Add a literature based node.

    :param g: RDFLib graph
    :param entry: the literature based data dictionary
    :param gene_node: The gene URIRef
    """
    source = entry.get("source", None)
    if source:
        if "PMID" in source:
            source = source.split(": ")[1]
            identifier = entry["id"]
            label = entry["disease_name"]
            source_url = f"https://pubmed.ncbi.nlm.nih.gov/{source}"
            g.add((URIRef(source_url), RDF.type, URIRef(NODE_TYPES["article"])))
            g.add(
                (
                    URIRef(source_url),
                    URIRef(PREDICATES["sio_refers_to"]),
                    URIRef(f"https://biodatafuse.org/identifiers/{identifier}"),
                )
            )
            g.add(
                (
                    URIRef(source_url),
                    URIRef(PREDICATES["sio_refers_to"]),
                    gene_node,
                )
            )
            g.add(
                (
                    URIRef(f"https://biodatafuse.org/identifiers/{identifier}"),
                    RDFS.label,
                    Literal(label, datatype=XSD.string),
                )
            )
            g.add(
                (
                    URIRef(f"https://biodatafuse.org/identifiers/{identifier}"),
                    RDF.type,
                    URIRef(NODE_TYPES["disease_node"]),
                )
            )


def generate_rdf(
    df: pd.DataFrame, base_uri: str, version_iri: str, author: str, orcid: str, metadata: dict, open_only: bool, load_ontology: bool
) -> Graph:
    """Generate an RDF graph from the provided DataFrame.

    :param df: DataFrame containing the data to be converted to RDF.
    :param base_uri: Base URI to use for the nodes in the graph.
    :param version_iri: Version IRI to use as the graph URI and add metadata.
    :param author: Author's name.
    :param orcid: Author's ORCID.
    :param metadata: Combined metadata for a BioDatafuse query.
    :return: RDF graph constructed from the DataFrame.
    """
    g = Graph()
    # Define the base URI and namespace bindings
    new_uris = {key: base_uri + value for key, value in URIS.items()}

    # Update g.bind only after constructing new URIs
    for key, new_value in new_uris.items():
        g.bind(key, new_value)

    for key, value in NAMESPACE_BINDINGS.items():
        g.bind(key, value)

    df = df.applymap(replace_na_none)

    for i, row in df.iterrows():
        # Unpack the relevant columns
        source_idx = row.get("identifier", None)
        source_namespace = row.get("identifier.source", None)
        target_idx = row.get("target", None)
        target_namespace = row.get("target.source", None)
        expression_data = row.get(BGEE_GENE_EXPRESSION_LEVELS_COL, None)
        experimental_process_data = row.get("PubChem_assays", None)
        processes_data = row.get("OpenTargets_go", None)
        compound_data = row.get("OpenTargets_gene_compounds", None)
        literature_based_data = row.get("literature_based_info", None)
        transporter_inhibitor_data = row.get("MolMeDB_transporter_inhibitor", None)
        stringdb_data = row.get("StringDB_ppi", None)
        # molmedb_data = row.get("MolMeDB_transporter_inhibitor", None)
        disease_data = []
        for source in [DISGENET_DISEASE_COL, OPENTARGETS_DISEASE_COL]:
            if open_only and source == DISGENET_DISEASE_COL:  # TODO implement open data only feature properly
                continue
            source_el = row.get(source)
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
        gene_node, protein_node = add_gene_protein_nodes(g, row)

        if len(disease_data) > 0:
            for j, disease in enumerate(disease_data):
                add_gene_disease_associations(
                    g, id_number, source_idx, gene_node, disease, new_uris, j
                )
        if expression_data:
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
            if row.get(source, None):
                for pathway_data in row[source]:
                    if pathway_data["pathway_id"]:
                        pathway_node = add_pathway_node(g=g, data=pathway_data, source=source)
                        g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node))
                        g.add((protein_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node))
                        g.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
                        g.add(
                            (
                                pathway_node,
                                URIRef(PREDICATES["sio_has_source"]),
                                URIRef(DATA_SOURCES[source]),
                            )
                        )

        if processes_data:
            for process_data in processes_data:
                go_cpf = add_go_cpf(g, process_data)
                if go_cpf:
                    g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), go_cpf))
                    g.add((go_cpf, URIRef(PREDICATES["sio_has_part"]), gene_node))

        if compound_data:
            for compound in compound_data:
                add_drug_node(g, compound, protein_node)

        if literature_based_data:
            if isinstance(literature_based_data, list):
                entries = literature_based_data
                for entry in entries:
                    add_literature_based_data(g, entry, gene_node)
            elif isinstance(literature_based_data, dict):
                entry = literature_based_data
                add_literature_based_data(g, entry, gene_node)

        if transporter_inhibitor_data:
            for entry in transporter_inhibitor_data:
                add_transporter_inhibitor_node(g, entry, base_uri)
        if stringdb_data and isinstance(stringdb_data, list):
            protein_name = f"{row.target}_protein"
            protein_name = f"{row.target}_protein"
            for entry in stringdb_data:
                if entry.get("Ensembl", None):
                    add_ppi_data(g, protein_node, protein_name, entry, base_uri, new_uris)

    # Add metadata to the RDF graph
    if metadata:
        add_metadata(
            g=g,
            version_iri=version_iri,
            author=author,
            orcid=orcid,
            metadata=metadata,
            graph_uri=version_iri,
        )
        # Load ontology or base RDF
    if load_ontology:
        with resources.path("pyBiodatafuse.resources", "biodatafuse.owl") as owl_path:
            temp_graph = Graph()  # Temporary graph for loading ontology
            temp_graph.parse(owl_path)
            # Iterate through classes and predicates in the ontology
            for subj, pred, obj in temp_graph:
                # Check if each element individually exists in g as part of any triple
                subj_exists = any(g.triples((subj, None, None)))
                pred_exists = any(g.triples((None, pred, None)))
                obj_exists = any(g.triples((None, None, obj)))
                # Add to graph only if any of these elements exist in g
                if subj_exists or pred_exists or obj_exists:
                    g.add((subj, pred, obj))
    return g
