# coding: utf-8

"""Python module to export a NetworkX graph to neo4j-compliant format and set the styling for Neo4j."""

import typing
from collections import defaultdict

import networkx as nx
from neo4j import GraphDatabase
from neomodel import (
    BooleanProperty,
    FloatProperty,
    RelationshipFrom,
    RelationshipTo,
    StringProperty,
    StructuredNode,
    StructuredRel,
    config,
    db,
)
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.saver import save_graph_to_graphml


# Predifinding relationship classes
class AssociatedWith(StructuredRel):
    """Relationship between a gene and a disease."""

    datasource = StringProperty()
    score = FloatProperty()
    ei = FloatProperty()


class PartOf(StructuredRel):
    """Relationship between a gene and a pathway."""

    datasource = StringProperty()


class InteractsWith(StructuredRel):
    """Relationship between a gene and a gene."""

    datasource = StringProperty()
    score = FloatProperty()


class Activates(StructuredRel):
    """Relationship between a compound and a gene."""

    datasource = StringProperty()


class HasSideEffect(StructuredRel):
    """Relationship between a compound and a side effect."""

    datasource = StringProperty()


class Treats(StructuredRel):
    """Relationship between a compound and a disease."""

    datasource = StringProperty()


class Inhibits(StructuredRel):
    """Relationship between a compound and a gene."""

    datasource = StringProperty()


class ExpressedBy(StructuredRel):
    """Relationship between a gene and a tissue."""

    datasource = StringProperty()
    expression_level = FloatProperty()
    developmental_stage = StringProperty()
    developmental_stage_id = StringProperty()
    confidence_level = StringProperty()
    confidence_level_id = StringProperty()


class UpstreamOf(StructuredRel):
    """Relationship between a molecular initiating event and an adverse outcome pathway."""

    datasource = StringProperty()


class DownstreamOf(StructuredRel):
    """Relationship between a key event and a key event."""

    datasource = StringProperty()


# Defining the nodes
class Gene(StructuredNode):
    """Gene node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    Ensembl = StringProperty()
    node_type = Cons.GENE_NODE_LABEL

    # outgoing relations
    interacts_with = RelationshipTo(
        Cons.GENE_NODE_LABEL, Cons.INTERACTS_WITH, model=InteractsWith
    )  # PPI
    part_of_pathway = RelationshipTo(
        Cons.PATHWAY_NODE_LABEL, Cons.PART_OF, model=PartOf
    )  # Gene -> Pathway
    part_of_biological_process = RelationshipTo(
        Cons.GO_BP_NODE_LABEL.replace(" ", ""), Cons.PART_OF, model=PartOf
    )  # Gene -> BiologicalProcess
    part_of_molecular_function = RelationshipTo(
        Cons.GO_MF_NODE_LABEL.replace(" ", ""), Cons.PART_OF, model=PartOf
    )  # Gene -> MolecularFunction
    part_of_cellular_component = RelationshipTo(
        Cons.GO_CC_NODE_LABEL.replace(" ", ""), Cons.PART_OF, model=PartOf
    )  # Gene -> CellularComponent
    associated_with_disease = RelationshipTo(
        Cons.DISEASE_NODE_LABEL, Cons.ASSOCIATED_WITH, model=AssociatedWith
    )  # Gene -> Disease
    expressed_by = RelationshipTo(
        Cons.ANATOMICAL_NODE_LABEL.replace(" ", ""), Cons.EXPRESSED_BY, model=ExpressedBy
    )  # Gene -> Tissue
    associated_with_aop = RelationshipTo(
        Cons.AOP_NODE_LABEL.replace(" ", ""), Cons.ASSOCIATED_WITH, model=AssociatedWith
    )  # Gene -> AO

    # incoming edges
    activates = RelationshipTo(
        Cons.COMPOUND_NODE_LABEL, Cons.ACTIVATES, model=Activates
    )  # Compound -> Gene
    inhibits = RelationshipTo(
        Cons.COMPOUND_NODE_LABEL, Cons.INHIBITS, model=Inhibits
    )  # Compound -> Gene
    # associated_with_compound = RelationshipTo(
    #     Cons.COMPOUND_NODE_LABEL, Cons.ASSOCIATED_WITH, model=AssociatedWith
    # )  # Compound -> Gene


class Disease(StructuredNode):
    """Disease node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    HPO = StringProperty()
    NCI = StringProperty()
    MONDO = StringProperty()
    ORDO = StringProperty()
    MESH = StringProperty()
    UMLS = StringProperty()
    disease_type = StringProperty()
    node_type = Cons.DISEASE_NODE_LABEL

    # incoming relations (Gene -> Disease)
    associated_with = RelationshipFrom(Gene, Cons.ASSOCIATED_WITH, model=AssociatedWith)


class Pathway(StructuredNode):
    """Pathway node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.PATHWAY_NODE_LABEL

    # incoming relations (Gene -> Pathway)
    part_of = RelationshipFrom(Gene, Cons.PART_OF, model=PartOf)


class BiologicalProcess(StructuredNode):
    """Biological process node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.GO_BP_NODE_LABEL

    # incoming relations (Gene -> BiologicalProcess)
    part_of = RelationshipFrom(Gene, Cons.PART_OF, model=PartOf)


class MolecularFunction(StructuredNode):
    """Molecular function node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.GO_MF_NODE_LABEL

    # incoming relations (Gene -> MolecularFunction)
    part_of = RelationshipFrom(Gene, Cons.PART_OF, model=PartOf)


class CellularComponent(StructuredNode):
    """Cellular component node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.GO_CC_NODE_LABEL

    # incoming relations (Gene -> CellularComponent)
    part_of = RelationshipFrom(Gene, Cons.PART_OF, model=PartOf)


class SideEffect(StructuredNode):
    """Side effect node."""

    name = StringProperty(required=True, unique_index=True, unique=True)
    datasource = StringProperty(required=True)
    node_type = Cons.SIDE_EFFECT_NODE_LABEL


class Compound(StructuredNode):
    """Compound node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    chembl_id = StringProperty()
    drugbank_id = StringProperty()
    compound_cid = StringProperty()
    clinical_trial_phase = StringProperty()
    is_approved = BooleanProperty()
    adverse_effect_count = FloatProperty()
    node_type = Cons.COMPOUND_NODE_LABEL

    # outgoing relations (Compound -> Gene, SideEffect, Disease)
    activates = RelationshipTo(Gene, Cons.ACTIVATES, model=Activates)
    has_side_effect = RelationshipTo(SideEffect, Cons.HAS_SIDE_EFFECT, model=HasSideEffect)
    treats = RelationshipTo(Disease, Cons.TREATS, model=Treats)
    inhibits = RelationshipTo(Gene, Cons.INHIBITS, model=Inhibits)

    # incoming relations (Gene -> Compound)
    associated_with = RelationshipFrom(Gene, Cons.ASSOCIATED_WITH, model=AssociatedWith)


class AnatomicalEntity(StructuredNode):
    """Anatomical entity node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.ANATOMICAL_NODE_LABEL

    # incoming relations (Gene -> AnatomicalEntity)
    part_of = RelationshipFrom(Gene, Cons.EXPRESSED_BY, model=ExpressedBy)


class AdverseOutcomePathway(StructuredNode):
    """Adverse outcome pathway node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.AOP_NODE_LABEL

    # incoming relations (AO -> KeyEvent)
    associated_with = RelationshipTo(Gene, Cons.ASSOCIATED_WITH, model=AssociatedWith)


class MolecularInitiatingEvent(StructuredNode):
    """Molecular initiating event node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    node_type = Cons.MIE_NODE_LABEL

    # outgoing relations (MolecularInitiatingEvent -> AOP)
    upstream_of = RelationshipTo(AdverseOutcomePathway, Cons.UPSTREAM_OF, model=UpstreamOf)

    # incoming relations (KeyEvent -> MIE)
    associated_with = RelationshipFrom(
        Cons.KEY_EVENT_NODE_LABEL.replace(" ", ""), Cons.ASSOCIATED_WITH, model=AssociatedWith
    )


class KeyEvent(StructuredNode):
    """Key event node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    organ = StringProperty()
    node_type = Cons.KEY_EVENT_NODE_LABEL

    # outgoing relations (KeyEvent -> MIE)
    upstream_of = RelationshipTo(MolecularInitiatingEvent, Cons.UPSTREAM_OF, model=UpstreamOf)
    associated_with = RelationshipFrom(
        Cons.AO_NODE_LABEL.replace(" ", ""), Cons.ASSOCIATED_WITH, model=AssociatedWith
    )

    # outgoing relations (KeyEvent -> KeyEvent)
    downstream_of = RelationshipTo(
        Cons.KEY_EVENT_NODE_LABEL.replace(" ", ""), Cons.DOWNSTREAM_OF, model=DownstreamOf
    )


class AdverseOutcome(StructuredNode):
    """Adverse outcome node."""

    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty()
    datasource = StringProperty(required=True)
    node_type = Cons.AO_NODE_LABEL

    # outgoing relations (AO -> KeyEvent)
    downstream_of = RelationshipTo(KeyEvent, Cons.ASSOCIATED_WITH, model=AssociatedWith)


def exporter(
    network,
    uri,
    username,
    password,
    neo4j_import_folder,
    network_name: str = "Network",
):
    """Import the network to neo4j.

    :param network: NetworkX network
    :param uri: URI for Neo4j database
    :param username: user name for Neo4j database
    :param password: password for Neo4j database
    :param neo4j_import_folder: exact path to neo4j database import folder
    :param network_name: network name given by users

    Usage example:
    >> network = nxGraph
    >> uri = "neo4j://localhost:7687"
    >> username = "neo4j"
    >> password = "biodatafuse"
    >> neo4j_import_folder = "../../neo4j-community-5.13.0/import/"
    >> network_name = "Network"
    >> exporter(
        network,
        uri,
        username,
        password,
        neo4j_import_folder, network_name
    )
    """
    # credentials
    uri_info = uri
    auth_info = (username, password)

    # write network to neo4j's import folder temporarily
    filename = network_name + ".graphml"
    save_graph_to_graphml(g=network, output_path=neo4j_import_folder + filename)

    # connect to database
    with GraphDatabase.driver(uri_info, auth=auth_info) as driver:
        driver.verify_connectivity()

    # clean database
    driver.execute_query(
        "MATCH (n) DETACH DELETE n",
        database_="neo4j",
    )

    # import network
    driver.execute_query(
        "CALL apoc.import.graphml($filename, {readLabels: true})",
        filename=filename,
        database_="neo4j",
    )

    # assign node types
    driver.execute_query(
        "MATCH (n) WITH COLLECT(DISTINCT n.labels) AS propertyValues, n "
        + "UNWIND propertyValues AS propValue "
        + "MATCH (n) "
        + "WHERE n.labels = propValue "
        + "WITH n, propValue AS newLabel "
        + "CALL apoc.create.addLabels(n, [newLabel]) YIELD node "
        + "RETURN node",
        database_="neo4j",
    )

    # assign edge types
    driver.execute_query(
        "MATCH (source)-[r]->(target) "
        + "WITH COLLECT(DISTINCT r.interaction) AS propertyValues, r "
        + "UNWIND propertyValues AS propValue "
        + "MATCH (source)-[r]->(target) "
        + "WHERE r.interaction = propValue "
        + "WITH r, source, target, propValue AS newType "
        + "CALL apoc.create.relationship(source, newType, {}, target) YIELD rel "
        + "DELETE r "
        + "RETURN rel ",
        database_="neo4j",
    )

    driver.close()


def connect_db(uri: str, username: str, password: str):
    """Connect to the Neo4j database."""
    # This requires a local community version installation of Neo4j
    my_driver = GraphDatabase().driver(uri, auth=(username, password))
    config.DRIVER = my_driver

    # Delete all nodes and relationships
    db.cypher_query("MATCH ()-[r]-() DELETE r")  # delete all relationships
    db.cypher_query("MATCH (n) DETACH DELETE n")  # delete all nodes


@typing.no_type_check
def load_graph(g: nx.MultiDiGraph, uri: str, username: str, password: str):
    """Load the BDF graph to Neo4j."""
    connect_db(uri, username, password)

    source_nodes = defaultdict(dict)
    for node, node_data in tqdm(g.nodes(data=True), desc="Loading nodes"):
        node_type = node_data[Cons.LABEL].lower()
        if node_type == Cons.GENE_NODE_LABEL.lower():
            n = Gene(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
                Ensembl=node_data.get(Cons.ENSEMBL, None),
            ).save()
        elif node_type == Cons.DISEASE_NODE_LABEL.lower():
            n = Disease(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
                MONDO=node_data.get(Cons.MONDO, None),
                UMLS=node_data.get(Cons.UMLS, None),
            ).save()
        elif node_type == Cons.PATHWAY_NODE_LABEL.lower():
            n = Pathway(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.GO_BP_NODE_LABEL.lower():
            n = BiologicalProcess(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.GO_MF_NODE_LABEL.lower():
            n = MolecularFunction(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.GO_CC_NODE_LABEL.lower():
            n = CellularComponent(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.SIDE_EFFECT_NODE_LABEL.lower():
            n = SideEffect(
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.COMPOUND_NODE_LABEL.lower():
            n = Compound(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
                chembl_id=node_data.get(Cons.CHEMBL_ID, None),
                drugbank_id=node_data.get(Cons.DRUGBANK_ID, None),
                compound_cid=node_data.get("compound_cid", None),
                clinical_trial_phase=node_data.get("clincal_trial_phase", None),
                is_approved=node_data.get("is_approved", None),
                adverse_effect_count=node_data.get("adverse_effect_count", None),
            ).save()
        elif node_type == Cons.ANATOMICAL_NODE_LABEL.lower():
            n = AnatomicalEntity(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.AOP_NODE_LABEL.lower():
            n = AdverseOutcomePathway(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.MIE_NODE_LABEL.lower():
            n = MolecularInitiatingEvent(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        elif node_type == Cons.KEY_EVENT_NODE_LABEL.lower():
            n = KeyEvent(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
                organ=node_data.get("organ", None),
            ).save()
        elif node_type == Cons.AO_NODE_LABEL.lower():
            n = AdverseOutcome(
                idx=node_data[Cons.ID],
                name=node_data[Cons.NAME],
                datasource=node_data[Cons.DATASOURCE],
            ).save()
        else:
            raise ValueError(f"Node type {node_type} not found in Neo4j template")

        source_nodes[node] = n

    edges = []
    # add edges
    for snode, tnode, data in tqdm(g.edges(data=True), desc="Loading edges"):
        source_node = source_nodes[snode]
        target_node = source_nodes[tnode]
        edge_type = data[Cons.LABEL].upper()

        try:
            if edge_type == Cons.ASSOCIATED_WITH:
                if target_node.node_type == Cons.DISEASE_NODE_LABEL:
                    source_node.associated_with_disease.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                elif target_node.node_type == Cons.AOP_NODE_LABEL:
                    source_node.associated_with_aop.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                elif source_node.node_type == Cons.KEY_EVENT_NODE_LABEL:
                    source_node.associated_with.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                else:
                    raise ValueError(
                        f"Edge type {edge_type} not found in Neo4j template - {source_node.node_type} -> {target_node.node_type}"
                    )

            elif edge_type == Cons.PART_OF:
                if target_node.node_type == Cons.PATHWAY_NODE_LABEL:
                    source_node.part_of_pathway.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                elif target_node.node_type == Cons.GO_BP_NODE_LABEL:
                    source_node.part_of_biological_process.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                elif target_node.node_type == Cons.GO_MF_NODE_LABEL:
                    source_node.part_of_molecular_function.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                elif target_node.node_type == Cons.GO_CC_NODE_LABEL:
                    source_node.part_of_cellular_component.connect(
                        target_node, {"datasource": data[Cons.DATASOURCE]}
                    )
                    edges.append([source_node, target_node, edge_type])
                else:
                    raise ValueError(
                        f"Edge type {edge_type} not found in Neo4j template - {source_node.node_type} -> {target_node.node_type}"
                    )

            elif edge_type == Cons.INTERACTS_WITH:
                source_node.interacts_with.connect(
                    target_node,
                    {"datasource": data[Cons.DATASOURCE], "score": data.get(Cons.SCORE, None)},
                )
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.ACTIVATES:
                source_node.activates.connect(target_node, {"datasource": data[Cons.DATASOURCE]})
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.HAS_SIDE_EFFECT:
                source_node.has_side_effect.connect(
                    target_node, {"datasource": data[Cons.DATASOURCE]}
                )
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.TREATS:
                source_node.treats.connect(target_node, {"datasource": data[Cons.DATASOURCE]})
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.INHIBITS:
                source_node.inhibits.connect(target_node, {"datasource": data[Cons.DATASOURCE]})
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.EXPRESSED_BY:
                source_node.expressed_by.connect(
                    target_node,
                    {
                        "datasource": data[Cons.DATASOURCE],
                        "expression_level": data.get(Cons.EXPRESSION_LEVEL, None),
                        "developmental_stage": data.get(Cons.DEVELOPMENTAL_STAGE_NAME, None),
                        "developmental_stage_id": data.get(Cons.DEVELOPMENTAL_STAGE_ID, None),
                        "confidence_level": data.get(Cons.CONFIDENCE_LEVEL_NAME, None),
                        "confidence_level_id": data.get(Cons.CONFIDENCE_ID, None),
                    },
                )
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.UPSTREAM_OF:
                source_node.upstream_of.connect(target_node, {"datasource": data[Cons.DATASOURCE]})
                edges.append([source_node, target_node, edge_type])

            elif edge_type == Cons.DOWNSTREAM_OF:
                source_node.downstream_of.connect(
                    target_node, {"datasource": data[Cons.DATASOURCE]}
                )
                edges.append([source_node, target_node, edge_type])
            else:
                raise ValueError(f"Edge type {edge_type} not found in Neo4j template")

        except AttributeError:
            print(snode, tnode, data)
            print(
                f"AttributeError: {edge_type} not found in Neo4j template - {source_node} -> {target_node}"
            )

    # Check to ensure that all edges and nodes are loaded
    assert len(source_nodes) == len(g.nodes)
    assert len(edges) == len(g.edges)
