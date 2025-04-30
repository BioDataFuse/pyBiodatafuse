# coding: utf-8

"""Python module to export a NetworkX graph to neo4j-compliant format and set the styling for Neo4j."""

from collections import defaultdict
from tqdm import tqdm
import networkx as nx
from neo4j import GraphDatabase

from pyBiodatafuse.graph.saver import save_graph_to_graphml

from neomodel import (
    StructuredNode,
    StringProperty,
    RelationshipTo,
    StructuredRel,
    FloatProperty,
    BooleanProperty,
    RelationshipFrom,
    config,
    db,
)


# Predifinding relationship classes
class AssociatedWith(StructuredRel):
    datasource = StringProperty()
    score = FloatProperty()
    ei = FloatProperty()


class PartOf(StructuredRel):
    datasource = StringProperty()


class InteractsWith(StructuredRel):
    datasource = StringProperty()
    score = FloatProperty()


class Activates(StructuredRel):
    datasource = StringProperty()


class HasSideEffect(StructuredRel):
    datasource = StringProperty()


class Treats(StructuredRel):
    datasource = StringProperty()


class Inhibits(StructuredRel):
    datasource = StringProperty()


# Defining the nodes
class Gene(StructuredNode):
    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    Ensembl = StringProperty()

    # outgoing relations
    interacts_with = RelationshipTo("Gene", "INTERACTS_WITH", model=InteractsWith)  # PPI
    part_of_pathway = RelationshipTo("Pathway", "PART_OF", model=PartOf)  # Gene -> Pathway
    part_of_biological_process = RelationshipTo(
        "BiologicalProcess", "PART_OF", model=PartOf
    )  # Gene -> BiologicalProcess
    part_of_molecular_function = RelationshipTo(
        "MolecularFunction", "PART_OF", model=PartOf
    )  # Gene -> MolecularFunction
    part_of_cellular_component = RelationshipTo(
        "CellularComponent", "PART_OF", model=PartOf
    )  # Gene -> CellularComponent
    associated_with = RelationshipTo(
        "Disease", "ASSOCIATED_WITH", model=AssociatedWith
    )  # Gene -> Disease


class Disease(StructuredNode):
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

    # incoming relations (Gene -> Disease)
    associated_with = RelationshipFrom(Gene, "ASSOCIATED_WITH", model=AssociatedWith)


class Pathway(StructuredNode):
    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)

    # incoming relations (Gene -> Pathway)
    part_of = RelationshipFrom(Gene, "PART_OF", model=PartOf)


class BiologicalProcess(StructuredNode):
    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)

    # incoming relations (Gene -> BiologicalProcess)
    part_of = RelationshipFrom(Gene, "PART_OF", model=PartOf)


class MolecularFunction(StructuredNode):
    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)

    # incoming relations (Gene -> MolecularFunction)
    part_of = RelationshipFrom(Gene, "PART_OF", model=PartOf)


class CellularComponent(StructuredNode):
    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)

    # incoming relations (Gene -> CellularComponent)
    part_of = RelationshipFrom(Gene, "PART_OF", model=PartOf)


class SideEffect(StructuredNode):
    name = StringProperty(required=True, unique_index=True, unique=True)
    datasource = StringProperty(required=True)


class Compound(StructuredNode):
    idx = StringProperty(required=True, unique_index=True, unique=True)
    name = StringProperty(required=True)
    datasource = StringProperty(required=True)
    chembl_id = StringProperty()
    drugbank_id = StringProperty()
    compound_cid = StringProperty()
    clinical_trial_phase = StringProperty()
    is_approved = BooleanProperty()
    adverse_effect_count = FloatProperty()

    # outgoing relations (Compound -> Gene, SideEffect, Disease)
    activates = RelationshipTo(Gene, "ACTIVATES", model=Activates)
    has_side_effect = RelationshipTo(SideEffect, "HAS_SIDE_EFFECT", model=HasSideEffect)
    treats = RelationshipTo(Disease, "TREATS", model=Treats)
    inhibits = RelationshipTo(Gene, "INHIBITS", model=Inhibits)

    # incoming relations (Gene -> Compound)
    associated_with = RelationshipFrom(Gene, "ASSOCIATED_WITH", model=AssociatedWith)


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


def load_graph(g: nx.MultiDiGraph, uri: str, username: str, password: str):
    """Loading the BDF graph to Neo4j."""

    connect_db(uri, username, password)

    source_nodes = defaultdict()
    for node, node_data in tqdm(g.nodes(data=True), desc="Loading nodes"):
        node_type = node_data["labels"].lower()
        if node_type == "gene":
            n = Gene(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
                Ensembl=node_data.get("Ensembl", None),
            ).save()
        elif node_type == "disease":
            n = Disease(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
                MONDO=node_data.get("MONDO", None),
                UMLS=node_data.get("UMLS", None),
            ).save()
        elif node_type == "pathway":
            n = Pathway(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
            ).save()
        elif node_type == "biological process":
            n = BiologicalProcess(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
            ).save()
        elif node_type == "molecular function":
            n = MolecularFunction(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
            ).save()
        elif node_type == "cellular component":
            n = CellularComponent(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
            ).save()
        elif node_type == "side effect":
            n = SideEffect(
                name=node_data["name"],
                datasource=node_data["datasource"],
            ).save()
        elif node_type == "compound":
            n = Compound(
                idx=node_data["id"],
                name=node_data["name"],
                datasource=node_data["datasource"],
                chembl_id=node_data.get("chembl_id", None),
                drugbank_id=node_data.get("drugbank_id", None),
                compound_cid=node_data.get("compound_cid", None),
                clinical_trial_phase=node_data.get("clincal_trial_phase", None),
                is_approved=node_data.get("is_approved", None),
                adverse_effect_count=node_data.get("adverse_effect_count", None),
            ).save()
        else:
            raise ValueError(f"Node type {node_type} not found in Neo4j template")

        source_nodes[node] = n

    # add edges
    for snode, tnode, data in tqdm(g.edges(data=True), desc="Loading edges"):
        source_node = source_nodes[snode]
        target_node = source_nodes[tnode]
        edge_type = data["label"].lower()

        try:
            if edge_type == "associated_with":
                source_node.associated_with.connect(target_node, {"datasource": data["datasource"]})

            elif edge_type == "part_of":
                if source_node.labels == "Pathway":
                    source_node.part_of_pathway.connect(
                        target_node, {"datasource": data["datasource"]}
                    )
                elif source_node.labels == "BiologicalProcess":
                    source_node.part_of_biological_process.connect(
                        target_node, {"datasource": data["datasource"]}
                    )
                elif source_node.labels == "MolecularFunction":
                    source_node.part_of_molecular_function.connect(
                        target_node, {"datasource": data["datasource"]}
                    )
                elif source_node.labels == "CellularComponent":
                    source_node.part_of_cellular_component.connect(
                        target_node, {"datasource": data["datasource"]}
                    )

            elif edge_type == "interacts_with":
                source_node.interacts_with.connect(
                    target_node,
                    {"datasource": data["datasource"], "score": data.get("score", None)},
                )
            elif edge_type == "activates":
                source_node.activates.connect(target_node, {"datasource": data["datasource"]})

            elif edge_type == "has_side_effect":
                source_node.has_side_effect.connect(target_node, {"datasource": data["datasource"]})

            elif edge_type == "treats":
                source_node.treats.connect(target_node, {"datasource": data["datasource"]})

            elif edge_type == "inhibits":
                source_node.inhibits.connect(target_node, {"datasource": data["datasource"]})

            else:
                raise ValueError(f"Edge type {edge_type} not found in Neo4j template")

        except AttributeError:
            print(
                f"AttributeError: {edge_type} not found in Neo4j template - {source_node} -> {target_node}"
            )
