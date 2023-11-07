# coding: utf-8

"""Python module to export a NetworkX graph to neo4j-compliant format and set the styling for Neo4j."""


import networkx as nx
from neo4j import GraphDatabase


def save_graph_to_neo4j_graphml(g: nx.MultiDiGraph, output_path: str):
    """Convert a NetworkX graph to Neo4J graphml file.

    :param g: the NetworkX graph object.
    :param output_path: the output path of the graphml file
    """
    nx.write_graphml(g, output_path, named_key_ids=True)


def export(
    network,
    uri,
    username,
    password,
    path_to_neo4j_import_folder,
    network_name: str = "Network",
):
    """Import the network to neo4j.

    :param network: NetworkX network
    :param uri: URI for Neo4j database
    :param username: user name for Neo4j database
    :param password: password for Neo4j database
    :param path_to_neo4j_import_folder: exact path to neo4j database import folder
    :param network_name: network name given by users

    Usage example:
    >> network = nxGraph
    >> uri = "neo4j://localhost:7687"
    >> username = "neo4j"
    >> paasword = "biodatafuse"
    >> path_to_neo4j_import_folder = "../../neo4j-community-5.13.0/import/"
    >> network_name = "Network"
    >> importNetworkToCytoscape(network, uri, username, password, path_to_neo4j_import_folder, network_name)
    """
    # credentials
    uri_info = uri
    auth_info = (username, password)

    # write network to neo4j's import folder temporarily
    filename = network_name + ".graphml"
    nx.write_graphml(network, path_to_neo4j_import_folder + filename)

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
        "MATCH (n) WITH COLLECT(DISTINCT n.node_type) AS propertyValues, n "
        + "UNWIND propertyValues AS propValue "
        + "MATCH (n) "
        + "WHERE n.node_type = propValue "
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

    # delete temporary file
    # os.remove(path_to_neo4j_import_folder + filename)
    driver.close()
