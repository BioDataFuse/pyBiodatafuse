"""Graph summary functions."""

import warnings
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import plotly.express as px
import seaborn as sns
from tabulate import tabulate

from pyBiodatafuse.analyzer.explorer.patent import get_patent_from_pubchem
from pyBiodatafuse.constants import COMPOUND_NAMESPACE_MAPPER
from pyBiodatafuse.graph.generator import build_networkx_graph, load_dataframe_from_pickle


class BioGraph(nx.MultiDiGraph):
    """BioGraph class to analyze the graph."""

    def __init__(self, graph=None, graph_path=None, graph_format="pickle", disease_df=None):
        """Initialize the BioGraph class.

        :param graph: networkx graph object
        :param graph_path: path to the graph file
        :param graph_format: format of the graph file
        :param disease_df: disease dataframe to build the graph
        :raises ValueError: if graph_format is not 'pickle' or 'gml'
        """
        if graph:
            self.graph = graph
        elif graph_path:
            if graph_format == "pickle" and disease_df:
                self.graph = build_networkx_graph(
                    load_dataframe_from_pickle(graph_path),
                    disease_df,
                )
            elif graph_format == "pickle" and not disease_df:
                warnings.warn(
                    "Disease dataframe not provided. Loading graph without disease data.",
                    stacklevel=2,
                )
                self.graph = nx.read_gpickle(graph_path)
            elif graph_format == "gml":
                self.graph = nx.read_gml(graph_path)
            elif graph_format == "graphml":
                self.graph = nx.read_graphml(graph_path)
            else:
                raise ValueError("graph_format must be either 'pickle', 'gml' or 'graphml'")

        self.node_count = self.count_nodes_by_type()
        self.edge_count = self.count_edge_by_type()
        self.node_source_count = self.count_nodes_by_data_source()
        self.edge_source_count = self.count_edge_by_data_source()
        self.graph_summary = self.get_graph_summary()

    def get_graph_summary(self) -> str:
        """Display graph summary."""
        stats = [
            ("Nodes", self.graph.number_of_nodes()),
            ("Edges", self.graph.number_of_edges()),
            ("Components", nx.number_weakly_connected_components(self.graph)),
            ("Network Density", "{:.2E}".format(nx.density(self.graph))),
        ]
        return tabulate(stats, tablefmt="html")

    def _plot_type_count(
        self, count_df: pd.DataFrame, interactive: bool = False, count_type: str = "Node"
    ) -> None:
        """Plot the type counts on barplot."""
        if count_type == "Node":
            plot_title = "Node Type Count"
            x_label = "Node Type"
            x_col = "node_type"
        elif count_type == "Edge":
            plot_title = "Edge Type Count"
            x_label = "Edge Type"
            x_col = "edge_type"
        else:
            raise ValueError("count_type must be either 'Node' or 'Edge'")

        if interactive:
            fig = px.bar(count_df, x=x_col, y="count")
            fig.update_layout(title=plot_title)
            fig.update_xaxes(title_text=x_label)
            fig.update_yaxes(title_text="Count")
            fig.show()
        else:
            plt.figure(figsize=(10, 6))
            sns.barplot(x=x_col, y="count", data=count_df)
            # counts on top of bar
            for i in range(count_df.shape[0]):
                count = count_df.iloc[i]["count"]
                plt.text(i, count, count, ha="center")
            plt.title(plot_title)
            plt.xlabel(x_label)
            plt.ylabel("Count")
            plt.tight_layout()
            plt.show()

    def count_nodes_by_type(
        self, plot: bool = False, interactive: bool = False
    ) -> Optional[pd.DataFrame]:
        """Count the differnent nodes type in the graph."""
        node_data = pd.DataFrame(self.graph.nodes(data=True), columns=["node", "data"])
        node_data["node_type"] = node_data["data"].apply(lambda x: x["labels"])
        node_count = node_data["node_type"].value_counts().reset_index()
        node_count = node_count.sort_values(by="count", ascending=False)

        if plot:
            self._plot_type_count(node_count, interactive, count_type="Node")
            return None

        return node_count

    def count_edge_by_type(
        self, plot: bool = False, interactive: bool = False
    ) -> Optional[pd.DataFrame]:
        """Count the different edge types in the graph."""
        edge_data = pd.DataFrame(self.graph.edges(data=True), columns=["source", "target", "data"])
        edge_data["edge_type"] = edge_data["data"].apply(lambda x: x["label"])
        edge_count = edge_data["edge_type"].value_counts().reset_index()
        edge_count = edge_count.sort_values(by="count", ascending=False)

        if plot:
            self._plot_type_count(edge_count, interactive, count_type="Edge")
            return None

        return edge_count

    def _plot_source_count(self, source_count_df: pd.DataFrame, count_type: str = "Node") -> None:
        """Plot count of nodes or edges by source."""
        if count_type == "Node":
            x_col = "node_type"
            x_color = "node_source"
        elif count_type == "Edge":
            x_col = "edge_type"
            x_color = "edge_source"

        fig = px.bar(source_count_df, x=x_col, y="count", color=x_color)
        fig.update_layout(
            title=f"{count_type} count by source",
            xaxis_title=f"{count_type} Type",
            yaxis_title="Count",
        )
        fig.show()

    def count_nodes_by_data_source(self, plot: bool = False) -> Optional[pd.DataFrame]:
        """Get the count of nodes by data source."""
        node_data = pd.DataFrame(self.graph.nodes(data=True), columns=["node", "data"])
        node_data["node_type"] = node_data["data"].apply(lambda x: x["labels"])
        node_data["node_source"] = node_data["data"].apply(lambda x: x["datasource"])
        node_source_count = (
            node_data.groupby(["node_type", "node_source"]).size().reset_index(name="count")
        )

        if plot:
            self._plot_source_count(node_source_count, count_type="Node")
            return None

        return node_source_count

    def count_edge_by_data_source(self, plot: bool = False) -> Optional[pd.DataFrame]:
        """Get the count of edges by data source."""
        edge_data = pd.DataFrame(self.graph.edges(data=True), columns=["source", "target", "data"])
        edge_data["edge_type"] = edge_data["data"].apply(lambda x: x["label"])
        edge_data["edge_source"] = edge_data["data"].apply(lambda x: x["datasource"])
        edge_source_count = (
            edge_data.groupby(["edge_type", "edge_source"]).size().reset_index(name="count")
        )
        edge_source_count = edge_source_count.sort_values(by="count", ascending=False)

        if plot:
            self._plot_source_count(edge_source_count, count_type="Edge")
            return None

        return edge_source_count

    def get_all_nodes_by_labels(self) -> Dict[str, Any]:
        """Get all nodes with their label type."""
        label_dict = {}  # type: Dict[str, Any]
        for node, data in self.graph.nodes(data=True):
            node_type = data["labels"]
            if node_type not in label_dict:
                label_dict[node_type] = []
            label_dict[node_type].append((node, data))

        return label_dict

    def get_all_nodes_by_type(self, label: str) -> list:
        """Get all nodes by specific label type."""
        label_dict = self.get_all_nodes_by_labels()
        return label_dict[label]

    def get_nodes_by_label(self, label: str) -> Optional[list]:
        """Get all nodes by specific label type."""
        if label in self.graph.nodes:
            return self.graph.nodes[label]

        return None

    def get_publications_for_genes(self):
        """Get publications for genes."""
        pass

    def get_patents_for_compounds(self):
        """Get patents for compounds."""
        # TODO: Function in test mode!!
        compound_nodes = self.get_all_nodes_by_type("Compound")

        t = []
        for node_idx, _ in compound_nodes:
            t.append(
                {
                    "target.source": COMPOUND_NAMESPACE_MAPPER[node_idx.split(":")[0]],
                    "target": node_idx.split(":")[1],
                }
            )

        df = pd.DataFrame(t)
        patent_dict = get_patent_from_pubchem(df)
        return patent_dict

    def node_in_graph(self, node_type: str):
        """Check if the node is in the graph."""
        possible_node_type = self.node_count["node_type"].to_list()

        assert node_type in possible_node_type, f"Node type {node_type} not in {possible_node_type}"

    def get_source_interactions(self):
        """Get interactions of a source."""
        pass

    def get_chemical_metatdata(self):
        """Get metadata of a chemical."""
        # """Adverse effects, Clinical trials,
        pass

    def get_subgraph(self, node_types: list):
        """Get subgraph of the graph.

        :param node_types: list of node types
        :raises AssertionError: if node type not in the graph
        :returns: subgraph with the given node types
        """
        assert all(
            [node_type in self.node_count["node_type"].to_list() for node_type in node_types]
        ), "Node type not in the graph. Please use one of the following node types: {}".format(
            self.node_count["node_type"].to_list()
        )

        nodes = [node for node, label in self.graph.nodes(data="labels") if label in node_types]
        subgraph = self.graph.subgraph(nodes)
        return subgraph
