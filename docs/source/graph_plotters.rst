Graph Plotters
==============

Here you can find the different graph loading, storing, and plotting functions in the package.

.. module:: pyBiodatafuse.graph

Graph Generator
~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.graph.generator.build_networkx_graph

.. autofunction:: pyBiodatafuse.graph.generator.load_dataframe_from_pickle

Graph Savers
~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.graph.saver.save_graph

.. autofunction:: pyBiodatafuse.graph.saver.save_graph_to_graphml

.. autofunction:: pyBiodatafuse.graph.saver.save_graph_to_edgelist

.. autofunction:: pyBiodatafuse.graph.saver.save_graph_to_tsv

.. autofunction:: pyBiodatafuse.graph.saver.save_cytoscape_json

Cytoscape
~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.graph.cytoscape.convert_graph_to_json

.. autofunction:: pyBiodatafuse.graph.cytoscape.load_graph

Neo4J
~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.graph.neo4j.exporter

.. autofunction:: pyBiodatafuse.graph.neo4j.connect_db

.. autofunction:: pyBiodatafuse.graph.neo4j.load_graph

RDF and GraphDB
~~~~~~~~~~~~~~~

.. autoclass:: pyBiodatafuse.graph.rdf.rdf.BDFGraph
   :members:

.. autoclass:: pyBiodatafuse.graph.rdf.graphdb.GraphDBManager
   :members:

.. autofunction:: pyBiodatafuse.graph.rdf.metadata.add_creator_node

.. autofunction:: pyBiodatafuse.graph.rdf.metadata.add_metadata

RDF Utilities
~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.graph.rdf.utils.replace_na_none

.. autofunction:: pyBiodatafuse.graph.rdf.utils.extract_curie

.. autofunction:: pyBiodatafuse.graph.rdf.utils.construct_uri

.. autofunction:: pyBiodatafuse.graph.rdf.utils.add_data_source_node

.. autofunction:: pyBiodatafuse.graph.rdf.utils.get_shapes

.. autofunction:: pyBiodatafuse.graph.rdf.utils.get_shacl_prefixes

.. autofunction:: pyBiodatafuse.graph.rdf.utils.get_node_label

.. autofunction:: pyBiodatafuse.graph.rdf.utils.discover_prefixes_from_graph
