# coding: utf-8

"""Test cases for dataset_provenance.py module."""

import unittest
from datetime import datetime

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.dataset_provenance import (
    DatasetProvenanceTracker,
    add_dataset_provenance_to_graph,
    create_dataset_node,
    get_datasource_for_column,
    record_datasource_from_metadata,
)


class TestDatasetProvenanceTracker(unittest.TestCase):
    """Test the DatasetProvenanceTracker class."""

    def test_record_dataset_usage(self):
        """Test recording dataset usage."""
        tracker = DatasetProvenanceTracker(base_uri="http://example.org/")
        
        tracker.record_dataset_usage(datasource="StringDB", version="12.0")
        tracker.record_dataset_usage(datasource="Bgee", version="15.1")
        
        datasets = tracker.get_used_datasets()
        
        self.assertIn("StringDB", datasets)
        self.assertIn("Bgee", datasets)
        self.assertEqual(datasets["StringDB"]["version"], "12.0")
        self.assertEqual(datasets["Bgee"]["version"], "15.1")
        self.assertEqual(datasets["StringDB"]["query_count"], 1)

    def test_record_dataset_usage_multiple_times(self):
        """Test recording the same dataset multiple times increments count."""
        tracker = DatasetProvenanceTracker(base_uri="http://example.org/")
        
        tracker.record_dataset_usage(datasource="StringDB")
        tracker.record_dataset_usage(datasource="StringDB")
        tracker.record_dataset_usage(datasource="StringDB")
        
        datasets = tracker.get_used_datasets()
        
        self.assertEqual(datasets["StringDB"]["query_count"], 3)


class TestCreateDatasetNode(unittest.TestCase):
    """Test the create_dataset_node function."""

    def test_create_dataset_node_basic(self):
        """Test creating a basic dataset node."""
        g = Graph()
        base_uri = "http://example.org/"
        
        dataset_uri = create_dataset_node(
            g=g,
            datasource="StringDB",
            base_uri=base_uri,
            version="12.0",
        )
        
        # Check the URI is correct
        self.assertEqual(str(dataset_uri), "http://example.org/dataset/stringdb")
        
        # Check it has the DCAT Dataset type
        self.assertIn(
            (dataset_uri, RDF.type, URIRef(Cons.DCAT_TYPES["dataset"])),
            g
        )
        
        # Check it has the VoID Dataset type
        self.assertIn(
            (dataset_uri, RDF.type, URIRef(Cons.VOID_TYPES["dataset"])),
            g
        )

    def test_create_dataset_node_with_version(self):
        """Test creating a dataset node with version info."""
        g = Graph()
        base_uri = "http://example.org/"
        
        dataset_uri = create_dataset_node(
            g=g,
            datasource="Bgee",
            base_uri=base_uri,
            version="15.1",
        )
        
        # Check version is recorded with PAV predicate
        version_triple = (
            dataset_uri,
            URIRef(Cons.PAV_PREDICATES["version"]),
            Literal("15.1"),
        )
        self.assertTrue(any(
            str(o) == "15.1" 
            for s, p, o in g.triples((dataset_uri, URIRef(Cons.PAV_PREDICATES["version"]), None))
        ))

    def test_create_dataset_node_with_endpoint(self):
        """Test creating a dataset node with endpoint URL."""
        g = Graph()
        base_uri = "http://example.org/"
        
        dataset_uri = create_dataset_node(
            g=g,
            datasource="WikiPathways",
            base_uri=base_uri,
            endpoint_url="https://sparql.wikipathways.org/sparql",
        )
        
        # Check the endpoint is recorded
        access_url_triples = list(g.triples((
            dataset_uri,
            URIRef(Cons.DCAT_PREDICATES["access_url"]),
            None
        )))
        self.assertTrue(len(access_url_triples) > 0)

    def test_create_dataset_node_with_timestamp(self):
        """Test creating a dataset node with query timestamp."""
        g = Graph()
        base_uri = "http://example.org/"
        query_time = datetime(2025, 1, 15, 10, 30, 0)
        
        dataset_uri = create_dataset_node(
            g=g,
            datasource="CompoundWiki",
            base_uri=base_uri,
            query_time=query_time,
        )
        
        # Check the timestamp is recorded
        timestamp_triples = list(g.triples((
            dataset_uri,
            URIRef(Cons.PAV_PREDICATES["retrieved_on"]),
            None
        )))
        self.assertTrue(len(timestamp_triples) > 0)


class TestAddDatasetProvenanceToGraph(unittest.TestCase):
    """Test the add_dataset_provenance_to_graph function."""

    def test_add_dataset_provenance(self):
        """Test adding dataset provenance to graph."""
        g = Graph()
        base_uri = "http://example.org/"
        graph_uri = "http://example.org/graph/v1.0"
        
        tracker = DatasetProvenanceTracker(base_uri)
        tracker.record_dataset_usage("StringDB", version="12.0")
        tracker.record_dataset_usage("Bgee", version="15.1")
        
        add_dataset_provenance_to_graph(g, base_uri, graph_uri, tracker)
        
        # Check that the graph has a DCAT Catalog type
        graph_resource = URIRef(graph_uri)
        self.assertIn(
            (graph_resource, RDF.type, URIRef(Cons.DCAT_TYPES["catalog"])),
            g
        )
        
        # Check that datasets are linked via PROV wasDerivedFrom
        derived_from_triples = list(g.triples((
            graph_resource,
            URIRef(Cons.PROV_PREDICATES["was_derived_from"]),
            None
        )))
        self.assertEqual(len(derived_from_triples), 2)


class TestGetDatasourceForColumn(unittest.TestCase):
    """Test the get_datasource_for_column function."""

    def test_known_columns(self):
        """Test mapping known columns to data sources."""
        self.assertEqual(
            get_datasource_for_column(Cons.BGEE_GENE_EXPRESSION_LEVELS_COL),
            Cons.BGEE
        )
        self.assertEqual(
            get_datasource_for_column(Cons.STRING_INTERACT_COL),
            Cons.STRING
        )
        self.assertEqual(
            get_datasource_for_column(Cons.COMPOUNDWIKI_COL),
            Cons.COMPOUNDWIKI
        )

    def test_unknown_column(self):
        """Test that unknown columns return None."""
        result = get_datasource_for_column("unknown_column")
        self.assertIsNone(result)


class TestRecordDatasourceFromMetadata(unittest.TestCase):
    """Test the record_datasource_from_metadata function."""

    def test_record_from_metadata_list(self):
        """Test recording datasources from BDF metadata format."""
        tracker = DatasetProvenanceTracker(base_uri="http://example.org/")
        
        metadata_list = [
            {
                "datasource": "StringDB",
                "metadata": {"source_version": "12.0"},
                "query": {
                    "url": "https://string-db.org/api",
                    "date": "2025-01-15 10:30:00"
                }
            },
            {
                "datasource": "Bgee",
                "metadata": {"source_version": "15.1"},
                "query": {
                    "url": "https://www.bgee.org/sparql/",
                    "date": "2025-01-15 10:35:00"
                }
            }
        ]
        
        record_datasource_from_metadata(tracker, metadata_list)
        
        datasets = tracker.get_used_datasets()
        
        self.assertIn("StringDB", datasets)
        self.assertIn("Bgee", datasets)
        self.assertEqual(datasets["StringDB"]["version"], "12.0")
        self.assertEqual(datasets["Bgee"]["version"], "15.1")


class TestRDFGraphProvenance(unittest.TestCase):
    """Integration test for provenance in the RDF graph."""

    def test_provenance_triples_count(self):
        """Test that provenance adds expected number of triples."""
        g = Graph()
        base_uri = "http://example.org/"
        graph_uri = "http://example.org/graph/v1.0"
        
        tracker = DatasetProvenanceTracker(base_uri)
        tracker.record_dataset_usage("StringDB", version="12.0")
        
        initial_len = len(g)
        add_dataset_provenance_to_graph(g, base_uri, graph_uri, tracker)
        final_len = len(g)
        
        # Should add multiple triples for the catalog + dataset
        self.assertGreater(final_len, initial_len)
        # At minimum: 1 catalog type, 1 catalog title, 1 derived_from, 1 source, 
        # plus dataset triples (type, void type, label, title, version, etc.)
        self.assertGreater(final_len - initial_len, 5)


if __name__ == "__main__":
    unittest.main()
