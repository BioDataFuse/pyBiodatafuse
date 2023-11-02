# -*- coding: utf-8 -*-

"""Trivial version test."""

import unittest

from pyBiodatafuse.annotators.disgenet import get_version_disgenet as disgenet_version
from pyBiodatafuse.annotators.opentargets import get_version as opentargets_version
from pyBiodatafuse.annotators.wikidata import get_version_wikidata as wikidata_version
from pyBiodatafuse.id_mapper import get_version_datasource_bridgedb as datasource_bridgedb_version
from pyBiodatafuse.id_mapper import get_version_webservice_bridgedb as webservice_bridgedb_version
from pyBiodatafuse.version import get_version


class TestVersion(unittest.TestCase):
    """Trivially test a version."""

    def test_version_type(self):
        """Test the version is a string.

        This is only meant to be an example test.
        """
        version = get_version()
        self.assertIsInstance(version, str)

    def test_opentargets_version_type(self):
        """Test the version from OpenTargets is a dictionary."""
        version = opentargets_version()
        self.assertIsInstance(version, dict)

    def test_wikidata_version_type(self):
        """Test the version from Wikidata is a dictionary."""
        version = wikidata_version()
        self.assertIsInstance(version, dict)

    def test_disgenet_version_type(self):
        """Test the version from DisGeNET is a dictionary."""
        version = disgenet_version()
        self.assertIsInstance(version, dict)

    def test_wb_bridgedb_version_type(self):
        """Test the version from BridgeDb is a dictionary."""
        version = webservice_bridgedb_version()
        self.assertIsInstance(version, dict)

    def test_db_bridgedb_version_type(self):
        """Test the version from datasources in BridgeDb is a list."""
        version = datasource_bridgedb_version()
        self.assertIsInstance(version, list)
