# -*- coding: utf-8 -*-

"""Trivial version test."""

import unittest

from pyBiodatafuse.annotators.opentargets import get_version as opentargets_version
from pyBiodatafuse.annotators.wikidata import get_version_wikidata as wikidata_version
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
