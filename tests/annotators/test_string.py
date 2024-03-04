#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Stringdb annotator."""

from unittest.mock import Mock, patch

import pandas as pd
import pytest

from pyBiodatafuse.annotators.stringdb import get_ppi, get_version_stringdb


@patch("pyBiodatafuse.annotators.stringdb.requests.get")
def test_get_version_stringdb(mock_requests_get):
    """Test the get_version_stringdb."""
    mock_requests_get.return_value.json.return_value = [
        {"string_version": "12.0", "stable_address": "https://version-12-0.string-db.org"}
    ]

    obtained_version = get_version_stringdb()

    expected_version = [
        {"string_version": "12.0", "stable_address": "https://version-12-0.string-db.org"}
    ]

    assert obtained_version == expected_version


@patch("pyBiodatafuse.annotators.stringdb.requests.post")
@patch("pyBiodatafuse.annotators.stringdb.requests.get")
def test_get_ppi(mock_requests_get, mock_requests_post, bridgedb_dataframe):
    """Test the get_ppi function."""
    mock_requests_post.side_effect = [
        Mock(
            content=(
                b"queryIndex\tstringId\tncbiTaxonId\ttaxonName\tpreferredName\tannotation\n"
                b"0\t9606.ENSP00000261007\t9606\tHomo sapiens\tCHRNA1\tAcetylcholine receptor"
                b"subunit alpha; After binding acetylcholine, the AChR responds by an extensive "
                b"change in conformation that affects all subunits and leads to opening of an "
                b"ion-conducting channel across the plasma membrane.\n"
                b"1\t9606.ENSP00000359224\t9606\tHomo sapiens\tALG14\tUDP-N-acetylglucosamine "
                b"transferase subunit ALG14 homolog; May be involved in protein N-glycosylation. "
                b"May play a role in the second step of the dolichol-linked oligosaccharide pathway. "
                b"May anchor the catalytic subunit ALG13 to the ER. Belongs to the ALG14 family.\n"
                b"2\t9606.ENSP00000417764\t9606\tHomo sapiens\tALG2\tAlpha-1,3/1,6-mannosyltransferase "
                b"ALG2; Mannosylates Man(2)GlcNAc(2)-dolichol diphosphate and Man(1)GlcNAc(2)-dolichol "
                b"diphosphate to form Man(3)GlcNAc(2)-dolichol diphosphate; Belongs to the "
                b"glycosyltransferase group 1 family. Glycosyltransferase 4 subfamily.\n"
            )
        ),
        Mock(
            content=(
                b"stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\t"
                b"ncbiTaxonId\tscore\tnscore\tfscore\tpscore\t"
                b"ascore\tescore\tdscore\ttscore\n9606.ENSP00000261007\t9606.ENSP00000359224\t"
                b"CHRNA1\tALG14\t9606\t0.543\t0\t0\t0\t0\t0\t0\t0.543\n9606.ENSP00000359224\t"
                b"9606.ENSP00000417764\tALG14\tALG2\t9606\t0.633\t0\t0\t0\t0.067\t0\t0.119\t0.589\n"
            )
        ),
    ]

    mock_requests_get.return_value.json.return_value = [
        {"string_version": "12.0", "stable_address": "https://version-12-0.string-db.org"}
    ]

    obtained_data, metadata = get_ppi(bridgedb_dataframe)

    expected_data = pd.Series(
        [
            [
                {"stringdb_link_to": "CHRNA1", "score": 0.543},
                {"stringdb_link_to": "ALG2", "score": 0.633},
            ],
            [{"stringdb_link_to": "ALG14", "score": 0.633}],
            [{"stringdb_link_to": "ALG14", "score": 0.543}],
        ]
    )
    expected_data.name = "stringdb"

    pd.testing.assert_series_equal(obtained_data["stringdb"], expected_data)


@pytest.fixture(scope="module")
def bridgedb_dataframe():
    """Reusable sample Pandas DataFrame to be used as input for the tests."""
    return pd.DataFrame(
        {
            "identifier": ["ALG14", "ALG2", "CHRNA1"],
            "identifier.source": ["HGNC", "HGNC", "HGNC"],
            "target": ["ENSG00000172339", "ENSG00000119523", "ENSG00000138435"],
            "target.source": ["Ensembl", "Ensembl", "Ensembl"],
        }
    )
