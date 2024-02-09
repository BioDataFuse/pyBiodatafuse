from pyBiodatafuse.annotators.stringdb import get_version_stringdb, get_ppi
import pytest
import pandas as pd
import numpy as np
from numpy import nan

def test_get_version_stringdb():

    obtained_version = get_version_stringdb();

    expected_version = [{'string_version': '12.0', 'stable_address': 'https://version-12-0.string-db.org'}]

    assert obtained_version ==  expected_version


def test_get_ppi(bridgedb_dataframe):

    obtained_data, metadata = get_ppi(bridgedb_dataframe)

    expected_data = pd.Series([[{'stringdb_link_to': 'CHRNA1', 'score': 0.543}, {'stringdb_link_to': 'ALG2', 'score': 0.633}],
                               [{'stringdb_link_to': 'ALG14', 'score': 0.633}],
                               [{'stringdb_link_to': 'ALG14', 'score': 0.543}]])
    expected_data.name = 'stringdb'

    pd.testing.assert_series_equal(obtained_data["stringdb"], expected_data)


@pytest.fixture(scope="module")
def bridgedb_dataframe():
    return pd.DataFrame({
        'identifier': ['ALG14', 'ALG2', 'CHRNA1'],
        'identifier.source': ['HGNC', 'HGNC', 'HGNC'],
        'target': ['ENSG00000172339', 'ENSG00000119523', 'ENSG00000138435'],
        'target.source': ['Ensembl', 'Ensembl', 'Ensembl']
    })

