from pyBiodatafuse.annotators.wikipathways import get_version_wikipathways, get_gene_wikipathway
import pytest
import pandas as pd
import numpy as np
from numpy import nan

def test_wikipathways_version():

    obtained_version = get_version_wikipathways();

    expected_version = {"wikipathways_version": 'WikiPathways RDF 20231210'}

    assert obtained_version ==  expected_version


def test_get_gene_wikipathway(bridgedb_dataframe):

    obtained_data, metadata = get_gene_wikipathway(bridgedb_dataframe)

    expected_data = pd.Series([[{'pathwayId': 'WP5153', 'pathwayLabel': 'N-glycan biosynthesis', 'pathwayGeneCount': '57'}],
                               [{'pathwayId': 'WP5153', 'pathwayLabel': 'N-glycan biosynthesis', 'pathwayGeneCount': '57'}],
                               [{'pathwayId': nan, 'pathwayLabel': nan, 'pathwayGeneCount': nan}]])
    expected_data.name = 'WikiPathways'

    pd.testing.assert_series_equal(obtained_data["WikiPathways"], expected_data)


@pytest.fixture(scope="module")
def bridgedb_dataframe():
    return pd.DataFrame({
        'identifier': ['ALG14', 'ALG2', 'CHRNA1'],
        'identifier.source': ['HGNC', 'HGNC', 'HGNC'],
        'target': ['199857', '85365', '1134'],
        'target.source': ['NCBI Gene', 'NCBI Gene', 'NCBI Gene']
    })

