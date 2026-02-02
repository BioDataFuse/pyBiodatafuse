import json
import os
from unittest.mock import patch

import pandas as pd
from IPython.display import Image, display

import pyBiodatafuse.annotators.aopwiki as aopwiki
import pyBiodatafuse.constants as Cons
from pyBiodatafuse import id_mapper

aopwiki_path = "tests/annotators/data/aop_mock_res_simple.json"

gene_list = ["4193"]

data_input = pd.DataFrame(gene_list, columns=["identifier"])

bridgedb_dataframe, bridgedb_metadata = id_mapper.bridgedb_xref(
    identifiers=data_input,
    input_species="Human",
    input_datasource="NCBI Gene",
    output_datasource="Ensembl",
)
aopwiki_simple_df, aopwiki_simple_metadata = aopwiki.get_aops(
    bridgedb_df=bridgedb_dataframe, pathway=False
)
aopwiki_simple_df.to_json(aopwiki_path)

aopwiki_df, aopwiki_metadata = aopwiki.get_aops(bridgedb_df=bridgedb_dataframe, pathway=True)
aopwiki_df.to_json(aopwiki_path.replace("_simple", ""))