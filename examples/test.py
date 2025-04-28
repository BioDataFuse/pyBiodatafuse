# Import modules
import pickle

import pandas as pd
from IPython.display import Image, display

from pyBiodatafuse import id_mapper
from pyBiodatafuse.annotators import (
    bgee,
    disgenet,
    minerva,
    molmedb,
    opentargets,
    pubchem,
    stringdb,
    wikipathways,
)
from pyBiodatafuse.constants import DISGENET_DISEASE_COL
from pyBiodatafuse.graph import generator
from pyBiodatafuse.graph.rdf import BDFGraph
from pyBiodatafuse.utils import (
    combine_sources,
    create_harmonized_input_file,
    create_or_append_to_metadata,
)

genes_of_interest = """PIK3CA
BRCA2
TP53"""


gene_list = genes_of_interest.split("\n")
data_input = pd.DataFrame(gene_list, columns=["identifier"])
data_input.head()
bridgedb_df, bridgedb_metadata = id_mapper.bridgedb_xref(
    identifiers=data_input,
    input_species="Human",
    input_datasource="HGNC", # Commented out
    output_datasource="All",
)

wikipathways_df, wikipathways_metadata = wikipathways.get_gene_wikipathways(bridgedb_df=bridgedb_df, query_interactions=True)
wikipathways_df.to_csv("wp.csv", index=False)

combined_df = combine_sources(bridgedb_df, [wikipathways_df])
combined_metadata = create_or_append_to_metadata(bridgedb_metadata, [wikipathways_metadata])

pygraph = generator.save_graph(
    combined_df=combined_df,
    combined_metadata=combined_metadata,
#    disease_compound=opentargets_df,
    graph_name="examples",
    graph_dir="./data",
)
