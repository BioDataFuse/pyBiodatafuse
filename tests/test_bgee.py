# coding: utf-8

"""Python file for testing the Bgee workflo."""

import pandas as pd
from pyBiodatafuse import id_mapper
from pyBiodatafuse.annotators import bgee


def main():
    """Run the workflow."""

    genes_of_interest = """AGRN"""
    gene_list = genes_of_interest.split("\n")

    anatomical_entities_of_interest = """
respiratory system
heart
brain
"""

    anatomical_entities_list = anatomical_entities_of_interest.split("\n")
    anatomical_entities_list = [anatEntity.strip()  for anatEntity in anatomical_entities_list if anatEntity.strip() != '']

    data_input = pd.DataFrame(gene_list, columns=["identifier"])

    bridgdb_df, bridgdb_metadata = id_mapper.bridgedb_xref(
        identifiers=data_input["identifier"],
        input_species="Human",
        input_datasource="HGNC",
        output_datasource="All",
    )

    anatomical_entities_df = pd.DataFrame(anatomical_entities_list, columns = ["AnatomicalEntityNames"])

    d = bgee.get_gene_literature(bridgdb_df, anatomical_entities_df)
    print(d)

if __name__ == "__main__":
    main()
