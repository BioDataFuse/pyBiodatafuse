# coding: utf-8

"""Python file for querying the DisGeNet database (https://www.disgenet.org/home/)."""

import pandas as pd
import requests
import datetime

from pyBiodatafuse.constants import DATA_DIR
from pyBiodatafuse.utils import (
    get_identifier_of_interest,
    collapse_data_sources,
)


def disgenetAnnotator(
    bridgedb_df: pd.DataFrame,
    api_key: str = "0209751bfa7b6a981a8f5fb5f062313067ecd36c",
    params: dict = {"source": "CURATED", "format": "json"},
) -> pd.DataFrame:
    """Query gene-disease associations from DisGeNET.

    @param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    @param api_key: DisGeNET API key (more details can be found at https://www.disgenet.org/api/#/Authorization)
    @param params: dictionary of parameters to be passed to the DisGeNET API (more details can be found at https://www.disgenet.org/api/#/gene)

    Usage example:
    >> api_key = "YOUR_KEY_HERE"
    >> params = {'source': 'CURATED', 'format': 'json'} # only curated data
    >> disgenetAnnotator(api_key, params)
    """

    # Extract the "target" values and join them into a single string separated by commas
    data_df = get_identifier_of_interest(bridgedb_df, "NCBI Gene")
    disgenet_input = ",".join(data_df["target"])

    # Set the DisGeNET API
    api_host = "https://www.disgenet.org/api"
    s = requests.Session()

    if not api_key:
        raise ValueError("Please provide a DisGeNET API key")

    # Add the API key to the requests headers
    s.headers.update({"Authorization": "Bearer %s" % api_key})
    # Record the start time
    start_time = datetime.datetime.now()

    # Split the targets into chunks of 99 or fewer
    targets_list = disgenet_input.split(",")
    chunk_size = 99
    chunks = [targets_list[i : i + chunk_size] for i in range(0, len(targets_list), chunk_size)]

    disgenet_output = []

    for chunk in chunks:
        # Join the chunked targets into a comma-separated string
        chunked_input = ",".join(chunk)
        # Get all the diseases associated with genes for the current chunk
        gda_response = s.get(f"{api_host}/gda/gene/{chunked_input}", params=params)
        chunk_output = gda_response.json()
        disgenet_output.extend(chunk_output)

    # Record the end time
    end_time = datetime.datetime.now()

    # Convert disgenet_output to a DataFrame
    disgenet_df = pd.DataFrame(disgenet_output)
    if "geneid" not in disgenet_df:
        return pd.DataFrame()
    else:
        # Drop the uniprotid column
        disgenet_df.drop("uniprotid", axis=1, inplace=True)
        # Add DisGeNET output as a new column to BridgeDb file
        disgenet_df.rename(columns={"geneid": "target", "gene_symbol": "identifier"}, inplace=True)
        disgenet_df["target"] = disgenet_df["target"].values.astype(str)

        # Save the output in a CSV file
        # disgenet_df.to_csv(disgenet_file_path, index=False)
        # print(f"The query output is saved as {disgenet_file_path}")

        selected_columns = [
            "gene_dsi",
            "gene_dpi",
            "gene_pli",
            "protein_class",
            "protein_class_name",
            "diseaseid",
            "disease_name",
            "disease_class",
            "disease_class_name",
            "disease_type",
            "disease_semantic_type",
            "score",
            "ei",
            "el",
            "year_initial",
            "year_final",
            "source",
        ]

        # Merge the two DataFrames based on 'geneid', 'gene_symbol', 'identifier', and 'target'
        merged_df = collapse_data_sources(
            data_df=data_df,
            source_namespace="NCBI Gene",
            target_df=disgenet_df,
            common_cols=["target", "identifier"],
            target_specific_cols=selected_columns,
            col_name="DisGeNET",
        )

        # Save the output in a CSV file
        # merged_df.to_csv(output_file_path, index=False)

        # print(f"The complete merged output is saved as {output_file_path}")

        """Metdata details"""
        # Get the current date and time
        current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Calculate the time elapsed
        time_elapsed = str(end_time - start_time)
        # Add version to metadata file
        version_response = s.get(api_host + "/version/")
        disgenet_version = version_response.json()
        # Add the datasource, query, query time, and the date to metadata
        disgenet_metadata = {
            "datasource": "DisGeNET",
            "metadata": {"source_version": disgenet_version},
            "query": {
                "size": len(disgenet_input.split(",")),
                "input_type": "NCBI Gene",
                "time": time_elapsed,
                "date": current_date,
                "url": gda_response.request.url,
            },
        }

        if s:
            s.close()

        return merged_df, disgenet_metadata