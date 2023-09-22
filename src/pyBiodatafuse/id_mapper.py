# coding: utf-8

"""Python file for mapping identifiers using BridgeDb."""

import requests
import pandas as pd
import datetime
import csv
import os

from .constants import BRIDGEDB_DIR, DATA_DIR
from .utils import create_or_append_to_metadata


def bridgedb_Xref(
    identifiers: pd.DataFrame,
    inputSpecies: str = None,
    inputDatasource: str = None,
    outputDatasource: list = None,
) -> pd.DataFrame:
    """Map input list using BridgeDb.

    @param identifiers: a dataframe with one column called identifier (the output of data_loader.py)
    @param inputSpecies: specify the species, for now only human would be supported
    @param inputDatasource: type of input identifier (more details can be found at https://www.bridgedb.org/pages/system-codes.html)
    @param outputDatasource: specify which type of identifiers you want to map your inputidentifiers with (more details can be found at https://www.bridgedb.org/pages/system-codes.html)

    Usage example:
    >> identifiers = path_to_output_of_data_loader
    >> inputSpecies = "Human"
    >> inputDatasource = "HGNC"
    >> outputDatasource = ["All"] or one or multiple source ["Ensembl"]
    """
    if inputSpecies is None:
        inputSpecies = "Human"

    if not inputSpecies:
        raise ValueError("Please provide the identifier datasource, e.g. HGNC")

    if outputDatasource is None:
        outputDatasource = [
            "RefSeq",
            "WikiGenes",
            "OMIM",
            "Uniprot-TrEMBL",
            "NCBI Gene",
            "Ensembl",
            "HGNC Accession Number",
            "PDB",
            "HGNC",
        ]

    data_sources = pd.read_csv(f"data/bridgedb/datasources.tsv", sep=",")
    input_source = data_sources.loc[
        data_sources["source"] == inputDatasource, "systemCode"
    ].iloc[0]

    if len(identifiers) != 0:
        post_con = (
            "\n".join(
                [f"{identifier}\t{input_source}" for identifier in identifiers]
            )
            + f"\t{input_source}\n"
        )

        # Setting up the query url
        url = "https://webservice.bridgedb.org"
        query_link = f"{url}/{inputSpecies}/xrefsBatch"

        # Record the start time
        start_time = datetime.datetime.now()

        # Getting the response to the query
        try:
            s = requests.post(url=query_link, data=post_con.encode())
        except Exception as e:
            print("Error:", e)
            return None

        # Extracting the content in the raw text format
        out = s.content.decode()
        lines = out.split("\n")

        # Record the end time
        end_time = datetime.datetime.now()

        # Processing each line and splitting values
        parsed_results = []

        for line in lines:
            if line:
                parts = line.split("\t")
                identifier = parts[0]
                identifier_source = parts[1]
                targets = parts[2].split(",")

            for target in targets:
                target_parts = target.split(":")
                target_source = target_parts[0]
                target_id = ":".join(target_parts[1:])

                parsed_results.append(
                    [identifier, identifier_source, target_id, target_source]
                )

        # Create a DataFrame
        bridgedb = pd.DataFrame(
            parsed_results,
            columns=["identifier", "identifier.source", "target", "target.source"],
        )

        # Replace 'target.source' values with complete source names from 'data_sources'
        bridgedb["target.source"] = bridgedb["target.source"].map(
            data_sources.set_index("systemCode")["source"]
        )

        # Subset based on the outputDatasource
        if not outputDatasource == "All":
            bridgedb = bridgedb[bridgedb["target.source"].isin(outputDatasource)]

        bridgedb = bridgedb.drop_duplicates()

        # Save the output in a CSV file
        # Specify the CSV file path for the BridgeDb output
        bridgedb_file_path = os.path.join(BRIDGEDB_DIR, "bridgedb_output.csv")

        bridgedb.to_csv(bridgedb_file_path, index=False)

        print(f"DataFrame saved to {bridgedb_file_path}")

        """Metdata details"""
        # Get the current date and time
        current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Calculate the time elapsed
        time_elapsed = str(end_time - start_time)
        # Add BridgeDb version to metadata file
        version_response = requests.get(url=f"{url}/config")

        # Check if the request was successful (status code 200)
        if version_response.status_code == 200:
            # Initialize an empty dictionary to store the data
            bridgedb_version = {}
            # Split the response content into lines and create a CSV reader
            lines = version_response.text.strip().split("\n")
            csv_reader = csv.reader(lines, delimiter="\t")

            # Iterate over the rows in the CSV and populate the dictionary
            for row in csv_reader:
                if len(row) == 2:
                    key, value = row
                    bridgedb_version[key] = value
        else:
            print(
                f"Failed to retrieve data. Status code: {version_response.status_code}"
            )

        # Add datasource version to metadata file
        datasource_response = requests.get(url=f"{url}/{inputSpecies}/properties")

        # Check if the request was successful (status code 200)
        if datasource_response.status_code == 200:
            datasource_version = datasource_response.text.strip().split("\n")
            datasource_version = [
                line.replace("\t", ": ") for line in datasource_version
            ]
        else:
            print(
                f"Failed to retrieve data. Status code: {datasource_response.status_code}"
            )

        # Add the datasource, query, query time, and the date to metadata
        bridgedb_metadata = {
            "datasource": "BridgeDb",
            "metadata": {
                "source_version": bridgedb_version,
                "data_version": datasource_version,
            },
            "query": {
                "size": len(identifiers),
                "input_type": inputDatasource,
                "time": time_elapsed,
                "date": current_date,
                "url": s.url,
                "request_string": f"{post_con.encode().decode('utf-8')}",
            },
        }

        create_or_append_to_metadata(
            bridgedb_metadata
        )  # Call the function from the metadata module

        print(f"The query metadata is saved in {DATA_DIR}\metadata.json")

        return bridgedb
