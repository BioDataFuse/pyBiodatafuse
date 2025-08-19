"""
update_datasources.py

This script updates the `resources/datasources.csv` file with pattern data fetched from the Bioregistry API.

Usage:
    python update_datasources.py

Requirements:
    - Python 3.9-3.10
    - pandas
    - bioregistry

Notes:
    - The Bioregistry API is accessed using the URL template:
      "https://bioregistry.io/api/registry/{prefix}?format=json".
    - If a prefix does not have a corresponding pattern in the API, an empty string is used.
    - This script is intended for occasional use to update the resource file.

Raises:
    - FileNotFoundError: If the `resources/datasources.csv` file is not found.
"""

import bioregistry
import pandas as pd


def update_datasources_file():
    """
    Updates the datasources.csv file with pattern data fetched from the Bioregistry API.

    Note:
        - The Bioregistry API is accessed via :func:`bioregistry.get_pattern`
        - If a prefix does not have a corresponding pattern in the API, an empty string is used.

    Raises:
        FileNotFoundError: If the datasources.csv file is not found.
    """

    # Access the datasources.csv file using importlib.resources
    csv_path = "resources/datasources.csv"
 
    df = pd.read_csv(csv_path)

    df["pattern"] = df["prefix"].map(bioregistry.get_pattern, na_action="ignore")

    # Write the updated data back to the CSV file
    df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    update_datasources_file()
