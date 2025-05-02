"""
update_datasources.py

This script updates the `resources/datasources.csv` file with pattern data fetched from the Bioregistry API.

Usage:
    python update_datasources.py

Requirements:
    - Python 3.9-3.10
    - pandas
    - requests

Notes:
    - The Bioregistry API is accessed using the URL template:
      "https://bioregistry.io/api/registry/{prefix}?format=json".
    - If a prefix does not have a corresponding pattern in the API, an empty string is used.
    - This script is intended for occasional use to update the resource file.

Raises:
    - FileNotFoundError: If the `resources/datasources.csv` file is not found.
    - requests.exceptions.RequestException: If there is an issue with the API request.
"""

import pandas as pd
import requests


def update_datasources_file():
    """
    Updates the datasources.csv file with pattern data fetched from the Bioregistry API.

    Note:
        - The Bioregistry API is accessed using the URL template:
          "https://bioregistry.io/api/registry/{prefix}?format=json".
        - If a prefix does not have a corresponding pattern in the API, an empty string is used.

    Raises:
        FileNotFoundError: If the datasources.csv file is not found.
        requests.exceptions.RequestException: If there is an issue with the API request.
    """
    url_bioregistry = "https://bioregistry.io/api/registry/{}?format=json"

    # Access the datasources.csv file using importlib.resources
    csv_path = "resources/datasources.csv"
    with open(csv_path) as f:
        # Read the CSV file using pandas
        df = pd.read_csv(f)

        # Add a new column "pattern" and fetch data from Bioregistry
        def fetch_pattern(prefix):
            if pd.notna(prefix):
                response = requests.get(url_bioregistry.format(prefix))
                if response.status_code == 200:
                    print(response.json().get("pattern", ""))
                    return response.json().get("pattern", "")
                else:
                    print("bad request for ", prefix)
            return ""

        df["pattern"] = df["prefix"].map(fetch_pattern)

        # Write the updated data back to the CSV file
        df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    update_datasources_file()
