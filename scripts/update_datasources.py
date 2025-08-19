# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "bioregistry",
#     "pandas",
# ]
# ///
#

"""
This script updates the `resources/datasources.csv` file with pattern data fetched from the Bioregistry.

This script can be run with ``uv``, which automatically creates
an environment with the appropriate dependencies using:

.. code-block:: console

    $ git clone https://github.com/BioDataFuse/pyBiodatafuse
    $ cd pyBiodatafuse
    $ uv run update_datasources.py

Notes:
    - Patterns are accesed via the Bioregistry package's :func:`bioregistry.get_pattern`
    - If a prefix does not have a corresponding pattern in the Bioregistry, an empty string is used.
    - This script is intended for occasional use to update the resource file.

Raises:
    - FileNotFoundError: If the `resources/datasources.csv` file is not found.
"""

import bioregistry
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent.resolve()
PATH = ROOT.joinpath("src", "pyBiodatafuse", "resources", "datasources.csv")


def update_datasources_file():
    """
    Updates the datasources.csv file with pattern data fetched from the Bioregistry.

    Note:
        - Patterns are accesed via the Bioregistry package's :func:`bioregistry.get_pattern`
        - If a prefix does not have a corresponding pattern in the Bioregistry, an empty string is used.

    Raises:
        FileNotFoundError: If the datasources.csv file is not found.
    """
    df = pd.read_csv(PATH)

    df["pattern"] = df["prefix"].map(bioregistry.get_pattern, na_action="ignore")

    # Write the updated data back to the CSV file
    df.to_csv(PATH, index=False)


if __name__ == "__main__":
    update_datasources_file()
