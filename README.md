<!--
<p align="center">
  <img src="https://github.com/BioDataFuse/pyBiodatafuse/raw/main/docs/source/logo.png" height="150">
</p>
-->

<h1 align="center">
  pyBioDataFuse
</h1>

<p align="center">
<!--     <a href="https://github.com/BioDataFuse/pyBiodatafuse/actions/workflows/tests.yml">
        <img alt="Tests" src="https://github.com/BioDataFuse/pyBiodatafuse/workflows/Tests/badge.svg" />
    </a> -->
    <a href="https://pypi.org/project/pyBiodatafuse">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/pyBiodatafuse" />
    </a>
    <a href="https://pypi.org/project/pyBiodatafuse">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/pyBiodatafuse" />
    </a>
    <a href="https://github.com/BioDataFuse/pyBiodatafuse/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/pyBiodatafuse" />
    </a>
<!--     <a href='https://pyBiodatafuse.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/pyBiodatafuse/badge/?version=latest' alt='Documentation Status' />
    </a> -->
    <a href="https://codecov.io/gh/BioDataFuse/pyBiodatafuse/branch/main">
        <img src="https://codecov.io/gh/BioDataFuse/pyBiodatafuse/branch/main/graph/badge.svg" alt="Codecov status" />
    </a>  
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img alt="Cookiecutter template from @cthoyt" src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" /> 
    </a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' />
    </a>
    <a href="https://github.com/BioDataFuse/pyBiodatafuse/blob/main/.github/CODE_OF_CONDUCT.md">
        <img src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" alt="Contributor Covenant"/>
    </a>
</p>


## 💪 Getting Started

> We introduce BioDataFuse, a query-based Python tool for seamless integration of biomedical databases. BioDataFuse establishes a modular framework for efficient data wrangling, enabling context-specific knowledge graph creation and supporting graph-based analyses. With a user-friendly interface, it enables users to dynamically create knowledge graphs from their input data. Supported by a robust Python package, pyBiodatafuse, this tool excels in data harmonization, aggregating diverse sources through modular queries. Moreover, BioDataFuse provides plugin capabilities for Cytoscape and Neo4j, allowing local graph hosting. Ongoing refinements enhance the graph utility through tasks like link prediction, making BioDataFuse a versatile solution for efficient and effective biological data integration.

To know more about the package, read our documentation [here](https://pybiodatafuse.readthedocs.io/en/latest/index.html).

## Creating your own graph

To generate your own graph, check out our tutorial notebook [in examples](examples).

<!-- ### Command Line Interface

The pyBiodatafuse command line tool is automatically installed. It can
be used from the shell with the `--help` flag to show all subcommands:

```shell
$ pyBiodatafuse --help
```

> TODO show the most useful thing the CLI does! The CLI will have documentation auto-generated
> by `sphinx`. -->

We support exporting of the graphs in Cytoscape, Neo4J and GraphDB. You can use the following functions:

```python
# on neo4j
neo4j.load_graph(pygraph, uri="bolt://localhost:7687", username="YOUR_USERNAME", password="YOUR_PASSWORD")  # change username and password

# on cytoscape
cytoscape.load_graph(pygraph, network_name="YOUR_CUSTOM_NAME")

# rdf ttl files
bdf = BDFGraph(
    base_uri="https://biodatafuse.org/YOUR_CUSTOM_NAME/",
    version_iri="https://biodatafuse.org/example/YOUR_CUSTOM_NAME.ttl",
    orcid="YOUR_ORCID",
    author="YOUR_NAME",
)

bdf.generate_rdf(combined_df, combined_metadata)  # Generate the RDF from the (meta)data files from the example runs
bdf.serialize(
    "YOUR_CUSTOM_NAME.ttl",
    format="ttl",
)
```

## 🚀 Installation

The most recent release can be installed from
[PyPI](https://pypi.org/project/pyBiodatafuse/) with:

```shell
$ pip install pyBiodatafuse
```

The most recent code and data can be installed directly from GitHub with:

```bash
$ pip install git+https://github.com/BioDataFuse/pyBiodatafuse.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/BioDataFuse/pyBiodatafuse/blob/master/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

### ⚖️ License

The code in this package is licensed under the MIT License.


### 📖 Citation

The work was started as part of the [Elixir BioHackathon 2023](https://github.com/elixir-europe/biohackathon-projects-2023/tree/main/17) integrating and bringing together multiple Core Data Resources together.
> Gadiya, Y., Ammar, A., Willighagen, E., Martinat, D., Sima, A. C., Balci, H., & Abbassi Daloii, T. (2023). BioHackEU23 report: Extending interoperability of experimental data using modular queries across biomedical resources. BioHackrXiv Preprints. https://doi.org/10.37044/osf.io/mhsqp

<!--
### 🎁 Support

This project has been supported by the following organizations (in alphabetical order):

- [Harvard Program in Therapeutic Science - Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/)

-->

<!--
### 💰 Funding

This project has been supported by the following grants:

| Funding Body                                             | Program                                                                                                                       | Grant           |
|----------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|-----------------|
| DARPA                                                    | [Automating Scientific Knowledge Extraction (ASKE)](https://www.darpa.mil/program/automating-scientific-knowledge-extraction) | HR00111990009   |
-->

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>

The final section of the README is for if you want to get involved by making a code contribution.

### Development Installation

To install in development mode, use the following:

```bash
$ git clone git+https://github.com/BioDataFuse/pyBiodatafuse.git
$ cd pyBiodatafuse
$ pip install -e .
```

### 🥼 Testing

After cloning the repository and installing `tox` with `pip install tox`, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
$ tox
```

Additionally, these tests are automatically re-run with each commit in a [GitHub Action](https://github.com/BioDataFuse/pyBiodatafuse/actions?query=workflow%3ATests).

### 📖 Building the Documentation

The documentation can be built locally using the following:

```shell
$ git clone git+https://github.com/BioDataFuse/pyBiodatafuse.git
$ cd pyBiodatafuse
$ tox -e docs
$ open docs/build/html/index.html
``` 

The documentation automatically installs the package as well as the `docs`
extra specified in the [`setup.cfg`](setup.cfg). `sphinx` plugins
like `texext` can be added there. Additionally, they need to be added to the
`extensions` list in [`docs/source/conf.py`](docs/source/conf.py).

### 📦 Making a Release

After installing the package in development mode and installing
`tox` with `pip install tox`, the commands for making a new release are contained within the `finish` environment
in `tox.ini`. Run the following from the shell:

```shell
$ tox -e finish
```

This script does the following:

1. Uses [Bump2Version](https://github.com/c4urself/bump2version) to switch the version number in the `setup.cfg`,
   `src/pyBiodatafuse/version.py`, and [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using [`build`](https://github.com/pypa/build)
3. Uploads to PyPI using [`twine`](https://github.com/pypa/twine). Be sure to have a `.pypirc` file configured to avoid the need for manual input at this
   step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use `tox -e bumpversion -- minor` after.
</details>
