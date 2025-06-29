"""

Module for managing GraphDB repositories via REST API.

This module provides functions to create, list, and manage GraphDB repositories,
as well as upload RDF data and execute SPARQL queries.
"""

import tempfile
from typing import Optional

import pandas as pd
import requests
from requests.exceptions import HTTPError


class GraphDBManager:
    """A class to manage GraphDB repositories via REST API."""

    @staticmethod
    def create_repository(
        base_url: str,
        repository_name: str = "default",
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        """
        Create a new repository in GraphDB.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_name: The name of the repository to create.
        :param username: The username for authentication.
        :param password: The password for authentication.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories"
        auth = (username, password) if username and password else None

        # Prepare the repository configuration file
        config_content = f"""
        @prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>.
        @prefix rep: <http://www.openrdf.org/config/repository#>.
        @prefix sr: <http://www.openrdf.org/config/repository/sail#>.
        @prefix sail: <http://www.openrdf.org/config/sail#>.
        @prefix graphdb: <http://www.ontotext.com/config/graphdb#>.

        [] a rep:Repository ;
            rep:repositoryID "{repository_name}" ;
            rdfs:label "{repository_name.capitalize()}" ;
            rep:repositoryImpl [
            rep:repositoryType "graphdb:SailRepository" ;
            sr:sailImpl [
                sail:sailType "graphdb:Sail" ;

                graphdb:read-only "false" ;

                # Inference and Validation
                graphdb:ruleset "rdfsplus-optimized" ;
                graphdb:disable-sameAs "true" ;
                graphdb:check-for-inconsistencies "false" ;

                # Indexing
                graphdb:entity-id-size "32" ;
                graphdb:enable-context-index "false" ;
                graphdb:enablePredicateList "true" ;
                graphdb:enable-fts-index "false" ;
                graphdb:fts-indexes ("default" "iri") ;
                graphdb:fts-string-literals-index "default" ;
                graphdb:fts-iris-index "none" ;

                # Queries and Updates
                graphdb:query-timeout "0" ;
                graphdb:throw-QueryEvaluationException-on-timeout "false" ;
                graphdb:query-limit-results "0" ;

                # Settable in the file but otherwise hidden in the UI and in the RDF4J console
                graphdb:base-URL "http://example.org/owlim#" ;
                graphdb:defaultNS "" ;
                graphdb:imports "" ;
                graphdb:repository-type "file-repository" ;
                graphdb:storage-folder "storage" ;
                graphdb:entity-index-size "10000000" ;
                graphdb:in-memory-literal-properties "true" ;
                graphdb:enable-literal-index "true" ;
            ]
            ].
        """.format(
            repository_name=repository_name
        )
        # Save the configuration content to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".ttl") as temp_file:
            temp_file.write(config_content.encode("utf-8"))
            temp_file_path = temp_file.name

        # Pass the temporary file path to the config argument
        with open(temp_file_path, "rb") as config_file:
            files = {"config": (temp_file_path, config_file, "text/turtle")}

            # Send the POST request with the configuration file
            response = requests.post(
                url,
                files=files,
                auth=auth,
                timeout=10,
            )
        response.raise_for_status()

    @staticmethod
    def list_repositories(
        base_url: str, username: Optional[str] = None, password: Optional[str] = None
    ):
        """
        List all repositories in the GraphDB instance.

        :param base_url: The base URL of the GraphDB instance.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :return: List of repositories as JSON.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories"
        auth = (username, password) if username and password else None
        response = requests.get(url, auth=auth, headers={"Accept": "application/json"}, timeout=10)
        response.raise_for_status()
        return response.json()

    @staticmethod
    def get_repository_info(
        base_url: str,
        repository_id: str,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        """
        Retrieve detailed information about a specific repository.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_id: The ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :return: Repository information as JSON.
        :raises HTTPError: If the request to retrieve repository info fails.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}"
        auth = (username, password) if username and password else None
        response = requests.get(url, auth=auth, headers={"Accept": "application/json"}, timeout=10)
        try:
            response.raise_for_status()
        except HTTPError as e:
            raise HTTPError(
                f"Failed to retrieve repository info: {e.response.status_code} - {e.response.text}"
            ) from e
        return response.json()

    @staticmethod
    def count_triples(
        base_url: str,
        repository_id: str,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        """
        Count the number of triples in a repository.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_id: The ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :return: Number of triples in the repository.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}/size"
        auth = (username, password) if username and password else None
        response = requests.get(url, auth=auth, timeout=10)
        response.raise_for_status()
        return dict(response.json())

    @staticmethod
    def restart_repository(
        base_url: str,
        repository_id: str,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        """
        Restart a specific repository.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_id: The ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}/restart"
        auth = (username, password) if username and password else None
        response = requests.post(url, auth=auth, timeout=10)
        response.raise_for_status()

    @staticmethod
    def delete_repository(
        base_url: str,
        repository_id: str,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        """
        Delete a specific repository.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_id: The ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}"
        auth = (username, password) if username and password else None
        response = requests.delete(url, auth=auth, timeout=10)
        response.raise_for_status()

    @staticmethod
    def upload_to_graphdb(
        base_url: str,
        repository_id: str,
        username: str,
        password: str,
        bdf_graph,
        file_format: str = "turtle",
    ):
        """
        Upload an RDF graph to a GraphDB repository.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_id: The ID of the repository to upload the graph to.
        :param username: The username for authentication.
        :param password: The password for authentication.
        :param bdf_graph: The RDF graph to upload.
        :param file_format: The format of the RDF graph (default is "turtle").
        :raises HTTPError: If the request to execute the query fails.
        """
        url = f"{base_url.rstrip('/')}/repositories/{repository_id}/statements"
        auth = (username, password)

        # Serialize the BDFGraph to the specified format
        rdf_data = bdf_graph.serialize(format=file_format)
        if isinstance(rdf_data, bytes):
            rdf_data = rdf_data.decode("utf-8")

        # Map format to Content-Type
        content_type_map = {
            "turtle": "text/turtle",
            "rdfxml": "application/rdf+xml",
            "ntriples": "application/n-triples",
            "jsonld": "application/ld+json",
        }
        content_type = content_type_map.get(file_format, "text/turtle")

        headers = {"Content-Type": content_type, "Accept": "application/json"}

        # Send the POST request to upload the RDF data
        response = requests.post(
            url,
            data=rdf_data,
            headers=headers,
            auth=auth,
            timeout=300,
        )
        if not response.ok:
            # Print error details for debugging
            print("GraphDB upload error:", response.status_code, response.reason)
            print("Response content:", response.text)
            print("RDF data preview:", rdf_data[:500])
            raise HTTPError(
                f"GraphDB upload failed: {response.status_code} {response.reason}\n{response.text}"
            )

    @staticmethod
    def query_graphdb(
        base_url: str,
        repository_name: str,
        username: str,
        password: str,
        query: str,
        response_format: str = "json",
    ):
        """
        Execute a SPARQL query on a GraphDB repository.

        :param base_url: The base URL of the GraphDB instance.
        :param repository_name: The name of the repository to query.
        :param username: The username for authentication.
        :param password: The password for authentication.
        :param query: The SPARQL query to execute.
        :param response_format: The format of the query response (default is "json").
        :return: Query results as a dictionary or pandas DataFrame.
        :raises HTTPError: If the request to execute the query fails.
        """
        endpoint = f"{base_url.rstrip('/')}/repositories/{repository_name}"
        headers = {"Accept": "application/sparql-results+json"}
        auth = (username, password)
        params = {"query": query}

        response = requests.get(endpoint, headers=headers, auth=auth, params=params, timeout=10)

        if response.status_code == 200:
            results = response.json()
            if response_format == "dataframe":
                # Convert SPARQL results to a pandas DataFrame
                columns = results["head"]["vars"]
                rows = [
                    {col: binding.get(col, {}).get("value", None) for col in columns}
                    for binding in results["results"]["bindings"]
                ]
                return pd.DataFrame(rows, columns=columns)
            return results
        else:
            raise HTTPError(
                f"Query failed with status code {response.status_code}: {response.text}"
            )
