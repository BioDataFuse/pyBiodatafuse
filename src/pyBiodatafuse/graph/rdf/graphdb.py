import tempfile

import pandas as pd
import requests
from requests.exceptions import HTTPError


class GraphDBManager:
    """
    A class to manage GraphDB repositories via REST API.
    """

    @staticmethod
    def create_repository(
        base_url: str,
        repository_name: str = "default",
        username: str = None,
        password: str = None,
    ):
        """
        Create a new repository in a GraphDB instance.

        :param base_url: Base URL of the GraphDB instance (e.g., http://localhost:7200).
        :param repository_name: Name of the repository to create. Defaults to "default".
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :raises HTTPError: If the HTTP request fails.
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
        files = {"config": (temp_file_path, open(temp_file_path, "rb"), "text/turtle")}

        # Send the POST request with the configuration file
        response = requests.post(
            url,
            files=files,
            auth=auth,
            timeout=10,
        )
        response.raise_for_status()

    @staticmethod
    def list_repositories(base_url: str, username: str = None, password: str = None):
        """
        List all repositories in the GraphDB instance.

        :param base_url: Base URL of the GraphDB instance.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :return: List of repositories as JSON.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories"
        auth = (username, password) if username and password else None
        response = requests.get(url, auth=auth, headers={"Accept": "application/json"})
        response.raise_for_status()
        return response.json()

    @staticmethod
    def get_repository_info(
        base_url: str, repository_id: str, username: str = None, password: str = None
    ):
        """
        Retrieve detailed information about a specific repository.

        :param base_url: Base URL of the GraphDB instance.
        :param repository_id: ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :return: Repository information as JSON.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}"
        auth = (username, password) if username and password else None
        response = requests.get(url, auth=auth, headers={"Accept": "application/json"})
        response.raise_for_status()
        return response.json()

    @staticmethod
    def count_triples(
        base_url: str, repository_id: str, username: str = None, password: str = None
    ):
        """
        Count the number of triples in a repository.

        :param base_url: Base URL of the GraphDB instance.
        :param repository_id: ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        :return: Number of triples in the repository.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}/size"
        auth = (username, password) if username and password else None
        response = requests.get(url, auth=auth)
        response.raise_for_status()
        return dict(response.json())

    @staticmethod
    def restart_repository(
        base_url: str, repository_id: str, username: str = None, password: str = None
    ):
        """
        Restart a specific repository.

        :param base_url: Base URL of the GraphDB instance.
        :param repository_id: ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}/restart"
        auth = (username, password) if username and password else None
        response = requests.post(url, auth=auth)
        response.raise_for_status()

    @staticmethod
    def delete_repository(
        base_url: str, repository_id: str, username: str = None, password: str = None
    ):
        """
        Delete a specific repository.

        :param base_url: Base URL of the GraphDB instance.
        :param repository_id: ID of the repository.
        :param username: Optional username for authentication.
        :param password: Optional password for authentication.
        """
        url = f"{base_url.rstrip('/')}/rest/repositories/{repository_id}"
        auth = (username, password) if username and password else None
        response = requests.delete(url, auth=auth)
        response.raise_for_status()

    @staticmethod
    def upload_to_graphdb(
        base_url: str,
        repository_id: str,
        username: str,
        password: str,
        bdf_graph,
        format: str = "turtle",
    ):
        """
        Upload a BDFGraph to the specified GraphDB repository.

        :param base_url: Base URL of the GraphDB instance (e.g., http://localhost:7200).
        :param repository_id: ID of the repository to upload the graph to.
        :param username: Username for authentication.
        :param password: Password for authentication.
        :param bdf_graph: The BDFGraph object to upload.
        :param format: Format of the RDF data (e.g., "turtle", "rdfxml").
        :raises HTTPError: If the HTTP request fails.
        """
        url = f"{base_url.rstrip('/')}/repositories/{repository_id}/statements"
        auth = (username, password)

        # Serialize the BDFGraph to the specified format
        rdf_data = bdf_graph.serialize(format=format)

        # Map format to Content-Type
        content_type_map = {
            "turtle": "text/turtle",
            "rdfxml": "application/rdf+xml",
            "ntriples": "application/n-triples",
            "jsonld": "application/ld+json",
        }
        content_type = content_type_map.get(format, "text/turtle")

        headers = {"Content-Type": content_type, "Accept": "application/json"}

        # Send the POST request to upload the RDF data
        response = requests.post(
            url,
            data=rdf_data,
            headers=headers,
            auth=auth,
            timeout=10,
        )
        response.raise_for_status()

    @staticmethod
    def query_graphdb(base_url, repository_name, username, password, query, response_format="json"):
        """
        Sends a SPARQL SELECT query to the specified GraphDB repository.

        Args:
            base_url (str): The base URL of the GraphDB instance.
            repository_name (str): The name of the repository to query.
            username (str): The username for authentication.
            password (str): The password for authentication.
            query (str): The SPARQL query to execute.
            response_format (str): The desired response format ("json" or "dataframe").

        Returns:
            dict or pandas.DataFrame: The query results in JSON format or as a DataFrame.

        Raises:
            Exception: If the query fails or the response is invalid.
        """

        endpoint = f"{base_url.rstrip('/')}/repositories/{repository_name}"
        headers = {"Accept": "application/sparql-results+json"}
        auth = (username, password)
        params = {"query": query}

        response = requests.get(endpoint, headers=headers, auth=auth, params=params)

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
            raise Exception(
                f"Query failed with status code {response.status_code}: {response.text}"
            )
