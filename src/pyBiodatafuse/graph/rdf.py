# coding: utf-8

"""Python module to produce an RDF graph from the property table."""

import logging

import pandas as pd
from bioregistry import get_iri, normalize_curie
from rdflib import Graph, Literal, Namespace, URIRef
from rdflib.namespace import RDF, RDFS, SKOS, XSD, DC
import numpy as np
from pyBiodatafuse.constants import (
    BGEE_GENE_EXPRESSION_LEVELS_COL,
    DISGENET_DISEASE_COL,
    NAMESPACE_BINDINGS,
    NODE_TYPES,
    PREDICATES,
    URIS,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def replace_na_none(item):
    # Replace 'na' strings
    if isinstance(item, str) and item.lower() in ['na', 'nan', 'none']:
        return None
    # Replace Python None (and handle NaNs)
    elif item is None or (isinstance(item, float) and pd.isna(item)):
        return None
    # Recursively handle lists
    elif isinstance(item, list):
        return [replace_na_none(sub_item) for sub_item in item]
    # Recursively handle dictionaries
    elif isinstance(item, dict):
        return {key: replace_na_none(value) for key, value in item.items()}
    # Recursively handle NumPy arrays
    elif isinstance(item, np.ndarray):
        return np.array([replace_na_none(sub_item) for sub_item in item], dtype=object)
    # Return the item as-is if no replacement needed
    else:
        return item


def add_gene_node(g: Graph, gene_base_node: str, id_number: str, row) -> URIRef:
    """Create and add a gene node to the RDF graph.

    :param g: RDF graph to which the gene node will be added.
    :param gene_base_node: Base URI for the gene node.
    :param id_number: Unique identifier for the gene node.
    :return: URIRef for the created gene node.
    """
    target = row['target']
    if target:
        gene_node = URIRef(f"http://identifiers.org/ensembl/{target}")
        g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((gene_node, RDFS.label, Literal(row['identifier'], datatype= XSD.string)))
        return gene_node
    else:
        return None


def add_disease_node(g: Graph, disease_data: dict) -> URIRef:
    """Create and add a disease node to the RDF graph.

    :param g: RDF graph to which the disease node will be added.
    :param disease_data: Dictionary containing disease information.
    :return: URIRef for the created disease node.
    """
    # UMLS IRIs not in Bioregistry
    disease_curie = disease_data.get("disease_umlscui")
    if disease_curie is None:
        disease_curie = disease_data.get("UMLS")
    disease_iri = f"https://www.ncbi.nlm.nih.gov/medgen/{disease_curie}"
    disease_node = URIRef(disease_iri)
    g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
    #TODO discuss mappings between AE and disease? g.add((disease_node, RDF.type, URIRef(NODE_TYPES["adverse_event_node"])))
    g.add(
        (disease_node, RDFS.label, Literal(disease_data.get("disease_name"), datatype=XSD.string))
    )
    for identifier_type in [
        "HPO",
        "NCI",
        "OMIM",
        "MONDO",
        "ORDO",
        "EFO",
        "DO",
        "MESH",
        "UMLS",
    ]:  # TODO constants
        if disease_data.get(identifier_type):
            curie_field = disease_data.get(identifier_type)
            if "," in curie_field:
                curies = [i for i in curie_field.split(", ")]
            else:
                curies = [curie_field]
            for curie in curies:
                disease_source_iri = get_iri(curie)
                if disease_source_iri == None:
                    if ":" in curie:
                        curie = curie.split(":")[1]
                    disease_source_iri = get_iri("obo:" + curie)
                    g.add((disease_node, SKOS.closeMatch, URIRef(disease_source_iri))) # Some of the data does not look like a skos:exactMatch
                else:
                    g.add((disease_node, SKOS.closeMatch, URIRef(disease_source_iri)))
        else:
            pass
    return disease_node


def add_gene_disease_association_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, new_uris: dict
) -> URIRef:
    """Create and add a gene-disease association node to the RDF graph.

    :param g: RDF graph to which the association node will be added.
    :param id_number: Unique identifier for the association node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the gene.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :return: URIRef for the created association node.
    """
    gene_disease_assoc_node = URIRef(
        f"{new_uris['gene_disease_association']}/{id_number}/{source_idx}_assoc_{disease_id}"
    )
    g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene_disease_association"])))
    return gene_disease_assoc_node


def add_score_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, score: float, new_uris: dict, i
) -> URIRef:
    """Create and add a score node to the RDF graph.

    :param g: RDF graph to which the score node will be added.
    :param id_number: Unique identifier for the score node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the score.
    :param score: Score value for the gene-disease association.
    :return: URIRef for the created score node.
    """
    score_node = URIRef(f"{new_uris['score_base_node']}/{id_number}/{i}/{source_idx}_{disease_id}")
    g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
    g.add(
        (
            score_node,
            URIRef(NAMESPACE_BINDINGS["sio"] + "has_value"),
            Literal(score, datatype=XSD.double),
        )
    )
    return score_node


def add_data_source_node(g: Graph) -> URIRef: #TODO use this function to create all sources at once (pubchem assays, pathway sites, etc)
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :return: URIRef for the created data source node.
    """
    #for source in sources (eg disgenet)
    data_source_name = Literal("DisGeNET", datatype=XSD.string)
    data_source_url = URIRef("https://disgenet.com/")
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url


def add_gene_disease_associations(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    disease_data: list,
    new_uris: dict,
    i: int,
) -> None:
    """Process and add gene-disease association data to the RDF graph.

    :param g: RDF graph to which the associations will be added.
    :param id_number: Unique identifier for the associations.
    :param source_idx: Source index for the associations.
    :param gene_node: URIRef of the gene node associated with the disease data.
    :param disease_data: List of dictionaries containing disease association information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    """

    data = disease_data
    disease_node = add_disease_node(g, data)
    if disease_node:
        gene_disease_assoc_node = URIRef(
            f"{new_uris['gene_disease_association']}/{id_number}/{i}/{source_idx}"
        )
        g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene_disease_association"])))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), gene_node))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), disease_node))
        if data.get("score"):
            score_node = add_score_node(
                g=g,
                id_number=id_number,
                source_idx=source_idx,
                score=data["score"],
                new_uris=new_uris,
                i=i,
                disease_id=disease_data.get("disease_umlscui"),
            )  # SEE
            g.add(
                (
                    gene_disease_assoc_node,
                    URIRef(PREDICATES["sio_has_measurement_value"]),
                    score_node,
                )
            )
        # evidence_source_node = add_evidence_source_node(g, id_number, source_idx, dat["disease_umlscui"], data["evidence_source"])
        # g.add((gene_disease_assoc_node, URIRef(PREDICATES['sio_has_source']), evidence_source_node))
        # TODO EVIDENCE?
        data_source_node = add_data_source_node(g)
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_has_source"]), data_source_node))


def add_gene_expression_data(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    expression_data: list,
    experimental_process_data: list,
    new_uris: dict,
) -> None:
    """Process and add gene expression data to the RDF graph.

    :param g: RDF graph to which the expression data will be added.
    :param id_number: Unique identifier for the expression data.
    :param source_idx: Source index for the expression data.
    :param gene_node: URIRef of the gene node associated with the expression data.
    :param expression_data: List of dictionaries containing gene expression information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    """
    for data in expression_data:
        if pd.isna(data["anatomical_entity_id"]):
            continue

        gene_expression_value_node = URIRef(
            f"{new_uris['gene_expression_value_base_node']}/{id_number}/{source_idx}"
        )
        developmental_stage_node = URIRef(get_iri(data["developmental_stage_id"].replace("_", ":")))
        anatomical_entity_node = URIRef(get_iri(data["anatomical_entity_id"].replace("_", ":")))
        gene_expression_value_node = URIRef(
            f"{new_uris['gene_expression_value_base_node']}/{id_number}/{source_idx}_{data['anatomical_entity_id']}"
        )

        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), gene_expression_value_node))
        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), anatomical_entity_node))
        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), developmental_stage_node))

        # g.add((gene_node, URIRef(PREDICATES['sio_is_part_of']), cellular_component))
        # g.add((gene_node, URIRef(PREDICATES['sio_is_represented_by']), gene_symbol))

        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_has_value"]),
                Literal(data["expression_level"], datatype=XSD.double),
            )
        )

        g.add((anatomical_entity_node, RDF.type, URIRef(NODE_TYPES["anatomical_entity_node"])))
        g.add(
            (
                anatomical_entity_node,
                RDFS.label,
                Literal(data["anatomical_entity_name"], datatype=XSD.string),
            )
        )
        g.add(
            (
                developmental_stage_node,
                RDFS.label,
                Literal(data["developmental_stage_name"], datatype=XSD.string),
            )
        )
        g.add((anatomical_entity_node, SKOS.exactMatch, anatomical_entity_node))
        # g.add((URIRef(new_uris["life_cycle_base_node"]), SKOS.exactMatch, developmental_stage_node))

        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add((gene_expression_value_node, URIRef(PREDICATES["sio_has_input"]), gene_node))
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_has_output"]),
                gene_expression_value_node,
            )
        )

        # Add experimental node
        for exp in experimental_process_data:
            experimental_process_node = add_experimental_process_node(g=g, id_number=id_number, source_idx=source_idx, data=exp, new_uris=new_uris,)
            if experimental_process_node:
                g.add((gene_node, URIRef(PREDICATES['sio_is_part_of']), experimental_process_node))
                g.add((experimental_process_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
        ## Input gene, anatomical entity
        if experimental_process_node:
            g.add((experimental_process_node, URIRef(PREDICATES['sio_has_input']), gene_node))
            g.add((experimental_process_node, URIRef(PREDICATES['sio_has_input']), anatomical_entity_node))

def add_tested_substance_node(
    g: Graph,
    inchi: str,
    smiles:str,
    compound_id: str,
    compound_name: str,):
    uri_cas = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + str(compound_id).strip('CID')
    tested_substance_node= URIRef((uri_cas))
    g.add((tested_substance_node, RDF.type, URIRef(NODE_TYPES['tested_substance_node'])))
    g.add((tested_substance_node, RDFS.label, Literal(compound_name, datatype= XSD.string)))
    g.add((tested_substance_node, URIRef(PREDICATES['chebi_inchi']), Literal(inchi, datatype=XSD.string)))
    g.add((tested_substance_node, URIRef(PREDICATES['chebi_smiles']), Literal(smiles, datatype=XSD.string)))
    g.add((tested_substance_node, URIRef(PREDICATES['cheminf_compound_id']), Literal(compound_id, datatype=XSD.string)))
    return tested_substance_node

def add_experimental_process_node(
    g: Graph,
    id_number: str,
    source_idx: str,
    data: list,
    new_uris: dict,):
    pubchem_node = URIRef('https://pubchem.ncbi.nlm.nih.gov/')
    g.add((pubchem_node, RDF.type, URIRef(NODE_TYPES['experimental_process_node'])))
    pubchem_assay_id = data['pubchem_assay_id']
    if pubchem_assay_id is not None:
        pubchem_assay_iri = 'https://pubchem.ncbi.nlm.nih.gov/bioassay/' + str(pubchem_assay_id).strip('AID')
        assay_type =  data['assay_type']
        inchi =  data['InChI']
        outcome =  data['outcome']
        compound_cid =  data['compound_cid']
        compound_name =  data['compound_name']
        smiles =  data['SMILES']
        experimental_process_node = URIRef((
            f"{new_uris['experimental_process_node']}/{id_number}/{source_idx}_{pubchem_assay_id}"))
        g.add((experimental_process_node, RDF.type, URIRef(NODE_TYPES["experimental_process_node"])))
        g.add((experimental_process_node, RDF.type, URIRef(NODE_TYPES["experimental_process_node"])))
        g.add((experimental_process_node, RDFS.label, Literal(assay_type, datatype= XSD.string)))
        g.add((experimental_process_node, DC.identifier, Literal(pubchem_assay_id, datatype=XSD.string)))
        # has tested substance
        tested_substance_node = add_tested_substance_node(g=g, inchi=inchi, smiles=smiles, compound_id=compound_cid, compound_name = compound_name)
        g.add((experimental_process_node, URIRef(PREDICATES['sio_has_input']), tested_substance_node))
        # has outcome
        g.add((experimental_process_node, URIRef(PREDICATES['sio_has_output']), Literal(outcome))) #TODO discuss whether to map to class, and which are the possible values?
        return experimental_process_node

def add_pathway_node(g: Graph,
    data: list,
    source : str):
    pathway_label = data['pathway_label']
    pathway_id = data['pathway_id']
    if pathway_id is not None:
        if source == 'WikiPathways':
            pathway_iri = 'https://www.wikipathways.org/pathways/' + str(pathway_id)
            iri_source = 'https://www.wikipathways.org/'

        if source == 'OpenTargets_reactome':
            pathway_iri = 'https://reactome.org/content/detail/' + str(pathway_id)
            iri_source = 'https://reactome.org/'
            source = 'Reactome'
        if source == 'MINERVA':
            pathway_iri = "https://minerva-net.lcsb.uni.lu/api/"  + str(pathway_id) #TODO just a placeholder
            iri_source = "https://minerva-net.lcsb.uni.lu/"
        pathway_node = URIRef(pathway_iri)
        g.add((pathway_node, RDF.type, URIRef(NODE_TYPES['pathway_node'])))
        g.add((pathway_node, RDFS.label, Literal(pathway_label, datatype=XSD.string)))
        g.add((pathway_node, URIRef(PREDICATES['sio_has_source']), URIRef(iri_source)))
        g.add((URIRef(iri_source), RDFS.label, Literal(source, datatype= XSD.string))) # TODO use function datasources  abvove to avoid redundance adding data sources
        gene_count = data.get('pathway_id', None)
        # TODO discuss gene count
        #if gene_count:
        #    g.add((pathway_node, ))
        return pathway_node
    else:
        return None

def add_compound_node(g, compound, gene_node):
    chembl_id = compound['chembl_id']
    drugbank_id = compound['drugbank_id']
    compound_name = compound['compound_name']
    compound_cid = compound['compound_cid'] #TODO report issue with CID not having the same format as in tested_substance_data
    is_approved = compound['is_approved']
    clincal_trial_phase = compound['clincal_trial_phase'] #TODO report issue typo
    clinical_phases = {
        "1.0": "http://purl.obolibrary.org/obo/OPMI_0000368",
        "2.0": "http://purl.obolibrary.org/obo/OPMI_0000369",
        '3.0': 'http://purl.obolibrary.org/obo/OPMI_0000370',
        '4.0': 'http://purl.obolibrary.org/obo/OPMI_0000371'
    }  # TODO add to constants

    relation = compound['relation'] #TODO right now the information is lost, as all classes (agonist, etc.) are turned into activates / inhibit. This complicates mappings. These two mappings are not technically 100% correct because we miss information on the MOA --is it direct or indirect (ant)agonism?
    relations = {"activates": "http://purl.obolibrary.org/obo/RO_0018027", #agonist
                 "inhibits": "http://purl.obolibrary.org/obo/RO_0018029" #antagonist
                 } #TODO add to constants (?)
    # adverse_effect_count = compound['adverse_effect_count']

    compound_node = URIRef(f'https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}')
    g.add((compound_node, RDFS.label, Literal(str(compound_name), datatype= XSD.string)))
    g.add((compound_node, RDF.type, URIRef(NODE_TYPES['tested_substance_node'])))
    if relation in relations.keys():
        relation_iri = relations[relation]
        g.add((compound_node, URIRef(relation_iri), gene_node))
    if clincal_trial_phase:
        clinical_phase_iri = clinical_phases[str(clincal_trial_phase)]
        g.add(
            (
                compound_node,
                URIRef("http://purl.obolibrary.org/obo/PATO_0000083"),
                URIRef(clinical_phase_iri),
            )
        )  # TODO add to constants
    for id in [drugbank_id, compound_cid]:
        # TODO add namespace mappings to constant to replace current if statements
        if id:
            iri = None
            if id == drugbank_id:
                iri = f'https://go.drugbank.com/drugs/{id}'
            if id == compound_cid:
                iri = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + str(compound_cid)
            if iri:
                id_node = URIRef(iri)  #TODO decide whether it is necesdsary for exactmatches to also contain annotations and properties
                g.add((compound_node, SKOS.exactMatch, id_node))
                g.add((id_node, SKOS.exactMatch, compound_node))
                g.add((id_node, RDFS.label, Literal(compound_name, datatype= XSD.string)))
                g.add((id_node, RDF.type, URIRef(NODE_TYPES['tested_substance_node'])))
            if type(compound['adverse_effect']) == list:
                for i in compound['adverse_effect']:
                    ae = i['name']
                    ae_node = add_ae_node(g, ae) # TODO file issue: For now we make up an internal IRI for each AE, maybe we could retrieve IRIs for AEs from OpenTargets?
                    if ae_node:
                        g.add(((compound_node, URIRef(PREDICATES['has_adverse_event']), ae_node)))                               
    return compound_node

def add_ae_node(g, ae):
    if ae:
        ae_f = str(ae).replace(" ", "_")
        ae_node = URIRef(f'https://biodatafuse.org/rdf/ae/{ae_f}') #TODO get IRIs from OpenTarget
        g.add((ae_node, RDFS.label, Literal(ae, datatype= XSD.string)))
        g.add((ae_node, RDF.type, URIRef(NODE_TYPES['adverse_event_node'])))
        return ae_node
    else:
        return None


# TODO EI, EL
# TODO CHECK REST OF THE RDF SCHEMA
# TODO REFACTOR, CLEAN
# TODO make sure all constants are in constants.py
# TODO


def generate_rdf(df: pd.DataFrame, BASEURI: str) -> Graph:
    """Generate RDF graph from the property table.

    :param df: DataFrame containing the data.
    :param BASEURI: String containing the base URI for the nodes in the graph.
    :return: RDF graph.
    """
    g = Graph()
    g.bind("dc", DC)
    new_uris = URIS
    for key, value in new_uris.items():
        new_value = BASEURI + value
        new_uris[key] = new_value
        # logger.info(f'Binding {key} to {new_value}')
        g.bind(key, new_value)

    # Get namespace bindings
    for key, value in NAMESPACE_BINDINGS.items():
        g.bind(key, value)
    df = df.applymap(replace_na_none)
    for i, row in df.iterrows():
        # Unpack the relevant columns
        source_idx = row["identifier"]
        source_namespace = row["identifier.source"]
        target_idx = row["target"]
        target_namespace = row["target.source"]
        expression_data = row[BGEE_GENE_EXPRESSION_LEVELS_COL]
        experimental_process_data = row['PubChem_assays']
        disease_data = []
        for source in [DISGENET_DISEASE_COL, "OpenTargets_diseases"]:  # TODO update constants
            source_el = row[source]
            if isinstance(source_el, list):
                disease_data += source_el

        # Check if any of the essential columns are NaN, and skip the iteration if so
        if (
            pd.isna(source_idx)
            or pd.isna(source_namespace)
            or pd.isna(target_idx)
            or pd.isna(target_namespace)
        ):
            # logger.warning(f"Skipping row {i} due to missing essential data.")
            continue

        # Normalize CURIEs and convert to IRIs
        source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
        target_curie = normalize_curie(f"{target_namespace}:{target_idx}")

        # Skip the row if CURIE normalization fails (e.g., if it's None)
        if not source_curie or not target_curie:
            # logger.warning(f"Skipping row {i} due to CURIE normalization failure.")
            continue

        id_number = f"{i:06d}"  # Create ID for this row
        gene_node = add_gene_node(g, new_uris["gene_base_node"], id_number, row)

        # Add gene-disease associations
        if len(disease_data) > 0:
            i = 0
            for disease in disease_data:
                add_gene_disease_associations(
                    g, id_number, source_idx, gene_node, disease, new_uris, i
                )
                i += 1

        # Add gene expression data
        add_gene_expression_data(g, id_number, source_idx, gene_node, expression_data, experimental_process_data, new_uris)

        # Add pathway nodes
        for source in ['WikiPathways', 'MINERVA', 'OpenTargets_reactome']:
            for pathway_data in row[source]:
                if pathway_data['pathway_id']:
                    pathway_node = add_pathway_node(g=g, data= pathway_data, source=source)
                    g.add((gene_node, URIRef(PREDICATES['sio_is_part_of']), pathway_node))
                    g.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))

        ## Source
        ## Add process data
        processes_data = _data = row['OpenTargets_go']
        for process_data in processes_data:
            go_cpf = add_go_cpf(g, process_data)
            if go_cpf:
                g.add((gene_node, URIRef(PREDICATES['sio_is_part_of']), go_cpf))
                g.add((go_cpf, URIRef(PREDICATES["sio_has_part"]), gene_node))
        ## Add compound data #TODO use constants
        compound_data = row['OpenTargets_compounds']
        for compound in compound_data:
            add_compound_node(g, compound, gene_node)

        # is compound active or inactive

        # Adverse effects
        # for adverse effect in ['OpenTargets_compounds'][i][j]['adverse_effect'] create ae node

    return g

def add_go_cpf(g, process_data) -> URIRef:
    go_types = {
        'C' : 'http://purl.obolibrary.org/obo/GO_0005575',
        'P' : "http://purl.obolibrary.org/obo/GO_0008150",
        'F' : "http://purl.obolibrary.org/obo/GO_0003674"
    } # TODO add to constants
    curie = process_data['go_id']
    if curie:
        iri = get_iri(curie)
        label = process_data['go_name']
        go_type = URIRef(go_types[process_data['go_type']]) if process_data['go_type'] in list(go_types.keys()) else None
        go_cpf = URIRef(iri)
        g.add((go_cpf, RDFS.label, Literal(label, datatype= XSD.string)))
        if go_type:
            g.add((go_cpf, RDF.type, go_type))
        return go_cpf
    else:
        return None
