from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
from KINEPIK_app.database_con import session_local
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# naming the blueprint
protein_bp = Blueprint("protein_bp", __name__, url_prefix="/api/proteins")

@protein_bp.route("/info", methods = ["GET"])
def instructions():
    '''This function hands out the user info about the possible parameters that can be given to the proteins queries. 
    The return returns the info from string/text to html which is then displayed to the user'''

    # information about the parameters
    protein_id_info = "protein_ids : uniprot ID, if multiple use comma as a separator"
    fields_info = "fields : kinase, mappedgene"
    # the format of the html page
    info = f"""
    <html>
        <body>
            <p>{protein_id_info}</p>
            <p>{fields_info}</p>
        </body>
    </html>
    """

    # return changes the format of the string to html that will be displayed to the user
    return Response(info, mimetype = "text/html")

@protein_bp.route("/results", methods=["GET"])
def get_protein():
    '''This function takes the url and based on the parameters given by the user, it will return information
    about the requested protein(s) in json format'''

    # creating a connection to the database
    session = session_local()
    # collecting the parameters given by the user in the url
    protein_ids = request.args.get("protein_ids")
    fields = request.args.get("fields")

    try:
        # if the user has not given any uniprot ids that are required specified error message will be shown
        if protein_ids == None:
            return f"Error: Uniprot ID required for the query"
        # if uniprot id(s) in the url this will move forward
        else:
            # collecting all the uniprot ids
            protein_ids = protein_ids.split(",")
            protein_infos = []
            # going through the list of uniprot ids given by the user
            for protein_id in protein_ids:
                protein_id = protein_id.upper()
                #print(protein_id)
                # collecting alll the info about the protein in the user's list from the Protein-table that is in the database
                all_protein_info = session.query(Protein).filter_by(id=protein_id).first()

                # if the protein is found from the database and user has given no parameters or all possible parameters, all possible info collected in dictionary format
                if all_protein_info and (fields == None or (fields and "kinase" in fields and "mappedgene" in fields)):
                    protein_info = [{
                            "UniprotID" : all_protein_info.id,
                            "UniprotName" : all_protein_info.name,
                            "KinaseInfo" : {
                                "IsKinase" : all_protein_info.kinase,
                                "KinaseGroup" : all_protein_info.kinase_group,
                                "KinaseFamily" : all_protein_info.kinase_family,
                                "KinaseSubfamily" : all_protein_info.kinase_subfamily
                            },
                            "Sequence" : all_protein_info.sequence,
                            "SequenceLength" : all_protein_info.length,
                            "Go-ids" : all_protein_info.go,
                            "Description" : all_protein_info.description,
                            "GeneInfo" : {
                                "MappedGene" : all_protein_info.gene,
                                "GeneSynonyms" : all_protein_info.gene_synonyms
                            },
                            "OrganismID" : all_protein_info.species
                        }]
                    # to have multiple proteins in one page, all the results are added to a list
                    protein_infos.append(protein_info)

                # if the user inputs fields=kinase and the protein is in the database, only basic info + kinase info shown
                elif all_protein_info and fields and "kinase" in fields:
                    protein_info = [{
                            "UniprotID" : all_protein_info.id,
                            "UniprotName" : all_protein_info.name,
                            "KinaseInfo" : {
                                "IsKinase" : all_protein_info.kinase,
                                "KinaseGroup" : all_protein_info.kinase_group,
                                "KinaseFamily" : all_protein_info.kinase_family,
                                "KinaseSubfamily" : all_protein_info.kinase_subfamily
                            },
                            "Sequence" : all_protein_info.sequence,
                            "SequenceLength" : all_protein_info.length,
                            "Go-ids" : all_protein_info.go,
                            "Description" : all_protein_info.description,
                            "OrganismID" : all_protein_info.species
                        }]
                    # to have multiple proteins in one page, all the results are added to a list
                    protein_infos.append(protein_info)

                # if the user inputs fields=mappedgene and the protein is in the database, only basic info + gene info shown
                elif all_protein_info and fields and "mappedgene" in fields:
                    protein_info = [{
                                "UniprotID" : all_protein_info.id,
                                "UniprotName" : all_protein_info.name,
                                "Sequence" : all_protein_info.sequence,
                                "SequenceLength" : all_protein_info.length,
                                "go-ids" : all_protein_info.go,
                                "Description" : all_protein_info.description,
                                "GeneInfo" : {
                                    "MappedGene" : all_protein_info.gene,
                                    "GeneSynonyms" : all_protein_info.gene_synonyms
                                },
                                "OrganismID" : all_protein_info.species
                            }]
                    # to have multiple proteins in one page, all the results are added to a list
                    protein_infos.append(protein_info)
            #print(protein_infos)

            # the list of the dictionaries is then transformed into a json format for the page
            return jsonify(protein_infos)

    # closing the connection to the database at the end
    finally:
        session.close()

