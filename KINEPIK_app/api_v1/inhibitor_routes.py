from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
from KINEPIK_app.database_con import session_local
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# naming the blueprint
inhibitor_bp = Blueprint("inhibitor_bp", __name__, url_prefix="/api/inhibitor")

@inhibitor_bp.route("/info", methods = ["GET"])
def instructions():
    '''This function hands out the user info about the possible parameters that can be given to the inhibitor queries. 
    The return returns the info from string/text to html which is then displayed to the user'''
    return

@inhibitor_bp.route("/all", methods = ["GET"])
def get_all_kinases():
    # connection to the database
    session = session_local()
    kinases_json = []
    phosphosites = request.args.get("phosphosites")
    if phosphosites is not None:
        phosphosites = int(phosphosites)

    try:
        if phosphosites == None or phosphosites == 0:
            all_kinase_info = session.query(Protein).filter_by(kinase=1).all()
            for kinase in all_kinase_info:
                kinase_info = [{
                        "UniprotID" : kinase.id,
                        "UniprotName" : kinase.name,
                        "KinaseInfo" : {
                            "IsKinase" : kinase.kinase,
                            "KinaseGroup" : kinase.kinase_group,
                            "KinaseFamily" : kinase.kinase_family,
                            "KinaseSubfamily" : kinase.kinase_subfamily
                            },
                        "Sequence" : kinase.sequence,
                        "SequenceLength" : kinase.length,
                        "Go-IDs" : kinase.go,
                        "Description" : kinase.description,
                        "GeneInfo" : {
                            "MappedGene" : kinase.gene,
                            "GeneSynonyms" : kinase.gene_synonyms
                            },
                        "OrganismID" : kinase.species                        }]
                kinases_json.append(kinase_info)
            return jsonify(kinases_json)
        
        elif phosphosites == 1:
            all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
            for kinase in all_kinase_info:
                kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = kinase.id).all()
                target_phosphosites_info = session.query(Interaction).filter_by(source = kinase.id).all()
                kinase_phos_list = []
                target_phos_list = []
                for k_phos in kinase_phosphosites_info:
                    kinase_phos_list.append(k_phos.phosphosite_id)
                for t_phos in target_phosphosites_info:
                    target_phos_list.append(t_phos.phosphosite_id)
                kinase_info = [{
                        "UniprotID" : kinase.id,
                        "UniprotName" : kinase.name,
                        "KinaseInfo" : {
                            "IsKinase" : kinase.kinase,
                            "KinaseGroup" : kinase.kinase_group,
                            "KinaseFamily" : kinase.kinase_family,
                            "KinaseSubfamily" : kinase.kinase_subfamily
                            },
                        "TargetPhosphosites": target_phos_list,
                        "PhosphositesOnKinase": kinase_phos_list,
                        "Sequence" : kinase.sequence,
                        "SequenceLength" : kinase.length,
                        "Go-IDs" : kinase.go,
                        "Description" : kinase.description,
                        "GeneInfo" : {
                            "MappedGene" : kinase.gene,
                            "GeneSynonyms" : kinase.gene_synonyms
                            },
                        "OrganismID" : kinase.species                        
                        }]
                kinases_json.append(kinase_info)

            return jsonify(kinases_json)

    finally:
        session.close()


@kinase_bp.route("/specific", methods = ["GET"])
def get_kinase():
    # connection to the database
    session = session_local()
    kinases_json = []
    kinase_ids = request.args.get("kinase_ids")
    kinase_ids = kinase_ids.split(",")
    phosphosites = request.args.get("phosphosites")

    try:
        if kinase_ids == None:
            return 
        else:
            if phosphosites == None:
                for kin in kinase_ids:
                    all_kinase_info = session.query(Protein).filter_by(id=kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
                    target_phosphosites_info = session.query(Interaction).filter_by(source = all_kinase_info.id).all()
                    kinase_phos_list = []
                    target_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    for t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = [{
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "TargetPhosphosites": target_phos_list,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "OrganismID" : all_kinase_info.species                        }]
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif phosphosites == "targets":
                for kin in kinase_ids:
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    target_phosphosites_info = session.query(Interaction).filter_by(source = all_kinase_info.id).all()
                    target_phos_list = []
                    for t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = [{
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "TargetPhosphosites": target_phos_list,
                            "OrganismID" : all_kinase_info.species                        }]
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif phosphosites == "sites":
                for kin in kinase_ids:
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
                    kinase_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    kinase_info = [{
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "OrganismID" : all_kinase_info.species                        }]
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
        
    finally:
        session.close()