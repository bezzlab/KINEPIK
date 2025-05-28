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
    pert_json = []
    pert_type = request.args.get("type")

    try:
        if pert_type == None or (pert_type and "small_molecule" in pert_type and "knockout" in pert_type):
            all_inhibitor_info = session.query(Perturbation).all()
            for inhibitor in all_inhibitor_info:
                if inhibitor.type == "small molecule":
                    inhibitor_info = [{
                            "PerturbationName" : inhibitor.name,
                            "Type" : inhibitor.type,
                            "PubChemCID" :  inhibitor.pubchem,
                            "SMILES" : inhibitor.smiles,
                            "Action" : inhibitor.action,
                            "Synonyms" : inhibitor.synonyms,
                            }]
                    pert_json.append(inhibitor_info)
                elif inhibitor.type == "CRISPR knockout":
                    inhibitor_info = [{
                            "PerturbationName" : inhibitor.name,
                            "Gene" : inhibitor.gene,
                            "Type" : inhibitor.type,
                            "Action" : inhibitor.action
                            }]
                    pert_json.append(inhibitor_info)
            return jsonify(pert_json)
        
        elif pert_type and "small_molecule" in pert_type:
            all_inhibitor_info = session.query(Perturbation).filter(Perturbation.type =="small molecule").all()
            for inhibitor in all_inhibitor_info:
                inhibitor_info = [{
                                "PerturbationName" : inhibitor.name,
                                "Type" : inhibitor.type,
                                "PubChemCID" :  inhibitor.pubchem,
                                "SMILES" : inhibitor.smiles,
                                "Action" : inhibitor.action,
                                "Synonyms" : inhibitor.synonyms,
                                }]
                pert_json.append(inhibitor_info)
            return jsonify(pert_json)
        
        elif pert_type and "knockout" in pert_type:
            all_inhibitor_info = session.query(Perturbation).filter(Perturbation.type =="CRISPR knockout").all()
            for inhibitor in all_inhibitor_info:
                inhibitor_info = [{
                                "PerturbationName" : inhibitor.name,
                                "Gene" : inhibitor.gene,
                                "Type" : inhibitor.type,
                                "Action" : inhibitor.action
                                }]
                pert_json.append(inhibitor_info)
            return jsonify(pert_json)

    finally:
        session.close()


@inhibitor_bp.route("/specific", methods = ["GET"])
def get_inhibitor():
    # connection to the database
    session = session_local()
    inhibitor_json = []
    inhibitors = request.args.get("inhibitors")
    if inhibitors is not None:
        inhibitors = inhibitors.split(",")

    try:
        if inhibitors == None:
            return "Inhibitor name required for results"
        else:
            for inhibitor in inhibitors:
                all_inhibitor_info = session.query(Perturbation_interaction).filter_by(perturbation = inhibitor).all()
                kinase_list = []
                for inhibit in all_inhibitor_info:
                    if inhibit.target is not None:
                        kinase_list.append(inhibit.target)
                pert_info = [{
                    "PerturbationName" : inhibitor,
                    "TargetKinases" : kinase_list
                    }]
                inhibitor_json.append(pert_info)
            
            return jsonify(inhibitor_json)
        
    finally:
        session.close()