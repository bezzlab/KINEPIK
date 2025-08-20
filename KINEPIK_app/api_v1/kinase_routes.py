from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
#from KINEPIK_app.database_con import session_local
from KINEPIK_app.database_con import getSessionForVersion
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# naming the blueprint
kinase_bp = Blueprint("kinase_bp", __name__, url_prefix="/api/<int:version>/kinases")

@kinase_bp.route("/info/<route>", methods = ["GET"])
def all_instructions(version,route):
    '''This function hands out the user info about the possible parameters that can be given to the kinases queries. 
    The return returns the info from string/text to html which is then displayed to the user'''
    # information about the parameters
    phosphosites_info_all = "phosphosites : 0 or 1"
    kinase_id_info = "kinase_ids : Uniprot ID, if multiple use comma (,) as a separator"
    phosphosites_info = "phosphosites : sites,targets"
    confidence_info = "confidence : 0 or 1"
    # the format of the all html page
    if route == "all":
        info = f"""
        <html>
            <body>
                <p>{phosphosites_info_all}</p>
            </body>
        </html>
        """
    # the format of thespecific html page
    if route == "specific":
        info = f"""
        <html>
            <body>
                <p>{kinase_id_info}</p>
                <p>{phosphosites_info}</p>
                <p>{confidence_info}</p>
            </body>
        </html>
        """

    # return changes the format of the string to html that will be displayed to the user
    return Response(info, mimetype = "text/html")

@kinase_bp.route("/all", methods = ["GET"])
def get_all_kinases(version):
    '''Function returns general information about all the kinases in the database. Phosphosite filter can be applied and then also kinase's phosphosites and its targets will be listed.
    Returns a json'''
    # connection to the database
    Session = getSessionForVersion(version)
    session = Session()
    kinases_json = []
    # user input
    phosphosites = request.args.get("phosphosites")
    if phosphosites is not None:
        phosphosites = int(phosphosites)

    try:
        # if user does not use phosphosite filter or inputs 0 there
        if phosphosites == None or phosphosites == 0:
            # collecting kinase info
            all_kinase_info = session.query(Protein).filter_by(kinase=1).all()
            # going through the protein info kinase by kinase and adding it to the json
            for kinase in all_kinase_info:
                kinase_info = {
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
                        "OrganismID" : kinase.species                        
                        }
                kinases_json.append(kinase_info)
            # returning the full json to user
            return jsonify(kinases_json)
        
        # if user wants phosphosites
        elif phosphosites == 1:
            # collecting kinase info
            all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
            # going through the protein info kinase by kinase and adding it to the json
            for kinase in all_kinase_info:
                # collecting phosphosite information from Phosphosite and Interaction tables based on the kinase id
                kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = kinase.id).all()
                target_phosphosites_info = session.query(Interaction).filter_by(source = kinase.id).all()
                # making lists out of the own and target phosphosites
                kinase_phos_list = []
                target_phos_list = []
                for k_phos in kinase_phosphosites_info:
                    kinase_phos_list.append(k_phos.phosphosite_id)
                for t_phos in target_phosphosites_info:
                    target_phos_list.append(t_phos.phosphosite_id)
                # adding all the info in to json
                kinase_info = {
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
                        }
                kinases_json.append(kinase_info)

            #returning the json to the user
            return jsonify(kinases_json)

    finally:
        session.close()


@kinase_bp.route("/specific", methods = ["GET"])
def get_kinase(version):
    '''Function returns general information about the kinases that user asks for. Phosphosite filter can be applied and then user can specify what type of phosphosites they want.
    Confidence filter will only show high confidence phosphosites. Returns a json'''
    # connection to the database
    Session = getSessionForVersion(version)
    session = Session()
    kinases_json = []
    # userr input
    kinase_ids = request.args.get("kinase_ids")
    kinase_ids = kinase_ids.split(",")
    phosphosites = request.args.get("phosphosites")
    phos_confidence = request.args.get("confidence")
    if phos_confidence is not None:
        phos_confidence = int(phos_confidence)

    try:
        # error handling if user did not provide a kinase id
        if kinase_ids == None:
            return "Kinase id required for results"
        else:
            # user did not input phosphosites or confidence
            if phosphosites == None and (phos_confidence == 0 or phos_confidence == None):
                # going through the list of kinases given by the user
                for kin in kinase_ids:
                    # making sure that the id format matches the database
                    kin = kin.upper()
                    # collecting data from the database
                    all_kinase_info = session.query(Protein).filter_by(id=kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
                    target_phosphosites_info = session.query(Interaction).filter(Interaction.source == all_kinase_info.id).filter(Interaction.target != all_kinase_info.id).all()
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

            elif phosphosites == None and phos_confidence == 1:
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id=kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).filter_by(confidence = "High").all()
                    target_phosphosites_info = session.query(Interaction,Phosphosite).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source == all_kinase_info.id).filter(Phosphosite.confidence == "High").all()
                    #target_phosphosites_info = session.query(Interaction).filter_by(source = all_kinase_info.id).filter_by(confidence = "high").all()
                    kinase_phos_list = []
                    target_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    for inter_id,t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "TargetPhosphosites": target_phos_list,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif phosphosites == "targets" and (phos_confidence == 0 or phos_confidence == None):
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    target_phosphosites_info = session.query(Interaction).filter(Interaction.source == all_kinase_info.id).filter(Interaction.target != all_kinase_info.id).all()
                    target_phos_list = []
                    for t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "TargetPhosphosites": target_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif phosphosites == "targets" and phos_confidence == 1:
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    target_phosphosites_info = session.query(Interaction,Phosphosite).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source == all_kinase_info.id).filter(Phosphosite.confidence == "High").all()
                    target_phos_list = []
                    for inter_id,t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "TargetPhosphosites": target_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            
            elif phosphosites == "sites" and (phos_confidence == 0 or phos_confidence == None):
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
                    kinase_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif phosphosites == "sites" and phos_confidence == 1:
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).filter_by(confidence = "High").all()
                    kinase_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif ("sites" in phosphosites and "targets" in phosphosites) and (phos_confidence == 0 or phos_confidence == None):
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id = kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
                    target_phosphosites_info = session.query(Interaction).filter(Interaction.source == all_kinase_info.id).filter(Interaction.target != all_kinase_info.id).all()
                    kinase_phos_list = []
                    target_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    for t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "TargetPhosphosites": target_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
            
            elif ("sites" in phosphosites and "targets" in phosphosites) and phos_confidence == 1:
                for kin in kinase_ids:
                    kin = kin.upper()
                    all_kinase_info = session.query(Protein).filter_by(id=kin).first()
                    kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).filter_by(confidence = "High").all()
                    target_phosphosites_info = session.query(Interaction,Phosphosite).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source == all_kinase_info.id).filter(Phosphosite.confidence == "High").all()
                    #target_phosphosites_info = session.query(Interaction).filter_by(source = all_kinase_info.id).filter_by(confidence = "high").all()
                    kinase_phos_list = []
                    target_phos_list = []
                    for k_phos in kinase_phosphosites_info:
                        kinase_phos_list.append(k_phos.phosphosite_id)
                    for inter_id,t_phos in target_phosphosites_info:
                        target_phos_list.append(t_phos.phosphosite_id)
                    kinase_info = {
                            "SourceUniprotID" : all_kinase_info.id,
                            "UniprotName" : all_kinase_info.name,
                            "TargetPhosphosites": target_phos_list,
                            "PhosphositesOnKinase": kinase_phos_list,
                            "OrganismID" : all_kinase_info.species                        
                            }
                    kinases_json.append(kinase_info)
                return jsonify(kinases_json)
        
    finally:
        session.close()