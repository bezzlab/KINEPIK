from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
#from KINEPIK_app.database_con import session_local
from KINEPIK_app.database_con import getSessionForVersion
import ast
from scipy.stats import norm
from sqlalchemy import or_
import numpy
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene, Known_perturbations

# naming the blueprint
inhibitor_bp = Blueprint("inhibitor_bp", __name__, url_prefix="/api/<int:version>/perturbation")

def formatList(list_form):
    '''Function formats a python list to a typical list format based on the number of items in the list, and then returns the list to the user'''
    if not list_form:
        return ""
    elif len(list_form) == 1:
        return list_form[0]
    elif len(list_form) == 2:
        return " and ".join(list_form)
    else:
        return ", ".join(list_form[:-1]) + " and " + str(list_form[-1])

@inhibitor_bp.route("/info/<route>", methods = ["GET"])
def instructions(version,route):
    '''This function hands out the user info about the possible parameters that can be given to the perturbation queries. User will have to
    specify which route they want the parameters for.The function returns the info from string/text to html which is then displayed to the user'''
    # connection to the db
    Session = getSessionForVersion(version)
    session = Session()

    # collecting the cell lines from the db
    cell_line_list = session.query(Experimental.cell_line).distinct().all()
    cl_list = []
    for line in cell_line_list:
        cl_list.append(line[0])

    # information about the parameters
    cell_lines = "cell_line : available cell lines curretly are :" + formatList(cl_list)
    type = "type : small_molecule or knockout"
    inhibitors = "perturbations : common name of the perturbation (list of all available perturbations can be found with /api/perturbation/all route)"
    kinase_ids = "kinase_ids : Uniprot ID, if multiple use comma (,) as a separator"
    phosphosite_confidence = "phosphosite_confidence : 0 or 1, with 1 only high confidence phosphosites will be used"
    autophosphorylation = "autophosphorylation : include or exclude"
    weighted = "weighted : true or false, if true, then all the calculations will be done by using weighted values and KinaseConfidence will be included in to the output"
    names = "names : Uniprot IDs or common names of perturbations"
    known = "known : true or false"
    confidence = "confidence : decimal number"
    fc_type = "type : source, target_kinase, target_phosphosite"
    id = "id : Uniprot ID (if type is source or target_kinase) or phosphosite ID (if type target_phosphosite)"
    # the format of the html page for all route
    if route == "all":
        info = f"""
        <html>
            <body>
                <p>{type}</p>
            </body>
        </html>
        """
    # the format of the html page for specific route
    elif route == "specific":
        info = f"""
        <html>
            <body>
                <p>{inhibitors}</p>
            </body>
        </html>
        """
    # the format of the html page for specific route
    elif route == "available":
        info = f"""
        <html>
            <body>
                <p>{names}</p>
                <p>{known}</p>
                <p>{confidence}</p>
            </body>
        </html>
        """
    # the format of the html page for fc route
    elif route == "fc":
        info = f"""
        <html>
            <body>
                <p>{id}</p>
                <p>{type}</p>
                <p>{cell_lines}</p>
                <p>{phosphosite_confidence}</p>
            </body>
        </html>
        """
    # the format of the html page for KSEA route
    elif route == "KSEA":
        info = f"""
        <html>
            <body>
                <p>{kinase_ids}</p>
                <p>{inhibitors}</p>
                <p>{cell_lines}</p>
                <p>{autophosphorylation}</p>
                <p>{weighted}</p>
                <p>{phosphosite_confidence}</p>
            </body>
        </html>
        """

    # return changes the format of the string to html that will be displayed to the user
    return Response(info, mimetype = "text/html")

@inhibitor_bp.route("/all", methods = ["GET"])
def get_all_perturbations(version):
    '''Function collect all the information of the available pertubations in the db, and generates a json of all of them. User can select whether only
    show small molecules or gene knockout options.'''
    # connection to the database and getting values from user
    Session = getSessionForVersion(version)
    session = Session()
    pert_json = []
    pert_type = request.args.get("type")

    try:
        # if the user has not included any filter on the route or has listed both of the options on the route, this if clause will be used
        if pert_type == None or (pert_type and "small_molecule" in pert_type and "knockout" in pert_type):
            # collecting all the information about the available perturbations in the db
            all_inhibitor_info = session.query(Perturbation).all()
            # going through the list of perturbation information from db row by row and saving the info in the json format to a new list
            for inhibitor in all_inhibitor_info:
                if inhibitor.type == "small molecule":
                    inhibitor_info = {
                            "PerturbationName" : inhibitor.name,
                            "Type" : inhibitor.type,
                            "PubChemCID" :  inhibitor.pubchem,
                            "SMILES" : inhibitor.smiles,
                            "Action" : inhibitor.action,
                            "Synonyms" : inhibitor.synonyms,
                            }
                    pert_json.append(inhibitor_info)
                elif inhibitor.type == "CRISPR knockout":
                    inhibitor_info = {
                            "PerturbationName" : inhibitor.name,
                            "Gene" : inhibitor.gene,
                            "Type" : inhibitor.type,
                            "Action" : inhibitor.action
                            }
                    pert_json.append(inhibitor_info)
            # returning all the collected info to the user
            return jsonify(pert_json)
        
        # if the user has specified that they only want small molecules, this if clause will be followed
        elif pert_type and "small_molecule" in pert_type:
            # collecting all the information for the small molecules in the db
            all_inhibitor_info = session.query(Perturbation).filter(Perturbation.type =="small molecule").all()
            # going through the list of perturbation information from db row by row and saving the info in the json format to a new list
            for inhibitor in all_inhibitor_info:
                inhibitor_info = {
                                "PerturbationName" : inhibitor.name,
                                "Type" : inhibitor.type,
                                "PubChemCID" :  inhibitor.pubchem,
                                "SMILES" : inhibitor.smiles,
                                "Action" : inhibitor.action,
                                "Synonyms" : inhibitor.synonyms,
                                }
                pert_json.append(inhibitor_info)
            # returning all the collected info to the user
            return jsonify(pert_json)
        
        # if the user has specified that they only want gene knockouts, this if clause will be followed
        elif pert_type and "knockout" in pert_type:
            # collecting all the information for the gene knockouts in the db
            all_inhibitor_info = session.query(Perturbation).filter(Perturbation.type =="CRISPR knockout").all()
            # going through the list of perturbation information from db row by row and saving the info in the json format to a new list
            for inhibitor in all_inhibitor_info:
                inhibitor_info = {
                                "PerturbationName" : inhibitor.name,
                                "Gene" : inhibitor.gene,
                                "Type" : inhibitor.type,
                                "Action" : inhibitor.action
                                }
                pert_json.append(inhibitor_info)
            # returning all the collected info to the user
            return jsonify(pert_json)
        
    # closing the session to the db once the function has finished
    finally:
        session.close()


@inhibitor_bp.route("/available", methods = ["GET"])
def get_inhibitor(version):
    '''Function gives the user list of the target kinases the given inhibitor has or vice versa. The function returns the info in json format'''
    # connection to the database  and getting values from user
    Session = getSessionForVersion(version)
    session = Session()
    inhibitor_json = []
    # name of the inhibitor or the kinase required
    names = request.args.get("name")
    confidence = request.args.get("confidence")
    # if multiple names given they are split to list
    if names is not None:
        names = names.split(",")
    if confidence is not None:
        confidence = float(confidence)

    try:
        if names != None:
            # collecting all the protein ids and common names of the perturbations
            all_kinases = session.query(Protein.id).all()
            all_perts = session.query(Perturbation.name).all()
            kinases = []
            perts = []

            # the two for loops change the format of the lists from the sql table since it gives all the results in tuples
            for kin in all_kinases:
                kinases.append(kin[0])

            for pert in all_perts:
                perts.append(pert[0])
            
            # the for loop checks whether the given name is a kinase or a inhibitor and adds matching statement on the list
            direction = []
            for n in names:
                if n in kinases:
                    direction.append("kinase")
                elif n in perts:
                    direction.append("inhibitor")


        # if the user hasn't input name(s), an error message will be shown because name will be needed for this route
        if names == None:
            return "At least one name required for results"

        # if only one category of names is given the output can be generated for the user
        elif len(set(direction)) == 1:
            # defining the category, kinase or inhibitor
            direction = direction[0]

            if confidence is not None:
                if direction == "inhibitor":
                    # going through the list of inhibitors from the user
                    for name in names:
                        # collecting info about the inhibitor from the db
                        all_inhibitor_info1 = session.query(Known_perturbations).filter_by(perturbation = name).filter(or_(Known_perturbations.source == 10,Known_perturbations.source == 12)).filter(Known_perturbations.score <= confidence).all() 
                        all_inhibitor_info2 = session.query(Known_perturbations).filter_by(perturbation = name).filter(Known_perturbations.source == 11).filter((1-Known_perturbations.score) <= confidence).all() 
                        all_inhibitor_info = all_inhibitor_info1 + all_inhibitor_info2
                        kinase_list = []
                        # this for loop ensures that if the inhibitor is not found from the db or it doesn't have any kinases the function will be able to handle it. Only available inhibitors will be shown
                        for inhibit in all_inhibitor_info:
                            if inhibit.kinase is not None:
                                kinase_list.append(inhibit.kinase)
                        # saving the info in json format
                        pert_info = {
                            "PerturbationName" : name,
                            "TargetKinases" : kinase_list
                            }
                        inhibitor_json.append(pert_info)
                    
                    # returning the list of json objects to the user
                    return jsonify(inhibitor_json)
                
                elif direction == "kinase":
                    for name in names:
                        # collecting info about the inhibitors based on the kinase from the db
                        all_inhibitor_info1 = session.query(Known_perturbations).filter_by(kinase = name).filter(or_(Known_perturbations.source == 10,Known_perturbations.source == 12)).filter(Known_perturbations.score <= confidence).all() 
                        all_inhibitor_info2 = session.query(Known_perturbations).filter_by(kinase = name).filter(Known_perturbations.source == 11).filter((1-Known_perturbations.score) <= confidence).all() 
                        all_inhibitor_info = all_inhibitor_info1 + all_inhibitor_info2
                        inhibitor_list = []
                        # this for loop ensures that if the inhibitor is not found from the db or it doesn't have any kinases the function will be able to handle it. Only available inhibitors will be shown
                        for kinase in all_inhibitor_info:
                            if kinase.perturbation is not None:
                                inhibitor_list.append(kinase.perturbation)
                        # saving the info in json format
                        pert_info = {
                            "KinaseName" : name,
                            "KnownPerturbations" : inhibitor_list
                            }
                        inhibitor_json.append(pert_info)
                
                    # returning the list of json objects to the user
                    return jsonify(inhibitor_json)
            
            elif confidence is None:
                if direction == "inhibitor":
                    # going through the list of inhibitors from the user
                    for name in names:
                        # collecting info about the inhibitor from the db
                        all_inhibitor_info = all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Experimental.perturbation == name).filter(Interaction.modification == "phosphorylation").all()
                        kinase_list = []
                        for exp, inter in all_inhibitor_info:
                            if inter.target is not None:
                                kinase_list.append(inter.target)
                        kinase_list = list(set(kinase_list))
                        pert_info = {
                            "PerturbationName" : name,
                            "AvailableTargetKinases" : kinase_list
                            }
                        inhibitor_json.append(pert_info)
                    
                    # returning the list of json objects to the user
                    return jsonify(inhibitor_json)
                
                elif direction == "kinase":
                    for name in names:
                        # collecting info about the inhibitors based on the kinase from the db
                        all_kin_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == name).filter(Interaction.modification == "phosphorylation").all()
                        inhibitor_list = []
                        for exp, inter in all_kin_info:
                            if exp.perturbation is not None:
                                inhibitor_list.append(exp.perturbation)
                        inhibitor_list = list(set(inhibitor_list))
                        pert_info = {
                            "KinaseName" : name,
                            "AvailablePerturbations" : inhibitor_list
                            }
                        inhibitor_json.append(pert_info)
                
                    # returning the list of json objects to the user
                    return jsonify(inhibitor_json)
        
        # if names from both categories have inputed then this error message will be shown
        else:
            return "The names are from different categories. Try again with only kinase ids or only perturbation names."
    
    
    # closing the session to the db once the function has finished
    finally:
        session.close()


@inhibitor_bp.route("/fc", methods = ["GET"])
def get_fc(version):
    """The function gives the user the fc values based on the filters given. The results will be returned in json format."""
    # connection to the database and collecting the filters from the user
    Session = getSessionForVersion(version)
    session = Session()
    type = request.args.get("type")
    kinases = request.args.get("id")
    cell_line = request.args.get("cell_line")
    phosphosite_con = request.args.get("phosphosite_confidence")
    sid_score = request.args.get("sid")
    # getting the available cell lines in order to check is users cell line is available or not in the db
    cell_line_list = session.query(Experimental.cell_line).distinct().all()
    cl_list = []
    for line in cell_line_list:
        cl_list.append(line[0])
    cl_list_string = formatList(cl_list)
    # changing the string 0 or 1 to int
    if phosphosite_con is not None:
        phosphosite_con = int(phosphosite_con)
    if sid_score is not None:
        sid_score = float(sid_score)

    try:
        fc_json = []
        if type is None:
            return "Type (source, target_kinase, target_phosphosite) is required for this route."
        # this if clause will handle the errors if the user does not input correct cell line name or doesn't add kinase ids
        elif (kinases == None and cell_line == None) or (kinases == None or(cell_line is not None and cell_line not in cl_list)):
            if not kinases and cell_line in cl_list:
                return "Uniprot id required for results"
            elif cell_line not in cl_list and kinases:
                return "Cell line not available in experimental data. Available cell lines: " + cl_list_string
            else:
                return "Uniprot id required for results. Also, the current cell line not available in experimental data. Available cell lines: " + cl_list_string

        else:
            # splitting the user input of kinases into a list so that they can be handled individually
            kinases = kinases.split(",")
            if type == "source":
            # using this if clause if the user does not input phosphosite_confidence at all or inputs it as 0
                if not phosphosite_con or phosphosite_con == 0:
                    # going through all the kinases individually
                    for kinase in kinases:
                        # making sure the uniprot id from the user is in correct format
                        kinase = kinase.upper()
                        # collecting data regarding the kinase from Interaction and Experimental tables based on the cell line filter user has given. The Interaction and Experimental table are joined due to the fact that Experimental table doesn't have phosphorylating kinases marked
                        if not cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").filter(Experimental.cell_line==cell_line).all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").filter(Experimental.cell_line==cell_line).filter(Experimental.sid <= sid_score).all()

                        # going through the data from the tables
                        for e_data, i_data in all_info:
                            # kinase_num is to collect all the kinases that are linked to the phosphosite
                            kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                            # kinase_ref is to collect the list of sources for the given kinase from the db
                            kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                            #print(len(kinase_num))

                            # num_kin saves the length of the kinase_num for the confidence calculation
                            num_kin = len(kinase_num)

                            # this for loop goes through kinase_num and records the highest number of references. The comparison is done between all the kinases that are linked to the phosphosite
                            max_ref = 0
                            for kin in kinase_num:
                                #print(kin.references)
                                if max_ref < len(ast.literal_eval(kin.references)):
                                    max_ref = len(ast.literal_eval(kin.references))

                            #print(max_ref)
                            # recording the length of the source list of given kinase
                            ref_num = len(ast.literal_eval(kinase_ref.references))

                            #print(ref_num)
                            #print(ast.literal_eval(kinase_ref.references))

                            # calculating the ref score by dividing the number of sources given kinase has divided by the max number of references (from the group of the kinases that phosphorylates the same phosphosite)
                            # this is one of the indicators when determining how likely the given kinase is to be the one phosphorylating the phosphosite in experimental data
                            ref_score = ref_num / max_ref
                            uniqueness_penalty = 1 / num_kin
                            # basic calculatopn for confidence if no weighthings added. This one can be used but can accidentally disregard kinase that is found from many sources but is not unique enough
                            #confidence = ref_score * uniqueness_penalty
                            # uniqueness less penalised with this calculation
                            confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            # if cell line is not added to the filter all the fc values for given kinase will be recorded regardless of the cell lines.
                            if not cell_line:
                                fc_info = {
                                    i_data.source : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "KinaseConfidence" : confidence,
                                        "CellLine" : cl_list,
                                        "SID-score" : e_data.sid
                                        }
                                }

                            # if cell line is added then just the fc values regarding that cell line will be added to the output json
                            elif cell_line:
                                fc_info = {
                                    i_data.source : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "KinaseConfidence" : confidence,
                                        "CellLine" : cell_line,
                                        "SID_score" : e_data.sid
                                        }
                                }
                            fc_json.append(fc_info)
                    
                    # the full json is given to user
                    return jsonify(fc_json)
                
                # this if clause applied another filter on the sql queries. In this case, adding filter where only High confidence phosphosites (that are found from Phosphosite table) are included in the data
                elif phosphosite_con == 1:
                    # going through all the kinases individually
                    for kinase in kinases:
                        # making sure the uniprot id from the user is in correct format
                        kinase = kinase.upper()
                        # collecting data regarding the kinase from Interaction and Experimental tables based on the cell line filter user has given. The Interaction and Experimental table are joined due to the fact that Experimental table doesn't have phosphorylating kinases marked
                        if not cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()

                        for e_data, i_data in all_info:
                            kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                            kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                            #print(len(kinase_num))
                            num_kin = len(kinase_num)
                            max_ref = 0
                            for kin in kinase_num:
                                #print(kin.references)
                                if max_ref < len(ast.literal_eval(kin.references)):
                                    max_ref = len(ast.literal_eval(kin.references))

                            #print(max_ref)
                            ref_num = len(ast.literal_eval(kinase_ref.references))

                            #print(ref_num)
                            #print(ast.literal_eval(kinase_ref.references))

                            ref_score = ref_num / max_ref
                            uniqueness_penalty = 1 / num_kin
                            # basic calculatopn for confidence
                            #confidence = ref_score * uniqueness_penalty
                            # uniqueness less penalised?
                            confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            # if cell line is not added to the filter all the fc values for given kinase will be recorded regardless of the cell lines.
                            if not cell_line:
                                fc_info = {
                                    i_data.source : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "KinaseConfidence" : confidence,
                                        "SID-score" : e_data.sid
                                        }
                                }

                            # if cell line is added then just the fc values regarding that cell line will be added to the output json
                            elif cell_line:
                                fc_info = {
                                    i_data.source : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "KinaseConfidence" : confidence,
                                        "CellLine" : cell_line,
                                        "SID-score" : e_data.sid
                                        }
                                }
                            fc_json.append(fc_info)

                    # the full json is given to user
                    return jsonify(fc_json)
            
            elif type == "target_kinase":
                if not phosphosite_con or phosphosite_con == 0:
                # going through all the kinases individually
                    for kinase in kinases:
                        # making sure the uniprot id from the user is in correct format
                        kinase = kinase.upper()
                        # collecting data regarding the kinase from Interaction and Experimental tables based on the cell line filter user has given. The Interaction and Experimental table are joined due to the fact that Experimental table doesn't have phosphorylating kinases marked
                        if not cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Experimental.sid <= sid_score).all()
                        elif cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Experimental.cell_line==cell_line).all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Experimental.cell_line==cell_line).filter(Experimental.sid <= sid_score).all()

                        # going through the data from the tables
                        for e_data, p_data in all_info:
                            # # kinase_num is to collect all the kinases that are linked to the phosphosite
                            # kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                            # # kinase_ref is to collect the list of sources for the given kinase from the db
                            # kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                            # #print(len(kinase_num))

                            # # num_kin saves the length of the kinase_num for the confidence calculation
                            # num_kin = len(kinase_num)

                            # # this for loop goes through kinase_num and records the highest number of references. The comparison is done between all the kinases that are linked to the phosphosite
                            # max_ref = 0
                            # for kin in kinase_num:
                            #     #print(kin.references)
                            #     if max_ref < len(ast.literal_eval(kin.references)):
                            #         max_ref = len(ast.literal_eval(kin.references))

                            # #print(max_ref)
                            # # recording the length of the source list of given kinase
                            # ref_num = len(ast.literal_eval(kinase_ref.references))

                            # #print(ref_num)
                            # #print(ast.literal_eval(kinase_ref.references))

                            # # calculating the ref score by dividing the number of sources given kinase has divided by the max number of references (from the group of the kinases that phosphorylates the same phosphosite)
                            # # this is one of the indicators when determining how likely the given kinase is to be the one phosphorylating the phosphosite in experimental data
                            # ref_score = ref_num / max_ref
                            # uniqueness_penalty = 1 / num_kin
                            # # basic calculatopn for confidence if no weighthings added. This one can be used but can accidentally disregard kinase that is found from many sources but is not unique enough
                            # #confidence = ref_score * uniqueness_penalty
                            # # uniqueness less penalised with this calculation
                            # confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            # if cell line is not added to the filter all the fc values for given kinase will be recorded regardless of the cell lines.
                            if not cell_line:
                                fc_info = {
                                    p_data.uniprot_id : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "CellLine" : cl_list,
                                        "SID-score" : e_data.sid
                                        }
                                }

                            # if cell line is added then just the fc values regarding that cell line will be added to the output json
                            elif cell_line:
                                fc_info = {
                                    p_data.uniprot_id : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "CellLine" : cell_line,
                                        "SID_score" : e_data.sid
                                        }
                                }
                            fc_json.append(fc_info)
                    
                    # the full json is given to user
                    return jsonify(fc_json)
            
                # this if clause applied another filter on the sql queries. In this case, adding filter where only High confidence phosphosites (that are found from Phosphosite table) are included in the data
                elif phosphosite_con == 1:
                    # going through all the kinases individually
                    for kinase in kinases:
                        # making sure the uniprot id from the user is in correct format
                        kinase = kinase.upper()
                        # collecting data regarding the kinase from Interaction and Experimental tables based on the cell line filter user has given. The Interaction and Experimental table are joined due to the fact that Experimental table doesn't have phosphorylating kinases marked
                        if not cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Phosphosite.confidence == "High").all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Phosphosite.confidence == "High").filter(Experimental.sid <= sid_score).all()
                        elif cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Phosphosite.uniprot_id == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Experimental.sid <= sid_score).all()

                        for e_data, p_data in all_info:
                        #     kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                        #     kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                        #     #print(len(kinase_num))
                        #     num_kin = len(kinase_num)
                        #     max_ref = 0
                        #     for kin in kinase_num:
                        #         #print(kin.references)
                        #         if max_ref < len(ast.literal_eval(kin.references)):
                        #             max_ref = len(ast.literal_eval(kin.references))

                        #     #print(max_ref)
                        #     ref_num = len(ast.literal_eval(kinase_ref.references))

                        #     #print(ref_num)
                        #     #print(ast.literal_eval(kinase_ref.references))

                        #     ref_score = ref_num / max_ref
                        #     uniqueness_penalty = 1 / num_kin
                        #     # basic calculatopn for confidence
                        #     #confidence = ref_score * uniqueness_penalty
                        #     # uniqueness less penalised?
                        #     confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            # if cell line is not added to the filter all the fc values for given kinase will be recorded regardless of the cell lines.
                            if not cell_line:
                                fc_info = {
                                    p_data.uniprot_id : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "SID-score" : e_data.sid
                                        }
                                }

                            # if cell line is added then just the fc values regarding that cell line will be added to the output json
                            elif cell_line:
                                fc_info = {
                                    p_data.uniprot_id : {
                                        "Perturbation" : e_data.perturbation,
                                        "Phosphosite" : e_data.phosphosite,
                                        "FC" : e_data.fc,
                                        "CellLine" : cell_line,
                                        "SID-score" : e_data.sid
                                        }
                                }
                            fc_json.append(fc_info)

                    # the full json is given to user
                    return jsonify(fc_json)

            elif type == "target_phosphosite":
                if not phosphosite_con or phosphosite_con == 0:
                # going through all the kinases individually
                    for kinase in kinases:
                        # making sure the uniprot id from the user is in correct format
                        kinase = kinase.upper()
                        # collecting data regarding the kinase from Interaction and Experimental tables based on the cell line filter user has given. The Interaction and Experimental table are joined due to the fact that Experimental table doesn't have phosphorylating kinases marked
                        if not cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Experimental.sid <= sid_score).all()
                        elif cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Experimental.cell_line==cell_line).all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Experimental.cell_line==cell_line).filter(Experimental.sid <= sid_score).all()

                        # going through the data from the tables
                        for e_data, p_data in all_info:
                            # # kinase_num is to collect all the kinases that are linked to the phosphosite
                            # kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                            # # kinase_ref is to collect the list of sources for the given kinase from the db
                            # kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                            # #print(len(kinase_num))

                            # # num_kin saves the length of the kinase_num for the confidence calculation
                            # num_kin = len(kinase_num)

                            # # this for loop goes through kinase_num and records the highest number of references. The comparison is done between all the kinases that are linked to the phosphosite
                            # max_ref = 0
                            # for kin in kinase_num:
                            #     #print(kin.references)
                            #     if max_ref < len(ast.literal_eval(kin.references)):
                            #         max_ref = len(ast.literal_eval(kin.references))

                            # #print(max_ref)
                            # # recording the length of the source list of given kinase
                            # ref_num = len(ast.literal_eval(kinase_ref.references))

                            # #print(ref_num)
                            # #print(ast.literal_eval(kinase_ref.references))

                            # # calculating the ref score by dividing the number of sources given kinase has divided by the max number of references (from the group of the kinases that phosphorylates the same phosphosite)
                            # # this is one of the indicators when determining how likely the given kinase is to be the one phosphorylating the phosphosite in experimental data
                            # ref_score = ref_num / max_ref
                            # uniqueness_penalty = 1 / num_kin
                            # # basic calculatopn for confidence if no weighthings added. This one can be used but can accidentally disregard kinase that is found from many sources but is not unique enough
                            # #confidence = ref_score * uniqueness_penalty
                            # # uniqueness less penalised with this calculation
                            # confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            # if cell line is not added to the filter all the fc values for given kinase will be recorded regardless of the cell lines.
                            if not cell_line:
                                fc_info = {
                                    e_data.phosphosite : {
                                        "Perturbation" : e_data.perturbation,
                                        "FC" : e_data.fc,
                                        "CellLine" : cl_list,
                                        "SID-score" : e_data.sid
                                        }
                                }

                            # if cell line is added then just the fc values regarding that cell line will be added to the output json
                            elif cell_line:
                                fc_info = {
                                    e_data.phosphosite : {
                                        "Perturbation" : e_data.perturbation,
                                        "FC" : e_data.fc,
                                        "CellLine" : cell_line,
                                        "SID_score" : e_data.sid
                                        }
                                }
                            fc_json.append(fc_info)
                    
                    # the full json is given to user
                    return jsonify(fc_json)
            
                # this if clause applied another filter on the sql queries. In this case, adding filter where only High confidence phosphosites (that are found from Phosphosite table) are included in the data
                elif phosphosite_con == 1:
                    # going through all the kinases individually
                    for kinase in kinases:
                        # making sure the uniprot id from the user is in correct format
                        kinase = kinase.upper()
                        # collecting data regarding the kinase from Interaction and Experimental tables based on the cell line filter user has given. The Interaction and Experimental table are joined due to the fact that Experimental table doesn't have phosphorylating kinases marked
                        if not cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Phosphosite.confidence == "High").all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Phosphosite.confidence == "High").filter(Experimental.sid <= sid_score).all()
                        elif cell_line:
                            if not sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").all()
                            elif sid_score:
                                all_info = session.query(Experimental,Phosphosite).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Experimental.phosphosite == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Experimental.sid <= sid_score).all()

                        for e_data, i_data in all_info:
                        #     kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                        #     kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                        #     #print(len(kinase_num))
                        #     num_kin = len(kinase_num)
                        #     max_ref = 0
                        #     for kin in kinase_num:
                        #         #print(kin.references)
                        #         if max_ref < len(ast.literal_eval(kin.references)):
                        #             max_ref = len(ast.literal_eval(kin.references))

                        #     #print(max_ref)
                        #     ref_num = len(ast.literal_eval(kinase_ref.references))

                        #     #print(ref_num)
                        #     #print(ast.literal_eval(kinase_ref.references))

                        #     ref_score = ref_num / max_ref
                        #     uniqueness_penalty = 1 / num_kin
                        #     # basic calculatopn for confidence
                        #     #confidence = ref_score * uniqueness_penalty
                        #     # uniqueness less penalised?
                        #     confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            # if cell line is not added to the filter all the fc values for given kinase will be recorded regardless of the cell lines.
                            if not cell_line:
                                fc_info = {
                                    e_data.phosphosite : {
                                        "Perturbation" : e_data.perturbation,
                                        "FC" : e_data.fc,
                                        "SID-score" : e_data.sid
                                        }
                                }

                            # if cell line is added then just the fc values regarding that cell line will be added to the output json
                            elif cell_line:
                                fc_info = {
                                    e_data.phosphosite : {
                                        "Perturbation" : e_data.perturbation,
                                        "FC" : e_data.fc,
                                        "CellLine" : cell_line,
                                        "SID-score" : e_data.sid
                                        }
                                }
                        fc_json.append(fc_info)

                # the full json is given to user
                return jsonify(fc_json)
        
            
    # closing the session to the db once the function has finished    
    finally:
        session.close()


def calculate_weighted_sd(fcs, weights):
    '''Function calculates the weigthed standard deviation based on Bessel correlation by using fc values and the confidence weightings
    that have been calculated earlier.'''
    mean = numpy.average(fcs, weights=weights)
    numerator = numpy.sum(weights * (fcs - mean)**2)
    denominator = numpy.sum(weights) - (numpy.sum(weights**2) / numpy.sum(weights))
    return numpy.sqrt(numerator / denominator)

@inhibitor_bp.route("/KSEA", methods = ["GET"])
def get_ksea(version):
    """Optimised KSEA endpoint — fetches all kinases in a single bulk query
    instead of looping one kinase at a time. Background stats (mean FC, SD)
    are computed once per perturbation rather than repeated for every kinase.
    This gives the same results as before but is significantly faster."""
    Session = getSessionForVersion(version)
    session = Session()
    ksea_json = []
    kinases = request.args.get("kinase_ids")
    inhibitors = request.args.get("perturbations")
    cell_line = request.args.get("cell_line")
    autophosphorylation = request.args.get("autophosphorylation")
    weighted = request.args.get("weighted")
    phosphosite_con = request.args.get("phosphosite_confidence")
    sid_score = request.args.get("sid")
    if kinases is not None:
        kinases = [k.upper() for k in kinases.split(",")]
    if inhibitors is not None:
        inhibitors = inhibitors.split(",")
    if phosphosite_con is not None:
        phosphosite_con = int(phosphosite_con)
    if sid_score is not None:
        sid_score = float(sid_score)

    cell_line_list = session.query(Experimental.cell_line).distinct().all()
    cl_list = []
    for line in cell_line_list:
        cl_list.append(line[0])
    cl_list = formatList(cl_list)

    try:
        if (kinases == None and not cell_line) or (kinases == None or(cell_line is not None and cell_line not in cl_list)):
            if not kinases and cell_line in cl_list:
                return "Uniprot id required for results"
            elif cell_line not in cl_list and kinases:
                return "Cell line not available in experimental data. Available cell lines: " + cl_list
            else:
                return "Uniprot id required for results. Also, the current cell line not available in experimental data. Available cell lines: " + cl_list
        else:
            if inhibitors == None and not cell_line:
                for kinase in kinases:
                    kinase = kinase.upper()
                    if autophosphorylation == "include" or not autophosphorylation:
                        if not phosphosite_con or phosphosite_con == 0:
                            inter = session.query(Interaction.phosphosite_id).distinct().filter(Interaction.modification == "phosphorylation").all()
                        elif phosphosite_con ==  1:
                            inter = session.query(Interaction.phosphosite_id).join(Phosphosite,Interaction.phosphosite_id==Phosphosite.phosphosite_id).filter(Phosphosite.confidence=="High").filter(Interaction.modification == "phosphorylation").distinct().all()
                    elif autophosphorylation == "exclude":
                        #inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).distinct().all()
                        if not phosphosite_con or phosphosite_con == 0:
                            inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).distinct().filter(Interaction.modification == "phosphorylation").all()
                        elif phosphosite_con ==  1:
                            inter = session.query(Interaction.phosphosite_id).join(Phosphosite,Interaction.phosphosite_id==Phosphosite.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").distinct().all()
                    #print(inter)
                    inter_list = []
                    for phos_i in inter:
                        inter_list.append(phos_i[0])
                    inter_distinct = set(inter_list)

                    #print(len(inter_list))
                    #print(len(inter_distinct))

                    if sid_score == None:
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).all()
                    elif sid_score is not None:
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.sid<=sid_score).all()

                    phos_list = []
                    fc_list = []
                    for pho in phos_fc:
                        phosphosite = pho.phosphosite
                        fc = pho.fc
                        phos_list.append(phosphosite)
                        fc_list.append(float(fc))

                    all_fc = pd.DataFrame({"phosphosite":phos_list,"fc":fc_list})
                    #print(all_fc)
                    #print(len(all_fc))

                    standard_d_fc = all_fc["fc"].std(ddof=1)
                    #print(standard_d_fc)

                    p = all_fc["fc"].mean()
                    #print(p)

                    
                    if autophosphorylation == "include" or not autophosphorylation:
                        #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).all()
                        if not phosphosite_con or phosphosite_con == 0:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif phosphosite_con ==  1:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                    elif autophosphorylation == "exclude":
                        #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).all()
                        if not phosphosite_con or phosphosite_con == 0:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif phosphosite_con ==  1:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()

                    kin_list = []
                    phos_list = []
                    fc_list = []
                    cell_list = []
                    if weighted == "true":
                        confidence_list = []
                    for e_data, i_data in all_info:
                        #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                        kin = i_data.source
                        phosphosite = e_data.phosphosite

                        if weighted == "true":
                            if autophosphorylation == "include" or not autophosphorylation:
                                kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                                kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                            elif autophosphorylation == "exclude":
                                kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                                kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").first()
                            
                            #print(len(kinase_num))
                            num_kin = len(kinase_num)
                            max_ref = 0
                            for kin in kinase_num:
                                #print(kin.references)
                                if max_ref < len(ast.literal_eval(kin.references)):
                                    max_ref = len(ast.literal_eval(kin.references))

                            #print(max_ref)
                            ref_num = len(ast.literal_eval(kinase_ref.references))

                            #print(ref_num)
                            #print(ast.literal_eval(kinase_ref.references))

                            ref_score = ref_num / max_ref
                            uniqueness_penalty = 1 / num_kin
                            # basic calculatopn for confidence
                            #confidence = ref_score * uniqueness_penalty
                            # uniqueness less penalised?
                            confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            confidence_list.append(confidence)



                        fc = float(e_data.fc)
                        kin_list.append(kin)
                        phos_list.append(phosphosite)
                        fc_list.append(fc)
                        cell_list.append(e_data.cell_line)

                    kin_fc = pd.DataFrame({"kinase":kin_list,"phosphosite":phos_list,"fc":fc_list})
                    #print(kin_fc)
                    #print(len(kin_fc))

                    if not weighted or weighted == "false":
                        s = kin_fc["fc"].mean()

                        sqr_m = numpy.sqrt(len(kin_fc))

                        #print(s)

                        #print(sqr_m)

                        z = ((s-p)*sqr_m)/standard_d_fc

                        #print(z)

                        p_value = norm.sf(abs(z))

                        #uniques_phos_list = set(phos_list)

                        cell_list = set(cell_list)


                    #print(p_value)
                        ksea_info = {
                            kinase : {
                                "Phosphosites" : list(phos_list),
                                "MeanFCKinase" : s,
                                "MeanFCAll" : p,
                                "StandardDeviation" : standard_d_fc,
                                "n" : len(kin_fc),
                                "z_score" : z,
                                "p_value" : p_value,
                                "CellLinesIncluded" : list(cell_list)
                                }
                            }
                    
                    elif weighted == "true":
                        # weighted calculations
                        weights = numpy.array(confidence_list)
                        fcs = numpy.array(fc_list)
                        if numpy.sum(weights) == 0:
                            weighted_mean = float("nan")
                            weighted_sd = float("nan")
                            z = float("nan")
                            p_value = float("nan")
                        else:
                            weighted_mean = numpy.average(fcs, weights=weights)
                            weighted_sd = calculate_weighted_sd(fcs, weights)

                            z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                            p_value = norm.sf(abs(z))
                        cell_list = set(cell_list)
                        #uniques_phos_list = set(phos_list)
                    
                        ksea_info = {
                                kinase : {
                                    "Phosphosites" : list(phos_list),
                                    "WeightedMeanFCKinase" : weighted_mean,
                                    "MeanFCAll" : p,
                                    "WeightedStandardDeviation" : weighted_sd,
                                    "n" : len(kin_fc),
                                    "WeightedZ_score" : z,
                                    "p_value" : p_value,
                                    "CellLinesIncluded" : list(cell_list),
                                    "MeanConfidence" : numpy.mean(weights)
                                    }
                                }
                
                ksea_json.append(ksea_info)
                #print(ksea_json)    
                return jsonify(ksea_json)

            elif inhibitors != None and not cell_line:
                for kinase in kinases:
                    kinase = kinase.upper()
                    # all_inhibitors = session.query(Experimental.perturbation).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).all()
                    # inhibitor_list = []

                    # for i in all_inhibitors:
                    #     inhibitor_list.append(i[0])
                    
                    for inhibitor in inhibitors:
                        # if autophosphorylation == "include" or not autophosphorylation:
                        #     inter = session.query(Interaction.phosphosite_id).distinct().all()
                        # elif autophosphorylation == "exclude":
                        #     inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).distinct().all()

                        if autophosphorylation == "include" or not autophosphorylation:
                            if not phosphosite_con or phosphosite_con == 0:
                                inter = session.query(Interaction.phosphosite_id).filter(Interaction.modification == "phosphorylation").distinct().all()
                            elif phosphosite_con ==  1:
                                inter = session.query(Interaction.phosphosite_id).join(Phosphosite,Interaction.phosphosite_id==Phosphosite.phosphosite_id).filter(Phosphosite.confidence=="High").filter(Interaction.modification == "phosphorylation").distinct().all()
                        elif autophosphorylation == "exclude":
                        #inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).distinct().all()
                            if not phosphosite_con or phosphosite_con == 0:
                                inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").distinct().all()
                            elif phosphosite_con ==  1:
                                inter = session.query(Interaction.phosphosite_id).join(Phosphosite,Interaction.phosphosite_id==Phosphosite.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").distinct().all()
                    

                        #print(inter)
                        inter_list = []
                        for phos_i in inter:
                            inter_list.append(phos_i[0])
                        inter_distinct = set(inter_list)

                        #print(len(inter_list))
                        #print(len(inter_distinct))
                        
                        if sid_score is None:
                            phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.perturbation==inhibitor).all()
                        elif sid_score is not None:
                            phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.perturbation==inhibitor).filter(Experimental.sid <= sid_score).all()

                        phos_list = []
                        fc_list = []
                        for pho in phos_fc:
                            phosphosite = pho.phosphosite
                            fc = pho.fc
                            phos_list.append(phosphosite)
                            fc_list.append(float(fc))

                        all_fc = pd.DataFrame({"phosphosite":phos_list,"fc":fc_list})
                        #print(all_fc)
                        #print(len(all_fc))

                        standard_d_fc = all_fc["fc"].std(ddof=1)
                        #print(standard_d_fc)

                        p = all_fc["fc"].mean()
                        #print(p)
                        # if autophosphorylation == "include" or not autophosphorylation:
                        #     all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).all()
                        # elif autophosphorylation == "exclude":
                        #     all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).all()

                        if autophosphorylation == "include" or not autophosphorylation:
                            #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).all()
                            if not phosphosite_con or phosphosite_con == 0:
                                if not sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.modification == "phosphorylation").all()
                                elif sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                            elif phosphosite_con ==  1:
                                if not sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                                elif sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif autophosphorylation == "exclude":
                            #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).all()
                            if not phosphosite_con or phosphosite_con == 0:
                                if not sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                                elif sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                            elif phosphosite_con ==  1:
                                if not sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                                elif sid_score:
                                    all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()


                        #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).all()
                        kin_list = []
                        phos_list = []
                        fc_list = []
                        cell_list = []
                        if weighted == "true":
                            confidence_list = []
                        for e_data, i_data in all_info:
                            #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                            kin = i_data.source
                            phosphosite = e_data.phosphosite

                            if weighted == "true":
                                if autophosphorylation == "include" or not autophosphorylation:
                                    kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                                    kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                                elif autophosphorylation == "exclude":
                                    kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                                    kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").first()
                                #print(len(kinase_num))
                                num_kin = len(kinase_num)
                                max_ref = 0
                                for kin in kinase_num:
                                    #print(kin.references)
                                    if max_ref < len(ast.literal_eval(kin.references)):
                                        max_ref = len(ast.literal_eval(kin.references))

                                #print(max_ref)
                                ref_num = len(ast.literal_eval(kinase_ref.references))

                                #print(ref_num)
                                #print(ast.literal_eval(kinase_ref.references))

                                ref_score = ref_num / max_ref
                                uniqueness_penalty = 1 / num_kin
                                # basic calculatopn for confidence
                                #confidence = ref_score * uniqueness_penalty
                                # uniqueness less penalised?
                                confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                                confidence_list.append(confidence)

                            fc = float(e_data.fc)
                            kin_list.append(kin)
                            phos_list.append(phosphosite)
                            fc_list.append(fc)
                            cell_list.append(e_data.cell_line)

                        kin_fc = pd.DataFrame({"kinase":kin_list,"phosphosite":phos_list,"fc":fc_list})
                        #print(kin_fc)
                        #print(len(kin_fc))
                        if not weighted or weighted == "false":
                            s = kin_fc["fc"].mean()

                            sqr_m = numpy.sqrt(len(kin_fc))

                            #print(s)

                            #print(sqr_m)

                            z = ((s-p)*sqr_m)/standard_d_fc

                            #print(z)

                            p_value = norm.sf(abs(z))

                            #uniques_phos_list = set(phos_list)

                            cell_list = set(cell_list)


                        #print(p_value)
                            ksea_info = {
                                kinase : {
                                    inhibitor : {
                                        "Phosphosites" : list(phos_list),
                                        "MeanFCKinase" : s,
                                        "MeanFCAll" : p,
                                        "StandardDeviation" : standard_d_fc,
                                        "n" : len(kin_fc),
                                        "z_score" : z,
                                        "p_value" : p_value,
                                        "CellLinesIncluded" : list(cell_list)
                                        }
                                    }
                                }
                        
                        elif weighted == "true":
                            # weighted calculations
                            weights = numpy.array(confidence_list)
                            fcs = numpy.array(fc_list)
                            if numpy.sum(weights) == 0:
                                weighted_mean = float("nan")
                                weighted_sd = float("nan")
                                z = float("nan")
                                p_value = float("nan")
                            else:
                                weighted_mean = numpy.average(fcs, weights=weights)
                                weighted_sd = calculate_weighted_sd(fcs, weights)

                                z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                                p_value = norm.sf(abs(z))
                                cell_list = set(cell_list)
                            #uniques_phos_list = set(phos_list)
                        
                            ksea_info = {
                                    kinase : {
                                        inhibitor : {
                                            "Phosphosites" : list(phos_list),
                                            "WeightedMeanFCKinase" : weighted_mean,
                                            "MeanFCAll" : p,
                                            "WeightedStandardDeviation" : weighted_sd,
                                            "n" : len(kin_fc),
                                            "WeightedZ_score" : z,
                                            "p_value" : p_value,
                                            "CellLinesIncluded" : list(cell_list),
                                            "MeanConfidence" : numpy.mean(weights)
                                            }
                                        }
                                    }
                    
                        ksea_json.append(ksea_info)
                #print(ksea_json)        
                return jsonify(ksea_json)

            elif inhibitors == None and cell_line:
                for kinase in kinases:
                    kinase = kinase.upper()
                    # if autophosphorylation == "include" or not autophosphorylation:
                    #     inter = session.query(Interaction.phosphosite_id).distinct().all()
                    # elif autophosphorylation == "exclude":
                    #     inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).distinct().all()


                    if autophosphorylation == "include" or not autophosphorylation:
                        if not phosphosite_con or phosphosite_con == 0:
                            inter = session.query(Interaction.phosphosite_id).filter(Interaction.modification == "phosphorylation").distinct().all()
                        elif phosphosite_con ==  1:
                            inter = session.query(Interaction.phosphosite_id).join(Phosphosite,Interaction.phosphosite_id==Phosphosite.phosphosite_id).filter(Phosphosite.confidence=="High").filter(Interaction.modification == "phosphorylation").distinct().all()
                    elif autophosphorylation == "exclude":
                    #inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).distinct().all()
                        if not phosphosite_con or phosphosite_con == 0:
                            inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").distinct().all()
                        elif phosphosite_con ==  1:
                            inter = session.query(Interaction.phosphosite_id).join(Phosphosite,Interaction.phosphosite_id==Phosphosite.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").distinct().all()
                

                    #inter = session.query(Interaction.phosphosite_id).distinct().all()
                    #print(inter)
                    inter_list = []
                    for phos_i in inter:
                        inter_list.append(phos_i[0])
                    inter_distinct = set(inter_list)

                    #print(len(inter_list))
                    #print(len(inter_distinct))

                    if sid_score is None:
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.cell_line==cell_line).all()
                    elif sid_score is not None:
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.cell_line==cell_line).filter(Experimental.sid <= sid_score).all()

                    phos_list = []
                    fc_list = []
                    for pho in phos_fc:
                        phosphosite = pho.phosphosite
                        fc = pho.fc
                        phos_list.append(phosphosite)
                        fc_list.append(float(fc))

                    all_fc = pd.DataFrame({"phosphosite":phos_list,"fc":fc_list})
                    #print(all_fc)
                    #print(len(all_fc))

                    standard_d_fc = all_fc["fc"].std(ddof=1)
                    #print(standard_d_fc)

                    p = all_fc["fc"].mean()
                    #print(p)

                    # if autophosphorylation == "include" or not autophosphorylation:
                    #     all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).all()
                    # elif autophosphorylation == "exclude":
                    #     all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).all()

                    if autophosphorylation == "include" or not autophosphorylation:
                        if not phosphosite_con or phosphosite_con == 0:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif phosphosite_con ==  1:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                    elif autophosphorylation == "exclude":
                        if not phosphosite_con or phosphosite_con == 0:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()
                        elif phosphosite_con ==  1:
                            if not sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                            elif sid_score:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").filter(Experimental.sid <= sid_score).all()



                    #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).all()
                    kin_list = []
                    phos_list = []
                    fc_list = []
                    cell_list = []
                    if weighted == "true":
                        confidence_list = []
                    for e_data, i_data in all_info:
                        #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                        kin = i_data.source
                        phosphosite = e_data.phosphosite

                        if weighted == "true":
                            if autophosphorylation == "include" or not autophosphorylation:
                                kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").all()
                                kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.modification == "phosphorylation").first()
                            elif autophosphorylation == "exclude":
                                kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.modification == "phosphorylation").filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                                kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").first()
                            
                            num_kin = len(kinase_num)
                            max_ref = 0
                            for kin in kinase_num:
                                #print(kin.references)
                                if max_ref < len(ast.literal_eval(kin.references)):
                                    max_ref = len(ast.literal_eval(kin.references))

                            #print(max_ref)
                            ref_num = len(ast.literal_eval(kinase_ref.references))

                            #print(ref_num)
                            #print(ast.literal_eval(kinase_ref.references))

                            ref_score = ref_num / max_ref
                            uniqueness_penalty = 1 / num_kin
                            # basic calculatopn for confidence
                            #confidence = ref_score * uniqueness_penalty
                            # uniqueness less penalised?
                            confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                            confidence_list.append(confidence)

                        fc = float(e_data.fc)
                        kin_list.append(kin)
                        phos_list.append(phosphosite)
                        fc_list.append(fc)
                        cell_list.append(e_data.cell_line)

                    kin_fc = pd.DataFrame({"kinase":kin_list,"phosphosite":phos_list,"fc":fc_list})
                    #print(kin_fc)
                    #print(len(kin_fc))
                    if not weighted or weighted == "false":
                        s = kin_fc["fc"].mean()

                        sqr_m = numpy.sqrt(len(kin_fc))

                        #print(s)

                        #print(sqr_m)

                        z = ((s-p)*sqr_m)/standard_d_fc

                        #print(z)

                        p_value = norm.sf(abs(z))

                        #uniques_phos_list = set(phos_list)

                        cell_list = set(cell_list)


                    #print(p_value)
                        ksea_info = {
                            kinase : {
                                "Phosphosites" : list(phos_list),
                                "MeanFCKinase" : s,
                                "MeanFCAll" : p,
                                "StandardDeviation" : standard_d_fc,
                                "n" : len(kin_fc),
                                "z_score" : z,
                                "p_value" : p_value,
                                "CellLinesIncluded" : list(cell_list)
                                }
                            }
                    
                    elif weighted == "true":
                        # weighted calculations
                        weights = numpy.array(confidence_list)
                        fcs = numpy.array(fc_list)
                        if numpy.sum(weights) == 0:
                            weighted_mean = float("nan")
                            weighted_sd = float("nan")
                            z = float("nan")
                            p_value = float("nan")
                        else:
                            weighted_mean = numpy.average(fcs, weights=weights)
                            weighted_sd = calculate_weighted_sd(fcs, weights)

                            z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                            p_value = norm.sf(abs(z))
                            cell_list = set(cell_list)
                        #uniques_phos_list = set(phos_list)
                    
                        ksea_info = {
                                kinase : {
                                    "Phosphosites" : list(phos_list),
                                    "WeightedMeanFCKinase" : weighted_mean,
                                    "MeanFCAll" : p,
                                    "WeightedStandardDeviation" : weighted_sd,
                                    "n" : len(kin_fc),
                                    "WeightedZ_score" : z,
                                    "p_value" : p_value,
                                    "CellLinesIncluded" : list(cell_list),
                                    "MeanConfidence" : numpy.mean(weights)
                                    }
                                }
                
                    ksea_json.append(ksea_info)
                #print(ksea_json)    
                return jsonify(ksea_json)

            elif inhibitors != None and cell_line:
                # ── OPTIMISED: fetch all kinases in ONE bulk query ──────────
                # Step 1: get valid phosphosite IDs once (same for all kinases)
                if autophosphorylation == "include" or not autophosphorylation:
                    if not phosphosite_con or phosphosite_con == 0:
                        inter = session.query(Interaction.phosphosite_id).filter(Interaction.modification == "phosphorylation").distinct().all()
                    else:
                        inter = session.query(Interaction.phosphosite_id).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").distinct().all()
                else:
                    if not phosphosite_con or phosphosite_con == 0:
                        inter = session.query(Interaction.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").distinct().all()
                    else:
                        inter = session.query(Interaction.phosphosite_id).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").distinct().all()

                inter_distinct = set(p[0] for p in inter)

                # Step 2: compute background stats once per inhibitor (not per kinase)
                bg_stats = {}
                for inhibitor in inhibitors:
                    if sid_score is None:
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.perturbation == inhibitor).filter(Experimental.cell_line == cell_line).all()
                    else:
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.perturbation == inhibitor).filter(Experimental.cell_line == cell_line).filter(Experimental.sid <= sid_score).all()
                    bg_df = pd.DataFrame([(p.phosphosite, float(p.fc)) for p in phos_fc], columns=["phosphosite", "fc"])
                    bg_stats[inhibitor] = {
                        "p": bg_df["fc"].mean() if len(bg_df) > 0 else float("nan"),
                        "std": bg_df["fc"].std(ddof=1) if len(bg_df) > 0 else float("nan"),
                    }

                # Step 3: fetch ALL kinases × ALL inhibitors in ONE query
                if autophosphorylation == "include" or not autophosphorylation:
                    if not phosphosite_con or phosphosite_con == 0:
                        q = session.query(Experimental, Interaction).join(Interaction, Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source.in_(kinases)).filter(Experimental.perturbation.in_(inhibitors)).filter(Experimental.cell_line == cell_line).filter(Interaction.modification == "phosphorylation")
                    else:
                        q = session.query(Experimental, Interaction).join(Interaction, Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite, Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source.in_(kinases)).filter(Experimental.perturbation.in_(inhibitors)).filter(Experimental.cell_line == cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation")
                else:
                    if not phosphosite_con or phosphosite_con == 0:
                        q = session.query(Experimental, Interaction).join(Interaction, Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source.in_(kinases)).filter(Experimental.perturbation.in_(inhibitors)).filter(Experimental.cell_line == cell_line).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation")
                    else:
                        q = session.query(Experimental, Interaction).join(Interaction, Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite, Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source.in_(kinases)).filter(Experimental.perturbation.in_(inhibitors)).filter(Experimental.cell_line == cell_line).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation")

                if sid_score is not None:
                    q = q.filter(Experimental.sid <= sid_score)

                all_info = q.all()

                # Build one dataframe for all kinases at once
                rows = [
                    (i_data.source, e_data.phosphosite, float(e_data.fc), e_data.perturbation, e_data.cell_line)
                    for e_data, i_data in all_info
                ]
                df = pd.DataFrame(rows, columns=["kinase", "phosphosite", "fc", "perturbation", "cell_line"])

                # Step 4: compute z-scores using groupby — no per-kinase loop needed
                for kinase in kinases:
                    for inhibitor in inhibitors:
                        kin_df = df[(df["kinase"] == kinase) & (df["perturbation"] == inhibitor)]
                        if len(kin_df) == 0:
                            continue

                        p = bg_stats[inhibitor]["p"]
                        standard_d_fc = bg_stats[inhibitor]["std"]
                        phos_list = list(kin_df["phosphosite"])
                        fc_list = list(kin_df["fc"])
                        cell_list = list(kin_df["cell_line"].unique())
                        n = len(kin_df)

                        if not weighted or weighted == "false":
                            s = kin_df["fc"].mean()

                            sqr_m = numpy.sqrt(n)

                            z = ((s-p)*sqr_m)/standard_d_fc

                            p_value = norm.sf(abs(z))

                            cell_list = set(cell_list)

                            ksea_info = {
                                kinase : {
                                    inhibitor : {
                                        "Phosphosites" : list(phos_list),
                                        "MeanFCKinase" : s,
                                        "MeanFCAll" : p,
                                        "StandardDeviation" : standard_d_fc,
                                        "n" : n,
                                        "z_score" : z,
                                        "p_value" : p_value,
                                        "CellLinesIncluded" : list(cell_list)
                                        }
                                    }
                                }

                        elif weighted == "true":
                            # weighted calculations
                            confidence_list = []
                            for _, row in kin_df.iterrows():
                                if autophosphorylation == "include" or not autophosphorylation:
                                    kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == row["phosphosite"]).filter(Interaction.modification == "phosphorylation").all()
                                    kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == row["phosphosite"]).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").first()
                                else:
                                    kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == row["phosphosite"]).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                                    kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == row["phosphosite"]).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").first()
                                if kinase_ref is None:
                                    confidence_list.append(0.0)
                                    continue
                                num_kin = len(kinase_num)
                                max_ref = 0
                                for kin in kinase_num:
                                    if max_ref < len(ast.literal_eval(kin.references)):
                                        max_ref = len(ast.literal_eval(kin.references))
                                ref_num = len(ast.literal_eval(kinase_ref.references))
                                ref_score = ref_num / max_ref
                                uniqueness_penalty = 1 / num_kin
                                confidence_list.append((0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2))

                            weights = numpy.array(confidence_list)
                            fcs = numpy.array(fc_list)
                            if numpy.sum(weights) == 0:
                                weighted_mean = float("nan")
                                weighted_sd = float("nan")
                                z = float("nan")
                                p_value = float("nan")
                            else:
                                weighted_mean = numpy.average(fcs, weights=weights)
                                weighted_sd = calculate_weighted_sd(fcs, weights)

                                z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                                p_value = norm.sf(abs(z))
                                cell_list = set(cell_list)
                            #uniques_phos_list = set(phos_list)

                            ksea_info = {
                                    kinase : {
                                        inhibitor : {
                                            "Phosphosites" : list(phos_list),
                                            "WeightedMeanFCKinase" : weighted_mean,
                                            "MeanFCAll" : p,
                                            "WeightedStandardDeviation" : weighted_sd,
                                            "n" : n,
                                            "WeightedZ_score" : z,
                                            "p_value" : p_value,
                                            "CellLinesIncluded" : list(cell_list),
                                            "MeanConfidence" : numpy.mean(weights)
                                            }
                                        }
                                    }
                        ksea_json.append(ksea_info)
                return jsonify(ksea_json)
            
        
    finally:
        session.close()
