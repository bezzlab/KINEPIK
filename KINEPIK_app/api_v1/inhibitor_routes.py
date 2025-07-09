from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
from KINEPIK_app.database_con import session_local
import ast
from scipy.stats import norm
import numpy
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# naming the blueprint
inhibitor_bp = Blueprint("inhibitor_bp", __name__, url_prefix="/api/inhibitor")

def formatList(list_form):
    if not list_form:
        return ""
    elif len(list_form) == 1:
        return list_form[0]
    elif len(list_form) == 2:
        return " and ".join(list_form)
    else:
        return ", ".join(list_form[:-1]) + " and " + str(list_form[-1])

@inhibitor_bp.route("/info", methods = ["GET"])
def instructions():
    '''This function hands out the user info about the possible parameters that can be given to the inhibitor queries. 
    The return returns the info from string/text to html which is then displayed to the user'''
    # information about the parameters
    session = session_local()
    cell_line_list = session.query(Experimental.cell_line).distinct().all()
    cl_list = []
    for line in cell_line_list:
        cl_list.append(line[0])

    cell_lines = "cell_line : available cell lines curretly are :" + formatList(cl_list)
    # the format of the html page
    info = f"""
    <html>
        <body>
            <p>{cell_lines}</p>
        </body>
    </html>
    """

    # return changes the format of the string to html that will be displayed to the user
    return Response(info, mimetype = "text/html")

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
                            "Actio" : inhibitor.action
                            }
                    pert_json.append(inhibitor_info)
            return jsonify(pert_json)
        
        elif pert_type and "small_molecule" in pert_type:
            all_inhibitor_info = session.query(Perturbation).filter(Perturbation.type =="small molecule").all()
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
            return jsonify(pert_json)
        
        elif pert_type and "knockout" in pert_type:
            all_inhibitor_info = session.query(Perturbation).filter(Perturbation.type =="CRISPR knockout").all()
            for inhibitor in all_inhibitor_info:
                inhibitor_info = {
                                "PerturbationName" : inhibitor.name,
                                "Gene" : inhibitor.gene,
                                "Type" : inhibitor.type,
                                "Action" : inhibitor.action
                                }
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
                pert_info = {
                    "PerturbationName" : inhibitor,
                    "TargetKinases" : kinase_list
                    }
                inhibitor_json.append(pert_info)
            
            return jsonify(inhibitor_json)
        
    finally:
        session.close()


@inhibitor_bp.route("/fc", methods = ["GET"])
def get_fc():
    # connection to the database
    session = session_local()
    fc_json = []
    kinases = request.args.get("kinase_id")
    cell_line = request.args.get("cell_line")
    phosphosite_con = request.args.get("phosphosite_confidence")
    if kinases is not None:
        kinases = kinases.split(",")
    cell_line_list = session.query(Experimental.cell_line).distinct().all()
    cl_list = []
    for line in cell_line_list:
        cl_list.append(line[0])
    cl_list = formatList(cl_list)
    if phosphosite_con is not None:
        phosphosite_con = int(phosphosite_con)

    try:
        if (kinases == None and not cell_line) or (kinases == None or(cell_line is not None and cell_line not in cl_list)):
            if not kinases and cell_line in cl_list:
                return "Uniprot id required for results"
            elif cell_line not in cl_list and kinases:
                return "Cell line not available in experimental data. Available cell lines: " + cl_list
            else:
                return "Uniprot id required for results. Also, the current cell line not available in experimental data. Available cell lines: " + cl_list

        else:
            if not phosphosite_con or phosphosite_con == 0:
                for kinase in kinases:
                    kinase = kinase.upper()
                    if not cell_line:
                        all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").all()
                    elif cell_line:
                        all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").filter(Experimental.cell_line==cell_line).all()

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

                        if not cell_line:
                            fc_info = {
                                i_data.source : {
                                    "Inhibitor" : e_data.perturbation,
                                    "Phosphosite" : e_data.phosphosite,
                                    "FC" : e_data.fc,
                                    "KinaseConfidence" : confidence
                                    }
                            }
                        elif cell_line:
                            fc_info = {
                                i_data.source : {
                                    "Inhibitor" : e_data.perturbation,
                                    "Phosphosite" : e_data.phosphosite,
                                    "FC" : e_data.fc,
                                    "KinaseConfidence" : confidence,
                                    "CellLine" : cell_line
                                    }
                            }
                        fc_json.append(fc_info)
                
                return jsonify(fc_json)
            

            elif phosphosite_con == 1:
                for kinase in kinases:
                    kinase = kinase.upper()
                    if not cell_line:
                        all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                    elif cell_line:
                        all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()

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

                        if not cell_line:
                            fc_info = {
                                i_data.source : {
                                    "Inhibitor" : e_data.perturbation,
                                    "Phosphosite" : e_data.phosphosite,
                                    "FC" : e_data.fc,
                                    "KinaseConfidence" : confidence
                                    }
                            }
                        elif cell_line:
                            fc_info = {
                                i_data.source : {
                                    "Inhibitor" : e_data.perturbation,
                                    "Phosphosite" : e_data.phosphosite,
                                    "FC" : e_data.fc,
                                    "KinaseConfidence" : confidence,
                                    "CellLine" : cell_line
                                    }
                            }
                        fc_json.append(fc_info)
                
                return jsonify(fc_json)
        
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
def get_ksea():
    session = session_local()
    ksea_json = []
    kinases = request.args.get("kinase_id")
    inhibitors = request.args.get("inhibitors")
    cell_line = request.args.get("cell_line")
    autophosphorylation = request.args.get("autophosphorylation")
    weighted = request.args.get("weighted")
    phosphosite_con = request.args.get("phosphosite_confidence")
    if kinases is not None:
        kinases = kinases.split(",")
    if inhibitors is not None:
        inhibitors = inhibitors.split(",")
    if phosphosite_con is not None:
        phosphosite_con = int(phosphosite_con)

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

                    phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).all()

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
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.modification == "phosphorylation").all()
                        elif phosphosite_con ==  1:
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                    elif autophosphorylation == "exclude":
                        #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).all()
                        if not phosphosite_con or phosphosite_con == 0:
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                        elif phosphosite_con ==  1:
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                    
                    kin_list = []
                    phos_list = []
                    fc_list = []
                    cell_list = []
                    if weighted == "True":
                        confidence_list = []
                    for e_data, i_data in all_info:
                        #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                        kin = i_data.source
                        phosphosite = e_data.phosphosite

                        if weighted == "True":
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

                    if not weighted or weighted == "False":
                        s = kin_fc["fc"].mean()

                        sqr_m = numpy.sqrt(len(kin_fc))

                        #print(s)

                        #print(sqr_m)

                        z = ((s-p)*sqr_m)/standard_d_fc

                        #print(z)

                        p_value = norm.sf(z)

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
                    
                    elif weighted == "True":
                        # weighted calculations
                        weights = numpy.array(confidence_list)
                        fcs = numpy.array(fc_list)

                        weighted_mean = numpy.average(fcs, weights=weights)
                        weighted_sd = calculate_weighted_sd(fcs, weights)

                        z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                        p_value = norm.sf(z)
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
                        
                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.perturbation==inhibitor).all()

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
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.modification == "phosphorylation").all()
                            elif phosphosite_con ==  1:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                        elif autophosphorylation == "exclude":
                            #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).all()
                            if not phosphosite_con or phosphosite_con == 0:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                            elif phosphosite_con ==  1:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                        


                        #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).all()
                        kin_list = []
                        phos_list = []
                        fc_list = []
                        cell_list = []
                        if weighted == "True":
                            confidence_list = []
                        for e_data, i_data in all_info:
                            #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                            kin = i_data.source
                            phosphosite = e_data.phosphosite

                            if weighted == "True":
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
                        if not weighted or weighted == "False":
                            s = kin_fc["fc"].mean()

                            sqr_m = numpy.sqrt(len(kin_fc))

                            #print(s)

                            #print(sqr_m)

                            z = ((s-p)*sqr_m)/standard_d_fc

                            #print(z)

                            p_value = norm.sf(z)

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
                        
                        elif weighted == "True":
                            # weighted calculations
                            weights = numpy.array(confidence_list)
                            fcs = numpy.array(fc_list)

                            weighted_mean = numpy.average(fcs, weights=weights)
                            weighted_sd = calculate_weighted_sd(fcs, weights)

                            z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                            p_value = norm.sf(z)
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

                    phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.cell_line==cell_line).all()

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
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.modification == "phosphorylation").all()
                        elif phosphosite_con ==  1:
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                    elif autophosphorylation == "exclude":
                        if not phosphosite_con or phosphosite_con == 0:
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                        elif phosphosite_con ==  1:
                            all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                    



                    #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.cell_line==cell_line).all()
                    kin_list = []
                    phos_list = []
                    fc_list = []
                    cell_list = []
                    if weighted == "True":
                        confidence_list = []
                    for e_data, i_data in all_info:
                        #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                        kin = i_data.source
                        phosphosite = e_data.phosphosite

                        if weighted == "True":
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
                    if not weighted or weighted == "False":
                        s = kin_fc["fc"].mean()

                        sqr_m = numpy.sqrt(len(kin_fc))

                        #print(s)

                        #print(sqr_m)

                        z = ((s-p)*sqr_m)/standard_d_fc

                        #print(z)

                        p_value = norm.sf(z)

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
                    
                    elif weighted == "True":
                        # weighted calculations
                        weights = numpy.array(confidence_list)
                        fcs = numpy.array(fc_list)

                        weighted_mean = numpy.average(fcs, weights=weights)
                        weighted_sd = calculate_weighted_sd(fcs, weights)

                        z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                        p_value = norm.sf(z)
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
                for kinase in kinases:
                    kinase = kinase.upper()
                    # all_inhibitors = session.query(Experimental.perturbation).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).all()
                    # inhibitor_list = []

                    # for i in all_inhibitors:
                    #     inhibitor_list.append(i[0])
                    
                    for inhibitor in inhibitors:
                        #inter = session.query(Interaction.phosphosite_id).distinct().all()
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

                        phos_fc = session.query(Experimental.phosphosite, Experimental.fc).filter(Experimental.phosphosite.in_(inter_distinct)).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).all()

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
                        #     all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).all()
                        # elif autophosphorylation == "exclude":
                        #     all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).all()

                        if autophosphorylation == "include" or not autophosphorylation:
                            if not phosphosite_con or phosphosite_con == 0:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).filter(Interaction.modification == "phosphorylation").all()
                            elif phosphosite_con ==  1:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                        elif autophosphorylation == "exclude":
                            if not phosphosite_con or phosphosite_con == 0:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Interaction.modification == "phosphorylation").all()
                            elif phosphosite_con ==  1:
                                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).join(Phosphosite,Experimental.phosphosite == Phosphosite.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).filter(Interaction.source != Interaction.target).filter(Phosphosite.confidence == "High").filter(Interaction.modification == "phosphorylation").all()
                        


                        #all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).filter(Experimental.perturbation==inhibitor).filter(Experimental.cell_line==cell_line).all()
                        kin_list = []
                        phos_list = []
                        fc_list = []
                        cell_list = []
                        if weighted == "True":
                            confidence_list = []
                        for e_data, i_data in all_info:
                            #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                            kin = i_data.source
                            phosphosite = e_data.phosphosite

                            if weighted == "True":
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
                        print(kin_fc)
                        print(len(kin_fc))
                        if not weighted or weighted == "False":
                            s = kin_fc["fc"].mean()

                            sqr_m = numpy.sqrt(len(kin_fc))

                            #print(s)

                            #print(sqr_m)

                            z = ((s-p)*sqr_m)/standard_d_fc

                            #print(z)

                            p_value = norm.sf(z)

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
                        
                        elif weighted == "True":
                            # weighted calculations
                            weights = numpy.array(confidence_list)
                            fcs = numpy.array(fc_list)

                            weighted_mean = numpy.average(fcs, weights=weights)
                            weighted_sd = calculate_weighted_sd(fcs, weights)

                            z = ((weighted_mean - p) * numpy.sqrt(len(fcs))) / weighted_sd
                            p_value = norm.sf(z)
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
            
        
    finally:
        session.close()