from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
from KINEPIK_app.database_con import session_local
import ast
from scipy.stats import norm
import numpy
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


@inhibitor_bp.route("/fc", methods = ["GET"])
def get_fc():
    # connection to the database
    session = session_local()
    fc_json = []
    kinases = request.args.get("kinase_id")
    if kinases is not None:
        kinases = kinases.split(",")

    try:
        if kinases == None:
            return "Uniprot id required for results"
        else:
            for kinase in kinases:
                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).all()
                for e_data, i_data in all_info:
                    kinase_num = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).all()
                    kinase_ref = session.query(Interaction).filter(Interaction.phosphosite_id == e_data.phosphosite).filter(Interaction.source == i_data.source).first()
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
                    confidence = ref_score * uniqueness_penalty
                    # uniqueness less penalised?
                    # confidence = (0.8 * ref_score + 0.2 * uniqueness_penalty) / (0.8 + 0.2)

                    fc_info = [{
                        i_data.source : {
                            "Inhibitor" : e_data.perturbation,
                            "Phosphosite" : e_data.phosphosite,
                            "FC" : e_data.fc,
                            "KinaseConfidence" : confidence
                            }
                    }]
                    fc_json.append(fc_info)
            
            return jsonify(fc_json)
        
    finally:
        session.close()


@inhibitor_bp.route("/KSEA", methods = ["GET"])
def get_ksea():
    session = session_local()
    ksea_json = []
    kinases = request.args.get("kinase_id")
    if kinases is not None:
        kinases = kinases.split(",")

    try:
        if kinases == None:
            return "Uniprot id required for results"
        else:
            for kinase in kinases:
                inter = session.query(Interaction.phosphosite_id).distinct().all()
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

                all_info = session.query(Experimental, Interaction).join(Interaction,Experimental.phosphosite == Interaction.phosphosite_id).filter(Interaction.source == kinase).all()
                kin_list = []
                phos_list = []
                fc_list = []
                for e_data, i_data in all_info:
                    #print("Kinase : " + i_data.source + ", phoshosite : " + e_data.phosphosite + ", fc : " + str(e_data.fc))
                    kin = i_data.source
                    phosphosite = e_data.phosphosite
                    fc = float(e_data.fc)
                    kin_list.append(kin)
                    phos_list.append(phosphosite)
                    fc_list.append(fc)

                kin_fc = pd.DataFrame({"kinase":kin_list,"phosphosite":phos_list,"fc":fc_list})
                #print(kin_fc)
                #print(len(kin_fc))
                s = kin_fc["fc"].mean()

                sqr_m = numpy.sqrt(len(kin_fc))

                #print(s)

                #print(sqr_m)

                z = ((s-p)*sqr_m)/standard_d_fc

                #print(z)

                p_value = norm.sf(abs(z))

                uniques_phos_list = set(phos_list)
                #print(p_value)
                ksea_info = [{
                    kinase : {
                        "Phosphosites" : list(uniques_phos_list),
                        "MeanFCKinase" : s,
                        "MeanFCAll" : p,
                        "StandardDeviation" : standard_d_fc,
                        "n" : len(kin_fc),
                        "z_score" : z,
                        "p_value" : p_value
                        }
                    }]
                ksea_json.append(ksea_info)
            
            return jsonify(ksea_json)
        
    finally:
        session.close()