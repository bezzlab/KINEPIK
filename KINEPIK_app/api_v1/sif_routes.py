from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
from KINEPIK_app.database_con import session_local
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# naming the blueprint
sif_bp = Blueprint("sif_bp", __name__, url_prefix="/api/sif")

@sif_bp.route("/info", methods = ["GET"])
def all_instructions():
    '''This function hands out the user info about the possible parameters that can be given to the kinases queries. 
    The return returns the info from string/text to html which is then displayed to the user'''
    # information about the parameters
    phosphosites_info = "phosphosites : 0 or 1"
    # the format of the html page
    info = f"""
    <html>
        <body>
            <p>{phosphosites_info}</p>
        </body>
    </html>
    """

    # return changes the format of the string to html that will be displayed to the user
    return Response(info, mimetype = "text/html")


@sif_bp.route("/all", methods = ["GET"])
def get_all_kinases():
    # connection to the database
    session = session_local()
    resolution = request.args.get("resolution")

    try:
        sif_list = []
        all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
        for kinase in all_kinase_info:
            inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
            effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(Effect.inhibit == 1).all()
            target_kin_list = []
            target_phos_list = []
            for inter in inter_info:
                target_kin_list.append(inter.target)
                target_phos_list.append(inter.phosphosite_id)
            
            for effect in effect_info:
                target_kin_list.append(effect.target)

            target_kin_list = set(target_kin_list)

            if resolution == "phosphosites":
                for target in target_phos_list:
                    sif_line = f"{kinase.id} phosphorylates {target}"
                    sif_list.append(sif_line)
                    source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                    sif_line = f"{kinase.id} gene {source_gene.gene}"
                    sif_list.append(sif_line)
                    target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                    sif_line = f"{target} on {target_kin.uniprot_id}"
                    sif_list.append(sif_line)
                    target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                    sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                    sif_list.append(sif_line)
            if resolution == "kinases":
                for target in target_kin_list:
                    sif_line = f"{kinase.id} phosphorylates {target}"
                    sif_list.append(sif_line)
                    source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                    sif_line = f"{kinase.id} gene {source_gene.gene}"
                    sif_list.append(sif_line)
                    target_gene = session.query(Protein).filter(Protein.id == target).first()
                    sif_line = f"{target} gene {target_gene.gene}"
                    sif_list.append(sif_line)

        sif_file = "\n".join(sif_list)

        return Response(sif_file, mimetype = "text")
        
    finally:
        session.close()


@sif_bp.route("/specific", methods = ["GET"])
def get_kinase():
    # connection to the database
    session = session_local()
    kinases = request.args.get("kinase_ids")
    kinases = kinases.split(",")
    resolution = request.args.get("resolution")

    try:
        sif_list = [] 
        for kinase in kinases:
            inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
            effect_info = session.query(Effect).filter(Effect.source == kinase).filter(Effect.inhibit == 1).all()
            target_kin_list = []
            target_phos_list = []
            for inter in inter_info:
                target_kin_list.append(inter.target)
                target_phos_list.append(inter.phosphosite_id)
            
            for effect in effect_info:
                target_kin_list.append(effect.target)

            target_kin_list = set(target_kin_list)

            if resolution == "phosphosites":
                for target in target_phos_list:
                    sif_line = f"{kinase} phosphorylates {target}"
                    sif_list.append(sif_line)
                    source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                    sif_line = f"{kinase} gene {source_gene.gene}"
                    sif_list.append(sif_line)
                    target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                    sif_line = f"{target} on {target_kin.uniprot_id}"
                    sif_list.append(sif_line)
                    target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                    sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                    sif_list.append(sif_line)
            if resolution == "kinases":
                for target in target_kin_list:
                    sif_line = f"{kinase} phosphorylates {target}"
                    sif_list.append(sif_line)
                    source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                    sif_line = f"{kinase} gene {source_gene.gene}"
                    sif_list.append(sif_line)
                    target_gene = session.query(Protein).filter(Protein.id == target).first()
                    sif_line = f"{target} gene {target_gene.gene}"
                    sif_list.append(sif_line)

        sif_file = "\n".join(sif_list)

        return Response(sif_file, mimetype = "text")
        
    finally:
        session.close()


@sif_bp.route("/attributes", methods = ["GET"])
def get_attributes():
    # connection to the database
    session = session_local()
    kinases = request.args.get("kinases")
    kinases = kinases.split(",")
    resolution = request.args.get("resolution")

    try:
        attribute_list = []
        attribute_list.append(f"ID Gene Uniprot Type")
        if kinases == "all": 
            all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
            for kinase in all_kinase_info:
                inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
                effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(Effect.inhibit == 1).all()
                target_kin_list = []
                target_phos_list = []
                for inter in inter_info:
                    target_kin_list.append(inter.target)
                    target_phos_list.append(inter.phosphosite_id)
                
                for effect in effect_info:
                    target_kin_list.append(effect.target)

                target_kin_list = set(target_kin_list)

                if resolution == "phosphosites":
                    line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
                    attribute_list.append(line)
                    for target in target_phos_list:
                        target_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").filter(Interaction.phosphosite_id==target).first()
                        gene = target.split("(")[0]
                        line = f"{target} {gene} {target_info.target} Phosphosite"
                        attribute_list.append(line)
                        # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                        # sif_line = f"{kinase} gene {source_gene.gene}"
                        # sif_list.append(sif_line)
                        # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                        # sif_line = f"{target} on {target_kin.uniprot_id}"
                        # sif_list.append(sif_line)
                        # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                        # sif_line = f"{kinase} gene {target_gene.gene}"
                        # sif_list.append(sif_line)
                if resolution == "kinases":
                    line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
                    for target in target_kin_list:
                        target_info = session.query(Protein).filter(Protein.id == target).first()
                        line = f"{target} {target_info.gene} {target} Kinase"
                        attribute_list.append(line)
                        # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                        # sif_line = f"{kinase} gene {source_gene.gene}"
                        # sif_list.append(sif_line)
                        # target_gene = session.query(Protein).filter(Protein.id == target).first()
                        # sif_line = f"{target} gene {target_gene.gene}"
                        # sif_list.append(sif_line)

            sif_file = "\n".join(attribute_list)

            return Response(sif_file, mimetype = "text")
        
        elif kinases != "all":
            for kinase in kinases:
                inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
                effect_info = session.query(Effect).filter(Effect.source == kinase).filter(Effect.inhibit == 1).all()
                target_kin_list = []
                target_phos_list = []
                for inter in inter_info:
                    target_kin_list.append(inter.target)
                    target_phos_list.append(inter.phosphosite_id)
                
                for effect in effect_info:
                    target_kin_list.append(effect.target)

                target_kin_list = set(target_kin_list)

                if resolution == "phosphosites":
                    for target in target_phos_list:
                        target_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").filter(Interaction.phosphosite_id==target).first()
                        gene = target.split("(")[0]
                        line = f"{target} {gene} {target_info.target} Phosphosite"
                        attribute_list.append(line)
                        # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                        # sif_line = f"{kinase} gene {source_gene.gene}"
                        # sif_list.append(sif_line)
                        # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                        # sif_line = f"{target} on {target_kin.uniprot_id}"
                        # sif_list.append(sif_line)
                        # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                        # sif_line = f"{kinase} gene {target_gene.gene}"
                        # sif_list.append(sif_line)
                if resolution == "kinases":
                    for target in target_kin_list:
                        target_info = session.query(Protein).filter(Protein.id == target).first()
                        line = f"{target} {target_info.gene} {target} Kinase"
                        attribute_list.append(line)
                        # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                        # sif_line = f"{kinase} gene {source_gene.gene}"
                        # sif_list.append(sif_line)
                        # target_gene = session.query(Protein).filter(Protein.id == target).first()
                        # sif_line = f"{target} gene {target_gene.gene}"
                        # sif_list.append(sif_line)

            sif_file = "\n".join(attribute_list)

        return Response(sif_file, mimetype = "text")
        
    finally:
        session.close()

# @kinase_bp.route("/specific", methods = ["GET"])
# def get_kinase():
#     # connection to the database
#     session = session_local()
#     kinases_json = []
#     kinase_ids = request.args.get("kinase_ids")
#     kinase_ids = kinase_ids.split(",")
#     phosphosites = request.args.get("phosphosites")
#     phos_confidence = request.args.get("confidence")
#     if phos_confidence is not None:
#         phos_confidence = int(phos_confidence)

#     try:
#         if kinase_ids == None:
#             return "Kinase id required for results"
#         else:
#             if phosphosites == None and (phos_confidence == 0 or phos_confidence == None):
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id=kin).first()
#                     kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
#                     target_phosphosites_info = session.query(Interaction).filter(Interaction.source == all_kinase_info.id).filter(Interaction.target != all_kinase_info.id).all()
#                     kinase_phos_list = []
#                     target_phos_list = []
#                     for k_phos in kinase_phosphosites_info:
#                         kinase_phos_list.append(k_phos.phosphosite_id)
#                     for t_phos in target_phosphosites_info:
#                         target_phos_list.append(t_phos.phosphosite_id)
#                     kinase_info = [{
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "TargetPhosphosites": target_phos_list,
#                             "PhosphositesOnKinase": kinase_phos_list,
#                             "OrganismID" : all_kinase_info.species                        }]
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)

#             elif phosphosites == None and phos_confidence == 1:
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id=kin).first()
#                     kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).filter_by(confidence = "High").all()
#                     target_phosphosites_info = session.query(Interaction,Phosphosite).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source == all_kinase_info.id).filter(Phosphosite.confidence == "High").all()
#                     #target_phosphosites_info = session.query(Interaction).filter_by(source = all_kinase_info.id).filter_by(confidence = "high").all()
#                     kinase_phos_list = []
#                     target_phos_list = []
#                     for k_phos in kinase_phosphosites_info:
#                         kinase_phos_list.append(k_phos.phosphosite_id)
#                     for inter_id,t_phos in target_phosphosites_info:
#                         target_phos_list.append(t_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "TargetPhosphosites": target_phos_list,
#                             "PhosphositesOnKinase": kinase_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
            
#             elif phosphosites == "targets" and (phos_confidence == 0 or phos_confidence == None):
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id = kin).first()
#                     target_phosphosites_info = session.query(Interaction).filter(Interaction.source == all_kinase_info.id).filter(Interaction.target != all_kinase_info.id).all()
#                     target_phos_list = []
#                     for t_phos in target_phosphosites_info:
#                         target_phos_list.append(t_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "TargetPhosphosites": target_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
            
#             elif phosphosites == "targets" and phos_confidence == 1:
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id = kin).first()
#                     target_phosphosites_info = session.query(Interaction,Phosphosite).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source == all_kinase_info.id).filter(Phosphosite.confidence == "High").all()
#                     target_phos_list = []
#                     for inter_id,t_phos in target_phosphosites_info:
#                         target_phos_list.append(t_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "TargetPhosphosites": target_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
            
            
#             elif phosphosites == "sites" and (phos_confidence == 0 or phos_confidence == None):
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id = kin).first()
#                     kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
#                     kinase_phos_list = []
#                     for k_phos in kinase_phosphosites_info:
#                         kinase_phos_list.append(k_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "PhosphositesOnKinase": kinase_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
            
#             elif phosphosites == "sites" and phos_confidence == 1:
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id = kin).first()
#                     kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).filter_by(confidence = "High").all()
#                     kinase_phos_list = []
#                     for k_phos in kinase_phosphosites_info:
#                         kinase_phos_list.append(k_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "PhosphositesOnKinase": kinase_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
            
#             elif ("sites" in phosphosites and "targets" in phosphosites) and (phos_confidence == 0 or phos_confidence == None):
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id = kin).first()
#                     kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).all()
#                     target_phosphosites_info = session.query(Interaction).filter(Interaction.source == all_kinase_info.id).filter(Interaction.target != all_kinase_info.id).all()
#                     kinase_phos_list = []
#                     target_phos_list = []
#                     for k_phos in kinase_phosphosites_info:
#                         kinase_phos_list.append(k_phos.phosphosite_id)
#                     for t_phos in target_phosphosites_info:
#                         target_phos_list.append(t_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "PhosphositesOnKinase": kinase_phos_list,
#                             "TargetPhosphosites": target_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
            
#             elif ("sites" in phosphosites and "targets" in phosphosites) and phos_confidence == 1:
#                 for kin in kinase_ids:
#                     kin = kin.upper()
#                     all_kinase_info = session.query(Protein).filter_by(id=kin).first()
#                     kinase_phosphosites_info = session.query(Phosphosite).filter_by(uniprot_id = all_kinase_info.id).filter_by(confidence = "High").all()
#                     target_phosphosites_info = session.query(Interaction,Phosphosite).join(Phosphosite, Interaction.phosphosite_id == Phosphosite.phosphosite_id).filter(Interaction.source == all_kinase_info.id).filter(Phosphosite.confidence == "High").all()
#                     #target_phosphosites_info = session.query(Interaction).filter_by(source = all_kinase_info.id).filter_by(confidence = "high").all()
#                     kinase_phos_list = []
#                     target_phos_list = []
#                     for k_phos in kinase_phosphosites_info:
#                         kinase_phos_list.append(k_phos.phosphosite_id)
#                     for inter_id,t_phos in target_phosphosites_info:
#                         target_phos_list.append(t_phos.phosphosite_id)
#                     kinase_info = {
#                             "SourceUniprotID" : all_kinase_info.id,
#                             "UniprotName" : all_kinase_info.name,
#                             "TargetPhosphosites": target_phos_list,
#                             "PhosphositesOnKinase": kinase_phos_list,
#                             "OrganismID" : all_kinase_info.species                        
#                             }
#                     kinases_json.append(kinase_info)
#                 return jsonify(kinases_json)
        
#     finally:
#         session.close()