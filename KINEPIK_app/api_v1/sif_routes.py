from flask import Blueprint,jsonify, json, request, Response
import pandas as pd
from KINEPIK_app.database_con import session_local
from sqlalchemy import or_
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# naming the blueprint
sif_bp = Blueprint("sif_bp", __name__, url_prefix="/api/sif")

@sif_bp.route("/info/<route>", methods = ["GET"])
def instructions(route):
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
def get_sif_all_kinases():
    # connection to the database
    session = session_local()
    resolution = request.args.get("resolution")

    try:
        sif_list = []
        all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
        for kinase in all_kinase_info:
            inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
            effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(or_(Effect.inhibit == 1,Effect.stimulate == 1)).all()
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
                    # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                    # sif_line = f"{kinase.id} gene {source_gene.gene}"
                    # sif_list.append(sif_line)
                    target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                    sif_line = f"{target} on {target_kin.uniprot_id}"
                    sif_list.append(sif_line)
                    # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                    # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                    # sif_list.append(sif_line)
            if resolution == "kinases":
                for target in target_kin_list:
                    sif_line = f"{kinase.id} phosphorylates {target}"
                    sif_list.append(sif_line)
                    # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                    # sif_line = f"{kinase.id} gene {source_gene.gene}"
                    # sif_list.append(sif_line)
                    # target_gene = session.query(Protein).filter(Protein.id == target).first()
                    # sif_line = f"{target} gene {target_gene.gene}"
                    # sif_list.append(sif_line)

        sif_file = "\n".join(sif_list)

        return Response(sif_file, mimetype = "text")
        
    finally:
        session.close()


@sif_bp.route("/specific", methods = ["GET"])
def get_sif_kinase():
    # connection to the database
    session = session_local()
    kinases = request.args.get("kinase_ids")
    kinases = kinases.split(",")
    resolution = request.args.get("resolution")

    try:
        sif_list = [] 
        for kinase in kinases:
            inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
            effect_info = session.query(Effect).filter(Effect.source == kinase).filter(or_(Effect.inhibit == 1,Effect.stimulate == 1)).all()
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
                    # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                    # sif_line = f"{kinase} gene {source_gene.gene}"
                    # sif_list.append(sif_line)
                    target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                    sif_line = f"{target} on {target_kin.uniprot_id}"
                    sif_list.append(sif_line)
                    # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                    # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                    # sif_list.append(sif_line)
            if resolution == "kinases":
                for target in target_kin_list:
                    sif_line = f"{kinase} phosphorylates {target}"
                    sif_list.append(sif_line)
                    # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
                    # sif_line = f"{kinase} gene {source_gene.gene}"
                    # sif_list.append(sif_line)
                    # target_gene = session.query(Protein).filter(Protein.id == target).first()
                    # sif_line = f"{target} gene {target_gene.gene}"
                    # sif_list.append(sif_line)

        sif_file = "\n".join(sif_list)

        return Response(sif_file, mimetype = "text")
        
    finally:
        session.close()

@sif_bp.route("/attributes", methods = ["GET"])
def get_attributes():
    # connection to the database
    session = session_local()
    kinases = request.args.get("kinases")
    resolution = request.args.get("resolution")
    type = request.args.get("type")



    try:
        attribute_list = []
        # if type == "IDs":
        #     attribute_list.append(f"ID Site_id")
        # elif type == "type":
        #     attribute_list.append(f"ID Type")
        #attribute_list.append(f"ID Gene Uniprot Type")
        if kinases == "all": 
            all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
            for kinase in all_kinase_info:
                inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
                effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(or_(Effect.inhibit == 1,Effect.stimulate == 1)).all()
                target_kin_list = []
                target_phos_list = []
                for inter in inter_info:
                    target_kin_list.append(inter.target)
                    target_phos_list.append(inter.phosphosite_id)
                
                for effect in effect_info:
                    target_kin_list.append(effect.target)

                target_kin_list = set(target_kin_list)

                if resolution == "phosphosites":
                    if type == "IDs":
                        for target in target_phos_list:
                            #sif_line = f"{kinase.id} phosphorylates {target}"
                            sif_line = f"{kinase.id} {kinase.gene}"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                            #sif_line = f"{target} on {target_kin.uniprot_id}"
                            sif_line = f"{target} {target}"
                            attribute_list.append(sif_line)
                            target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                            # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                            sif_line = f"{target_kin.uniprot_id} {target_gene.gene}"
                            attribute_list.append(sif_line)
                    elif type == "type":
                        for target in target_phos_list:
                            sif_line = f"{kinase.id} Kinase"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                            sif_line = f"{target} Phosphosite"
                            attribute_list.append(sif_line)
                            target_data = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                            if target_data.kinase == 0:
                                target_type = "Protein"
                            elif target_data.kinase == 1:
                                target_type = "Kinase"
                            sif_line = f"{target_kin.uniprot_id} {target_type}"
                            attribute_list.append(sif_line)
                if resolution == "kinases":
                    if type == "IDs":
                        for target in target_phos_list:
                            #sif_line = f"{kinase.id} phosphorylates {target}"
                            sif_line = f"{kinase.id} {kinase.gene}"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            target_kin = session.query(Protein).filter(Protein.id == target).first()
                            #sif_line = f"{target} on {target_kin.uniprot_id}"
                            sif_line = f"{target} {target_kin.gene}"
                            attribute_list.append(sif_line)
                            # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                            # # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                            # sif_line = f"{target_kin.uniprot_id} {target_gene.gene}"
                            # attribute_list.append(sif_line)
                    elif type == "type":
                        for target in target_phos_list:
                            sif_line = f"{kinase.id} Kinase"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                            # sif_line = f"{target} Phosphosite"
                            # attribute_list.append(sif_line)
                            target_data = session.query(Protein).filter(Protein.id == target).first()
                            if target_data.kinase == 0:
                                target_type = "Protein"
                            elif target_data.kinase == 1:
                                target_type = "Kinase"
                            sif_line = f"{target} {target_type}"
                            attribute_list.append(sif_line)
            attribute_list = set(attribute_list)
            sif_file = "\n".join(attribute_list)

            return Response(sif_file, mimetype = "text")
        
        elif kinases != "all":
            kinases = kinases.split(",")
            for kinase in kinases:
                inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
                if resolution == "kinases":
                    effect_info = session.query(Effect).filter(Effect.source == kinase).filter(or_(Effect.inhibit == 1,Effect.stimulate==1)).all()
                target_kin_list = []
                target_phos_list = []
                for inter in inter_info:
                    target_kin_list.append(inter.target)
                    target_phos_list.append(inter.phosphosite_id)
                
                if resolution == "kinases":
                    for effect in effect_info:
                        target_kin_list.append(effect.target)

                target_kin_list = set(target_kin_list)

                if resolution == "phosphosites":
                    if type == "IDs":
                        for target in target_phos_list:
                            #sif_line = f"{kinase.id} phosphorylates {target}"
                            source = session.query(Protein).filter(Protein.id == kinase).first()
                            sif_line = f"{kinase} {source.gene}"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                            #sif_line = f"{target} on {target_kin.uniprot_id}"
                            sif_line = f"{target} {target}"
                            attribute_list.append(sif_line)
                            target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                            # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                            sif_line = f"{target_kin.uniprot_id} {target_gene.gene}"
                            attribute_list.append(sif_line)
                    elif type == "type":
                        for target in target_phos_list:
                            sif_line = f"{kinase} Kinase"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                            sif_line = f"{target} Phosphosite"
                            attribute_list.append(sif_line)
                            target_data = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                            if target_data.kinase == 0:
                                target_type = "Protein"
                            elif target_data.kinase == 1:
                                target_type = "Kinase"
                            sif_line = f"{target_kin.uniprot_id} {target_type}"
                            attribute_list.append(sif_line)
                    elif type == "known":
                        sif_line = f"{kinase} Kinase"
                        attribute_list.append(sif_line)
                
                
                if resolution == "kinases":
                    if type == "IDs":
                        for target in target_phos_list:
                            #sif_line = f"{kinase.id} phosphorylates {target}"
                            source = session.query(Protein).filter(Protein.id == kinase).first()
                            sif_line = f"{kinase} {source.gene}"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            target_kin = session.query(Protein).filter(Protein.id == target).first()
                            #sif_line = f"{target} on {target_kin.uniprot_id}"
                            sif_line = f"{target} {target_kin.gene}"
                            attribute_list.append(sif_line)
                            # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
                            # # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
                            # sif_line = f"{target_kin.uniprot_id} {target_gene.gene}"
                            # attribute_list.append(sif_line)
                    elif type == "type":
                        for target in target_phos_list:
                            sif_line = f"{kinase} Kinase"
                            attribute_list.append(sif_line)
                            # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
                            # sif_line = f"{kinase.id} gene {source_gene.gene}"
                            # attribute_list.append(sif_line)
                            # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
                            # sif_line = f"{target} Phosphosite"
                            # attribute_list.append(sif_line)
                            target_data = session.query(Protein).filter(Protein.id == target).first()
                            if target_data.kinase == 0:
                                target_type = "Protein"
                            elif target_data.kinase == 1:
                                target_type = "Kinase"
                            sif_line = f"{target} {target_type}"
                            attribute_list.append(sif_line)
            attribute_list = set(attribute_list)
            sif_file = "\n".join(attribute_list)

        return Response(sif_file, mimetype = "text")
        
    finally:
        session.close()

# @sif_bp.route("/attributes", methods = ["GET"])
# def get_attributes():
#     # connection to the database
#     session = session_local()
#     kinases = request.args.get("kinases")
#     kinases = kinases.split(",")
#     resolution = request.args.get("resolution")

#     try:
#         attribute_list = []
#         attribute_list.append(f"ID Gene Uniprot Type")
#         if kinases == "all": 
#             all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
#             for kinase in all_kinase_info:
#                 inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
#                 effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(or_(Effect.inhibit == 1,Effect.stimulate == 1)).all()
#                 target_kin_list = []
#                 target_phos_list = []
#                 for inter in inter_info:
#                     target_kin_list.append(inter.target)
#                     target_phos_list.append(inter.phosphosite_id)
                
#                 for effect in effect_info:
#                     target_kin_list.append(effect.target)

#                 target_kin_list = set(target_kin_list)

#                 if resolution == "phosphosites":
#                     line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                     attribute_list.append(line)
#                     for target in target_phos_list:
#                         target_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").filter(Interaction.phosphosite_id==target).first()
#                         gene = target.split("(")[0]
#                         line = f"{target} {gene} {target_info.target} Phosphosite"
#                         attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
#                         # sif_line = f"{target} on {target_kin.uniprot_id}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
#                         # sif_line = f"{kinase} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)
#                 if resolution == "kinases":
#                     line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                     for target in target_kin_list:
#                         target_info = session.query(Protein).filter(Protein.id == target).first()
#                         line = f"{target} {target_info.gene} {target} Kinase"
#                         attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target).first()
#                         # sif_line = f"{target} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)

#             sif_file = "\n".join(attribute_list)

#             return Response(sif_file, mimetype = "text")
        
#         elif kinases != "all":
#             for kinase in kinases:
#                 inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
#                 effect_info = session.query(Effect).filter(Effect.source == kinase).filter(or_(Effect.inhibit == 1,Effect.stimulate == 1)).all()
#                 target_kin_list = []
#                 target_phos_list = []
#                 for inter in inter_info:
#                     target_kin_list.append(inter.target)
#                     target_phos_list.append(inter.phosphosite_id)
                
#                 for effect in effect_info:
#                     target_kin_list.append(effect.target)

#                 target_kin_list = set(target_kin_list)

#                 if resolution == "phosphosites":
#                     for target in target_phos_list:
#                         target_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").filter(Interaction.phosphosite_id==target).first()
#                         gene = target.split("(")[0]
#                         line = f"{target} {gene} {target_info.target} Phosphosite"
#                         attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
#                         # sif_line = f"{target} on {target_kin.uniprot_id}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
#                         # sif_line = f"{kinase} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)
#                 if resolution == "kinases":
#                     for target in target_kin_list:
#                         target_info = session.query(Protein).filter(Protein.id == target).first()
#                         line = f"{target} {target_info.gene} {target} Kinase"
#                         attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target).first()
#                         # sif_line = f"{target} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)

#             sif_file = "\n".join(attribute_list)

#         return Response(sif_file, mimetype = "text")
        
#     finally:
#         session.close()








# from flask import Blueprint,jsonify, json, request, Response
# import pandas as pd
# from sqlalchemy import or_
# from KINEPIK_app.database_con import session_local
# from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Tissue_location, Experimental, Protein, Gene

# # naming the blueprint
# sif_bp = Blueprint("sif_bp", __name__, url_prefix="/api/sif")

# @sif_bp.route("/info", methods = ["GET"])
# def all_instructions():
#     '''This function hands out the user info about the possible parameters that can be given to the kinases queries. 
#     The return returns the info from string/text to html which is then displayed to the user'''
#     # information about the parameters
#     phosphosites_info = "phosphosites : 0 or 1"
#     # the format of the html page
#     info = f"""
#     <html>
#         <body>
#             <p>{phosphosites_info}</p>
#         </body>
#     </html>
#     """

#     # return changes the format of the string to html that will be displayed to the user
#     return Response(info, mimetype = "text/html")


# @sif_bp.route("/all", methods = ["GET"])
# def get_all_kinases():
#     # connection to the database
#     session = session_local()
#     resolution = request.args.get("resolution")

#     try:
#         sif_list = []
#         all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
#         for kinase in all_kinase_info:
#             inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
#             if resolution == "kinases":
#                 effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(or_(Effect.inhibit == 1,Effect.stimulate==1)).all()
#             target_kin_list = []
#             target_phos_list = []
#             for inter in inter_info:
#                 target_kin_list.append(inter.target)
#                 target_phos_list.append(inter.phosphosite_id)

#             if resolution == "kinases":
#                 for effect in effect_info:
#                     target_kin_list.append(effect.target)

#             #target_kin_list = set(target_kin_list)

#             if resolution == "phosphosites":
#                 for target in target_phos_list:
#                     sif_line = f"{kinase.id} Kinase phosphorylates {target} Phosphosite"
#                     sif_list.append(sif_line)
#                     # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
#                     # sif_line = f"{kinase.id} gene {source_gene.gene}"
#                     # sif_list.append(sif_line)
#                     # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
#                     # sif_line = f"{target} on {target_kin.uniprot_id}"
#                     # sif_list.append(sif_line)
#                     # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
#                     # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
#                     # sif_list.append(sif_line)
#             if resolution == "kinases":
#                 for target in target_kin_list:
#                     target_info = session.query(Protein).filter_by(id = target).all()
#                     if target_info.kinase == 0:
#                         target_type = "Protein"
#                     elif target_info.kinase == 1:
#                         target_type = "Kinase"
#                     sif_line = f"{kinase.id} Kinase phosphorylates {target} {target_type}"
#                     sif_list.append(sif_line)
#                     # source_gene = session.query(Protein).filter(Protein.id == kinase.id).first()
#                     # sif_line = f"{kinase.id} gene {source_gene.gene}"
#                     # sif_list.append(sif_line)
#                     # target_gene = session.query(Protein).filter(Protein.id == target).first()
#                     # sif_line = f"{target} gene {target_gene.gene}"
#                     # sif_list.append(sif_line)

#         sif_file = "\n".join(sif_list)

#         return Response(sif_file, mimetype = "text")
        
#     finally:
#         session.close()


# @sif_bp.route("/specific", methods = ["GET"])
# def get_kinase():
#     # connection to the database
#     session = session_local()
#     kinases = request.args.get("kinase_ids")
#     kinases = kinases.split(",")
#     resolution = request.args.get("resolution")

#     try:
#         sif_list = [] 
#         for kinase in kinases:
#             inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
#             if resolution == "kinases":
#                 effect_info = session.query(Effect).filter(Effect.source == kinase).filter(or_(Effect.inhibit == 1,Effect.stimulate==1)).all()
#             target_kin_list = []
#             target_phos_list = []
#             for inter in inter_info:
#                 target_kin_list.append(inter.target)
#                 target_phos_list.append(inter.phosphosite_id)
            
#             if resolution == "kinases":
#                 for effect in effect_info:
#                     target_kin_list.append(effect.target)

#             target_kin_list = set(target_kin_list)

#             if resolution == "phosphosites":
#                 for target in target_phos_list:
#                     kin_info = session.query(Protein).filter(Protein.id == kinase).first()
#                     target_info = session.query(Interaction).filter(Interaction.source == kinase).filter(Interaction.phosphosite_id == target).first()
#                     target_gene = session.query(Protein).filter(Protein.id == target_info.target).first()
#                     if target_gene.kinase == 0:
#                         type = "Protein"
#                     elif target_gene.kinase == 1:
#                         type = "Kinase"
#                     sif_line = f"{kinase} Kinase {kin_info.gene} phosphorylates {target} Phosphosite on {target_info.target} {type} {target_gene.gene}"
#                     sif_list.append(sif_line)
#                     # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                     # sif_line = f"{kinase} gene {source_gene.gene}"
#                     # sif_list.append(sif_line)
#                     # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
#                     # sif_line = f"{target} on {target_kin.uniprot_id}"
#                     # sif_list.append(sif_line)
#                     # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
#                     # sif_line = f"{target_kin.uniprot_id} gene {target_gene.gene}"
#                     # sif_list.append(sif_line)
#             if resolution == "kinases":
#                 for target in target_kin_list:
#                     target_info = session.query(Protein).filter_by(id = target).first()
#                     if target_info.kinase == 0:
#                         target_type = "Protein"
#                     elif target_info.kinase == 1:
#                         target_type = "Kinase"
#                     sif_line = f"{kinase} Kinase phosphorylates {target} {target_type}"
#                     sif_list.append(sif_line)
#                     # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                     # sif_line = f"{kinase} gene {source_gene.gene}"
#                     # sif_list.append(sif_line)
#                     # target_gene = session.query(Protein).filter(Protein.id == target).first()
#                     # sif_line = f"{target} gene {target_gene.gene}"
#                     # sif_list.append(sif_line)

#         sif_file = "\n".join(sif_list)

#         return Response(sif_file, mimetype = "text")
        
#     finally:
#         session.close()


# @sif_bp.route("/attributes", methods = ["GET"])
# def get_attributes():
#     # connection to the database
#     session = session_local()
#     kinases = request.args.get("kinases")
#     resolution = request.args.get("resolution")
#     type = request.args.get("type")

#     try:
#         attribute_list = []
#         if type == "gene":
#             attribute_list.append(f"ID Gene")
#         elif type == "on_kinase":
#             attribute_list.append(f"ID Uniprot")
#         #attribute_list.append(f"ID Gene Uniprot Type")
#         if kinases == "all": 
#             all_kinase_info = session.query(Protein).filter_by(kinase = 1).all()
#             for kinase in all_kinase_info:
#                 inter_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").all()
#                 if resolution == "kinases":
#                     effect_info = session.query(Effect).filter(Effect.source == kinase.id).filter(or_(Effect.inhibit == 1,Effect.stimulate==1)).all()
#                 target_kin_list = []
#                 target_phos_list = []
#                 for inter in inter_info:
#                     target_kin_list.append(inter.target)
#                     target_phos_list.append(inter.phosphosite_id)
                
#                 if resolution == "kinases":
#                     for effect in effect_info:
#                         target_kin_list.append(effect.target)

#                 target_kin_list = set(target_kin_list)

#                 if resolution == "phosphosites":
#                     if type == "gene":
#                         #line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                         #attribute_list.append(line)
#                         for target in target_phos_list:
#                             kin_info = session.query(Protein).filter_by(id = kinase.id).first()
#                             #gene = target.split("(")[0]
#                             line = f"{kinase.id} {kin_info.gene}"
#                             if line not in attribute_list:
#                                 attribute_list.append(line)
#                     if type == "on_kinase":
#                         #line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                         #attribute_list.append(line)
#                         for target in target_phos_list:
#                             target_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").filter(Interaction.phosphosite_id==target).first()
#                             #gene = target.split("(")[0]
#                             line = f"{target} {target_info.target}"
#                             if line not in attribute_list:
#                                 attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
#                         # sif_line = f"{target} on {target_kin.uniprot_id}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
#                         # sif_line = f"{kinase} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)
#                 if resolution == "kinases":
#                     if type == "gene":
#                     #line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                         for target in target_kin_list:
#                             target_info = session.query(Protein).filter(Protein.id == target).first()
#                             line = f"{target} {target_info.gene}"
#                             attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target).first()
#                         # sif_line = f"{target} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)

#             sif_file = "\n".join(attribute_list)

#             return Response(sif_file, mimetype = "text")
        
#         elif kinases != "all":
#             kinases = kinases.split(",")
#             if type == "gene":
#                 attribute_list.append(f"ID Gene")
#             elif type == "on_kinase":
#                 attribute_list.append(f"ID Uniprot")
#             for kinase in kinases:
#                 inter_info = session.query(Interaction).filter_by(source = kinase).filter(Interaction.modification == "phosphorylation").all()
#                 if resolution == "kinases":
#                     effect_info = session.query(Effect).filter(Effect.source == kinase).filter(or_(Effect.inhibit == 1,Effect.stimulate==1)).all()
#                 target_kin_list = []
#                 target_phos_list = []
#                 for inter in inter_info:
#                     target_kin_list.append(inter.target)
#                     target_phos_list.append(inter.phosphosite_id)
                
#                 if resolution == "kinases":
#                     for effect in effect_info:
#                         target_kin_list.append(effect.target)

#                 target_kin_list = set(target_kin_list)

#                 if resolution == "phosphosites":
#                     if type == "gene":
#                         #line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                         #attribute_list.append(line)
#                         for target in target_phos_list:
#                             kin_info = session.query(Protein).filter_by(id = kinase).all()
#                             #gene = target.split("(")[0]
#                             line = f"{kinase} {kin_info.gene}"
#                             attribute_list.append(line)
#                     if type == "on_kinase":
#                         #line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                         attribute_list.append(line)
#                         for target in target_phos_list:
#                             target_info = session.query(Interaction).filter_by(source = kinase.id).filter(Interaction.modification == "phosphorylation").filter(Interaction.phosphosite_id==target).first()
#                             #gene = target.split("(")[0]
#                             line = f"{target} {target_info.target}"
#                             attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_kin = session.query(Phosphosite).filter(Phosphosite.phosphosite_id == target).first()
#                         # sif_line = f"{target} on {target_kin.uniprot_id}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target_kin.uniprot_id).first()
#                         # sif_line = f"{kinase} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)
#                 if resolution == "kinases":
#                     if type == "gene":
#                     #line = f"{kinase.id} {kinase.gene} {kinase.id} Kinase"
#                         for target in target_kin_list:
#                             target_info = session.query(Protein).filter(Protein.id == target).first()
#                             line = f"{target} {target_info.gene}"
#                             attribute_list.append(line)
#                         # source_gene = session.query(Protein).filter(Protein.id == kinase).first()
#                         # sif_line = f"{kinase} gene {source_gene.gene}"
#                         # sif_list.append(sif_line)
#                         # target_gene = session.query(Protein).filter(Protein.id == target).first()
#                         # sif_line = f"{target} gene {target_gene.gene}"
#                         # sif_list.append(sif_line)

#             sif_file = "\n".join(attribute_list)

#         return Response(sif_file, mimetype = "text")
        
#     finally:
#         session.close()
