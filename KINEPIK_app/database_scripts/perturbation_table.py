import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy import select, distinct
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Base,Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Subcell_location, Experimental, Protein, Gene

# creating the connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

# use these lines if you want to drop the old tables and recreate them again
# Perturbation.__table__.drop(engine)
# Base.metadata.create_all(engine)

def collectPerturbations():
    '''Function collects perturbations from the Perturbation_interaction and Experimental sql tables and returns a list with unique 
    perturbations'''
    # selecting all the perturbations from the tables that have perturbations
    inter_perts = session.query(Perturbation_interaction.perturbation).distinct().all()
    exp_perts = session.query(Experimental.perturbation).distinct().all()
    
    # combining all of them into one list to be handled further
    all_perts = inter_perts + exp_perts

    pert_list = []
    #org_pert_list = []
    # go through all the perturbations in the list
    for pert in all_perts:
        if pert[0] != None:
            # name = pert[0].lower().strip()
            name = pert[0].strip()
        # if the protein is not already in the pert_list, it will be added to it
            if name not in pert_list:
                pert_list.append(name)
                #org_pert_list.append(name)
    #pert_list = list(set(org_pert_list))
    # making sure there aren't any duplicates
    pert_list = list(set(pert_list))
    return pert_list 

# save protein list from the function to the variable
perturbations = collectPerturbations()
#print(perturbations)
#print(len(perturbations))

# change the list to pubchem cids
def nameToDF(pert_list):
    '''The function takes the list of the common names of the perturbations and creates a new list with their CIDs. The CIDs are collected
    from pubchem API endpoint by using the common names. Then with the CIDs other info about the perturbations are collected from pubhem.
    The function return a df with all the needed info for the sql table'''
    # creating the empty column lists that will be populated
    name_col = []
    type_col = []
    pubchem_col = []
    smiles_col = []
    synonyms_col = []
    gene_col = []
    action_col = []

    # going through the perturbations name by name
    for name in pert_list:
        #print(name[0])
        common_name = name #.capitalize()
        #print(common_name)
        # creating the url to pubchem API endpoint that will give the CID for the common name
        id_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+common_name+"/cids/TXT"
        # to make it look like a normal internet search
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
        #id_file = "id_data.txt"

        # requesting the infomation from the url
        response = requests.get(id_url, headers=headers)

        # if the connection works
        action = None # the null response that will be changed if results are found
        if response.status_code == 200 and "Status: 404" not in response:
            gene = None
            # save the response
            pubchem = response.text.strip()
            pubchem = pubchem.split("\n")[0]
            #print(pubchem)
            type = "small molecule"
            # url to pubchem to get synonyms, using the CID
            syn_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+pubchem+"/synonyms/json"
            syn_response = requests.get(syn_url, headers=headers)
            # if the connection works the synonyms are recorded from the json
            if syn_response.status_code == 200:
                syn_json = syn_response.json()
                synonyms = str(syn_json["InformationList"]["Information"][0]["Synonym"])
            # url to pubchem to get the SMILES for specific CID
            smiles_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+pubchem+"/property/smiles/TXT"
            smiles_response = requests.get(smiles_url,headers=headers)
            # if the connection works the SMILES for the CID is recorded
            if smiles_response.status_code == 200:
                smiles = smiles_response.text.strip()
                
            # url to pubchem to get the action with the CID
            action_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"+pubchem+"/JSON/"
            action_response = requests.get(action_url,headers=headers)
            # if the connection works the json is read
            if action_response.status_code == 200:
                action_json = action_response.json()
                action_sectionss = action_json["Record"]["Section"]
                # if Pharmacology and Biochemistry section found from TOCHeadings its content in String is read and if one of the words is found from the text then that action term is assigned to it
                for section in action_sectionss:
                    try:
                        if "Pharmacology and Biochemistry" == section.get("TOCHeading"):
                            action_section = section["Section"][0]["Information"]
                            action_text = str(action_section[0]["Value"]["StringWithMarkup"][0]["String"])
                            if "inhibit" in action_text or "antagonist" in action_text:
                                action = "inhibit"
                            elif "activate" in action_text or "agonist" in action_text:
                                action = "activate"
                            elif "modulate" in action_text:
                                action = "modulate"
                            elif "degrade" in action_text:
                                action = "degrade"
                            elif "cytotoxic" in action_text:
                                action = "cytotoxic"
                            else:
                                # if non of the words found then the action is unknown
                                action = "unknown"
                            break
                    # if Pharmacology and Biochemistry section not found the action is Not found
                    except:
                        action = "Not found"

        else:
            # if the modification is not small molecule the name will be read and then classified based on its name to knockout/knockdown and the gene will be recorded
            pubchem = None 
            synonyms = None
            smiles = None
            # spliting the common_name cause the first part of is is the gene and the second the modification type
            gene = common_name.split("_")[0].upper()
            # reading the modification and based on it recording its type and action
            if "_sg" in common_name:
                type = "CRISPR knockout"
                action = "knockout"
            elif "_sh" in common_name:
                type = "shRNA knockdown"
                action = "knockdown"
            elif "_si" in common_name:
                type = "siRNA knockdown"
                action = "knockdown"
            else:
                type = None
                action = None
        #print(action)
        
        # changing the null actions to not found
        if action is None:
            action = "not found"

        # adding all the collected info to the columns
        name_col.append(common_name)
        type_col.append(type)
        pubchem_col.append(pubchem)
        smiles_col.append(smiles)
        synonyms_col.append(synonyms)
        gene_col.append(gene)
        action_col.append(action)
    
    # adding all the columns together to create df which is returned to the user
    df = pd.DataFrame({"Name":name_col,"Type":type_col,"PubChem_CID":pubchem_col,"SMILES":smiles_col,"Synonyms":synonyms_col,"Gene":gene_col,"Action":action_col})

    return df

# getting the perturbation table with the nameToDF function and the list of common names of the perturbations
perturbation_table = nameToDF(perturbations) 
#print(perturbation_table)
#print(len(perturbation_table))


def addPertToSQL(perturbation_table):
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in perturbation_table.iterrows():
        perturbation = Perturbation(
            name = row["Name"],
            type = row["Type"],
            pubchem = row["PubChem_CID"],
            smiles = row["SMILES"],
            synonyms = row["Synonyms"],
            gene = row["Gene"],
            action = row["Action"])
        # the info will be added to the batch
        batch.append(perturbation)
        
        # when the number of inputs reaches the batch size the whole batch will be added to the sql table and the batch will be emptied after that
        if (i + 1) % batch_size == 0:
            session.add_all(batch)
            session.commit()
            batch = []

    # if there are some inputs left in the batch after the for loop those will be added to the sql table too
    if batch:
        session.add_all(batch)
        session.commit()

    # closing the connection to the db
    session.close()

# adding the Perturbation info to the db
addPertToSQL(perturbation_table)