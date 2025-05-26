import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy import select, distinct
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Interaction, Phosphosite, Effect, Perturbation, Perturbation_interaction, Cell_line, Experimental, Protein, Gene

# connecting to the old database and creating a session for it
engine1 = create_engine("sqlite+pysqlite:///Data/chepro.db", echo=True)
con = engine1.connect()

Session1 = sessionmaker(bind=engine1)
session1 = Session1()

# reading the old sql table for perturbagens
df = pd.read_sql_table("Perturbagen",con)
#print(df)

#name_list = df["name"].tolist()

# save perturbation name list from the df to the variable
perturbations = df["name"].tolist()
#perturbations = list(set(perturbations))
#print(perturbations)
#print(len(perturbations))

# change the list to pubchem cids
def nameToDF(pert_list):
    '''The function takes the list of perturbation common names and collects all the needed info for 
    the sql table through various APIs. The function returns all the info in df'''
    # creating empty lists for each column that will be populated later
    name_col = []
    type_col = []
    pubchem_col = []
    smiles_col = []
    synonyms_col = []
    action_col = []

    # going through the list name by name and collecting the info for each column
    for name in pert_list:
        #print(name[0])
        # capitalising the names so that they all match
        common_name = name
        #print(common_name)
        # connecting to pubchem to get the pubchem cid with the common name
        id_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+common_name+"/cids/TXT"
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
        #id_file = "id_data.txt"

        # requesting the infomation from the url
        response = requests.get(id_url, headers=headers)

        # if the connection works
        if response.status_code == 200 and "Status: 404" not in response.text:
            # the perturbation is not a gene and thus None assigned to that variable
            gene = None
            # save the response
            pubchem = response.text.strip()
            #print(pubchem)
            # if the connection works the perturbation's type is small molecule
            type = "small molecule"
            # collecting the synonymns through pubchem
            syn_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+pubchem+"/synonyms/json"
            syn_response = requests.get(syn_url, headers=headers)
            # if the connection works
            if syn_response.status_code == 200:
                # reading the json output from API
                syn_json = syn_response.json()
                # accessing the synonym list in the json and saving them to a variable
                synonyms = str(syn_json["InformationList"]["Information"][0]["Synonym"])
            # collecting the smiles through pubchem API
            smiles_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+pubchem+"/property/smiles/TXT"
            smiles_response = requests.get(smiles_url,headers=headers)
            # if the connection works the smiles is read from the API output which was in text format
            if smiles_response.status_code == 200:
                smiles = smiles_response.text.strip()

        # since it was known that all the perturbations in discoverx were inhibitors the value was assigned to action automatically
        action = "inhibit"

        # adding all the variables to the columns
        name_col.append(common_name)
        type_col.append(type)
        pubchem_col.append(pubchem)
        smiles_col.append(smiles)
        synonyms_col.append(synonyms)
        action_col.append(action)

    # once all info collected, all of the lists are combined into a df and it is returned to the user
    df = pd.DataFrame({"Name":name_col,"Type":type_col,"PubChem_CID":pubchem_col,"SMILES":smiles_col,"Synonyms":synonyms_col,"Action":action_col})

    return df

# collecting the info with the nameToDF function
perturbation_table = nameToDF(perturbations) 
#print(perturbation_table)
#print(len(perturbation_table))

# creating a new connection and session to the new database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

Session = sessionmaker(bind=engine)
session = Session()

# collecting all the already existing perturbation names from the perturbation table
old_perts_sql = session.query(Perturbation.name).distinct().all()
#print(old_perts_sql)
# making a list out of the old perturbation names
old_perts = []
for pert in old_perts_sql:
    #print(pert)
    old_perts.append(pert[0])
#print(old_perts)

# empty list that will be populated with the perturbations that where found from both tables
matched = []

# going through the list of old perturbation names name by name
for name in old_perts:
    # if the name is found from both new and old lists the data in the table will be updated
    if name in perturbations:
        #print(name)
        # collecting the info from Perturbation table and if it exists then the row is updated
        perturbation = session.query(Perturbation).filter_by(name=name).first()
        if perturbation:
            # action from the new df is selected based on the name and then added to the sql table
            #synonym = perturbation_table.loc[perturbation_table["Name"]==name, "Synonyms"].iloc[0]
            action = perturbation_table.loc[perturbation_table["Name"]==name, "Action"].iloc[0]
            #print(synonym)
            perturbation.action = action
            session.commit()
            # adding the perturbations that match to matched list so that they can be removed from the df
            matched.append(perturbation)
            

# after the names have been updated the name with its row will be removed from the df
for match in matched:
    perturbation_table = perturbation_table[perturbation_table["Name"] != match.name]
#print(perturbation_table)

def addPertToSQL(perturbation_table):
    '''The function takes the df for Perturbation sql table. It creates a connection the the db and populates it 
    based on the give df'''
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

# the new perturbations are added to the sql table
addPertToSQL(perturbation_table)