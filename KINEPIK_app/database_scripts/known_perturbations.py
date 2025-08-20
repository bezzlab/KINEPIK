import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy import select, distinct
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Known_perturbations

# connecting to the old database and creating a session for it
engine1 = create_engine("sqlite+pysqlite:///chepro.db", echo=True)
con = engine1.connect()

Session1 = sessionmaker(bind=engine1)
session1 = Session1()

# reading the old sql table for perturbagens
df = pd.read_sql_table("PK_relationship",con)
#print(df)

def nameToUniprot(df):
    '''The function takes the df and changes the proteins to uniprot ids. All this is done by connecting to uniprot
    API endpoint. The function returns a df with the changed protein names.'''
    # empty columns lists that will be populated
    name_col = []
    perturbation_col = []
    uniprot_col = []
    source_col = []
    score_col = []
    # going throuhg the list protein by protein
    for i, row in df.iterrows():
        protein = row["kinase"]
        # creating the url to uniprot API endpoint by using the protein name
        id_url = "https://rest.uniprot.org/uniprotkb/search?query=gene:"+protein+"+AND+organism_id:9606+AND+reviewed:true&format=json&fields=accession,protein_name,organism_name"
        # to make it look like a normal internet search
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }

        # requesting the infomation from the url
        response = requests.get(id_url, headers=headers)
        #print(protein)
        try:
            # if the connection works uniprot id and species id collected from json and added to the lists
            if response.status_code == 200: 
                protein_json = response.json()
                results = protein_json["results"][0]
                uniprot = results["primaryAccession"]
                #print(sp_id, "\t", uniprot)
                name_col.append(protein)
                uniprot_col.append(uniprot)
                perturbation_col.append(row["perturbagen"])
                score_col.append(row["score"])
                if row["source"] == "kuster":
                    source_col.append(10)
                elif row["source"] == "discoverx":
                    source_col.append(11)
                elif row["source"] == "vendor":
                    source_col.append(12)
        # adding None values if no info found
        except:
            name_col.append(protein)
            uniprot_col.append(None)
            perturbation_col.append(row["perturbagen"])
            score_col.append(row["score"])
            if row["source"] == "kuster":
                source_col.append(10)
            elif row["source"] == "discoverx":
                source_col.append(11)
            elif row["source"] == "vendor":
                source_col.append(12)
    
    # adding the lists together to create the df that is returned
    df = pd.DataFrame({"Name":name_col, "Uniprot": uniprot_col, "Perturbation": perturbation_col, "Score": score_col, "Source": source_col})
    return df

# using the nameToUniprot to change the protein names to uniprot ids to match the sql table. 
known_perturbations = nameToUniprot(df)
print(known_perturbations)


# creating a new connection and session to the new database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

Session = sessionmaker(bind=engine)
session = Session()


def addKnownPertToSQL(known_perturbation_table):
    '''The function takes the df for Known_perturbations sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in known_perturbation_table.iterrows():
        known_perturbation = Known_perturbations(
            kinase = row["Uniprot"],
            perturbation = row["Perturbation"],
            source = row["Source"],
            score = row["Score"])
        # the info will be added to the batch
        batch.append(known_perturbation)

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
addKnownPertToSQL(known_perturbations)