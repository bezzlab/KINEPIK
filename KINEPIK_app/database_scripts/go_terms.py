import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Protein
import time

# create connection and the session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

# select all the proteins from the database
sql_proteins = session.query(Protein.id).distinct().all()


# create empty lists for the columns that will be populated
go_col = []
name_col = []

# the for loop will go through all the proteins that we collected from the database
for protein in sql_proteins:
    # to make it look like a normal internet search
    headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
    
    # the url for the api endpoint that will give the needed data. Since the protein is in tuple [0] is added
    url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=UniProtKB:"+protein[0]+"&aspect=biological_process"
        
    # requesting the infomation from the url
    response = requests.get(url, headers=headers)

    # if the connection works
    if response.status_code == 200:
        
        # saving the response in json format
        protein_data = response.json()
        
        # collecting the list of dictionaries that is found under "results"-key
        go_data = protein_data["results"]

        # empty list to be populated with the go terms. This will also clear the list from the previous protein
        go_id_list = []
        # the for loop goes through all the go terms in the api endpoint json
        for result in go_data:
            # only collecting the go ids that are manually curated and not IEA(inferred electronic annotation) which are created automaticcally by computational methods
            if result["goEvidence"] != "IEA":
                go_id = result["goId"]
                go_id_list.append(go_id)
        
        # adding both the go id and the uniprot id to the lists
        name_col.append(protein[0])
        go_col.append(go_id_list)
        
            
    # if the connection does not work the status code returned
    else:
        print(f"Error: {response.status_code} with protein {protein}")
    
    # short breaks to make sure that uniprot api does not close the connection due to hundreds of queries
    time.sleep(0.2) 

#print(len(name_col))
#print(len(go_col))

# combining the lists into a df
go_df = pd.DataFrame({
    "Uniprot":name_col,
    "GO-ids":go_col
})

#print(go_df)

# going through the new dataframe row by row and collecting the results and updating the go id column based on the uniprot id
for idx, row in go_df.iterrows():
    uniprot = row["Uniprot"]
    go_term = str(row["GO-ids"])
    protein_sql = session.query(Protein).filter_by(id=uniprot).first()
    if protein_sql:
        protein_sql.go = go_term
        session.commit()