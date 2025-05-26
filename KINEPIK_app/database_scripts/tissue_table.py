import sqlalchemy
import pandas as pd
import requests
from sqlalchemy import create_engine
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Tissue_location,Protein

# creating the connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

Session = sessionmaker(bind=engine)
session = Session()

def fromAPItoPandas(url):
    '''The function takes the url of the API endpoint and changes the format of it from json to df'''
    # to make it look like a normal internet search
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }

    # requesting the infomation from the url
    response = requests.get(url, headers=headers)

    # if the connection works
    if response.status_code == 200:
        # save the response as json and then change it to df
        x = response.json()

        df = pd.DataFrame(x)
        return df
    
    # if the connection does not work print error message with the status code
    else:
        return f"Error: {response.status_code}"



def combineTissueData(protein_list):
    '''Function takes a list of uniprot protein ids and collects tissue information from Omnipath regarding the tissue the proteins is linked to. 
    The function returns a df with the info about the tissues'''
    # empty df that will be populated later
    all_data = pd.DataFrame()
    # going through all the proteins in the list
    for protein in protein_list:
        # the protein is added to the query to access specific info. Since the protein is tuple [0] is added to access the protein id
        url = "https://omnipathdb.org/annotations?format=json&databases=HPA_tissue&proteins="+protein[0]
        # using fromAPItoPandas function with the above url to get the API respunse to df format
        one_data = fromAPItoPandas(url)
        #print(HPA_loc_data)
        # if the connection was sunccessful the given df is pivoted since in Omnipath the info is in multiple rows. Theses rows are combined to one row woth pivot_table()
        if isinstance(one_data, pd.DataFrame) and not one_data.empty:
            one_data = one_data.pivot_table(
                index=['uniprot', 'genesymbol', 'entity_type', 'source', 'record_id'],
                columns='label',
                values='value',
                aggfunc='first'
            ).reset_index()

            #print(one_data)
            #print(one_data["level"].unique())
            #print(one_data["status"].unique())

            # the new df with one protein is added to the all_data df which will eventually have all the info
            all_data = pd.concat([all_data,one_data])
            # since the two dfs can have same protein multiple times they are dropped and also the index is reset
            all_data = all_data.reset_index(drop=True).drop_duplicates()
        
        # if the connection to API was not successful warning message is printed
        else:
            print(f"Failed with protein: {protein}")
    
    # the df with all the proteins is returned
    return all_data

# collecting all the protein id from the Protein table
protein_list = session.query(Protein.id).distinct().all()
#print(protein_list)
# running the combineTissueData function with the above protein list
tissue_df = combineTissueData(protein_list)
#print("here")
#print(tissue_df["uniprot"].unique())

def speciesID (length,sp):
    '''The function takes the lenght of the wanted df as int and the species id as string. It will create a new df that has a species column that is as 
    long as the original df'''
    species_col = []
    for i in range(length):
        species_col.append(sp)
    df = pd.DataFrame({"species":species_col})
    return df

# creating the species df with the speciesID function 
sp_df = speciesID(len(tissue_df),"9606")
# collecting the wanted columns from the tissue_df and combining them with the species column
tissue_table = tissue_df[["uniprot","source","level","organ","status","tissue"]]
tissue_table = pd.concat([tissue_table,sp_df],axis=1)
#print(tissue_table)

def addTissueToSQL(tissue_table):
    '''The function takes the df for Tissue_location sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session  
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in tissue_table.iterrows():
        tissue = Tissue_location(
            protein = row["uniprot"],
            tissue = row["tissue"],
            level = row["level"],
            status = row["status"],
            organ = row["organ"],
            species = row["species"],
            references = row["source"]
        )
        # the info will be added to the batch
        batch.append(tissue)

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

# adding the Tissue_location info to the db
addTissueToSQL(tissue_table)