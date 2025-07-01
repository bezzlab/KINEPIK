import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy import select, distinct
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Base, Subcell_location

# # use this part if new tables needs to be created, otherwise will append the existing tables
# engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

# Session = sessionmaker(bind=engine)
# session = Session()

# # note by deleting cell line table you will also delete the manually added cell lines
# Subcell_location.__table__.drop(engine)
# Base.metadata.create_all(engine)

def locCSVtoDF(file_name, sheet_name, reference):
    '''The function takes the file name, sheet name and the refernece id. It reads the file and based on the info in it, it will create a df that has most of
    the needed info for the Subcell_location sql table. The function returns a df'''
    # reading the file into a df
    data = pd.read_excel(file_name, sheet_name=sheet_name)

    # recording the length of the df
    length = len(data)

    # by using the length the reference and cell line columns are created
    ref = []
    cell_line = []
    for z in range(length):
        ref.append(reference)
        # the sheet name is used as cell line name since the dataset uses cell lines as the sheet names
        cell_line.append(sheet_name)
        
    #print(data)

    # creating empty column lists that will be populated
    name_col = []
    neighborhood_col = []
    compartment_col = []
    n_ci_col = []
    c_ci_col = []

    # going through the whole df based on the length
    for i in range(length):
        # collecting the info about the protein, neighborhood and compartment from the df
        protein = data.iloc[i].loc["Protein"]
        neighborhood = data.loc[data["Protein"]==protein, "Neighborhood Class"].iloc[0]
        compartment = data.loc[data["Protein"]==protein, "Compartment Class"].iloc[0]
        # if the neighborhood is secretory it will collect the CI value from the Secrtetory column
        if neighborhood == "Secretory":
            n_ci = data.loc[data["Protein"]==protein, "Secretory"].iloc[0]
            # same logic will be applied to the compartments where based on the name in compartment the CI value will be collected
            if compartment == "S1":
                c_ci = data.loc[data["Protein"]==protein, "S1"].iloc[0]
            elif compartment == "S2":
                c_ci = data.loc[data["Protein"]==protein, "S2"].iloc[0]
            elif compartment == "S3":
                c_ci = data.loc[data["Protein"]==protein, "S3"].iloc[0]
            elif compartment == "S4":
                c_ci = data.loc[data["Protein"]==protein, "S4"].iloc[0]
            else:
                # None will be given as value if there are no strong signals from aby of the compartments
                c_ci = None

        # same logic will be applied to other types of neighborhoods and compartments as above
        elif neighborhood == "Nuclear":
            n_ci = data.loc[data["Protein"]==protein, "Nuclear"].iloc[0]
            if compartment == "N1":
                c_ci = data.loc[data["Protein"]==protein, "N1"].iloc[0]
            elif compartment == "N2":
                c_ci = data.loc[data["Protein"]==protein, "N2"].iloc[0]
            elif compartment == "N3":
                c_ci = data.loc[data["Protein"]==protein, "N3"].iloc[0]
            elif compartment == "N4":
                c_ci = data.loc[data["Protein"]==protein, "N4"].iloc[0]
            else:
                c_ci = None

        elif neighborhood == "Cytosol":
            n_ci = data.loc[data["Protein"]==protein, "Cytosol"].iloc[0]
            if compartment == "C1":
                c_ci = data.loc[data["Protein"]==protein, "C1"].iloc[0]
            elif compartment == "C2":
                c_ci = data.loc[data["Protein"]==protein, "C2"].iloc[0]
            elif compartment == "C3":
                c_ci = data.loc[data["Protein"]==protein, "C3"].iloc[0]
            elif compartment == "C4":
                c_ci = data.loc[data["Protein"]==protein, "C4"].iloc[0]
            else:
                c_ci = None

        elif neighborhood == "Mitochondria":
            n_ci = data.loc[data["Protein"]==protein, "Mitochondria"].iloc[0]
            if compartment == "M1":
                c_ci = data.loc[data["Protein"]==protein, "M1"].iloc[0]
            elif compartment == "M2":
                c_ci = data.loc[data["Protein"]==protein, "M2"].iloc[0]
            else:
                c_ci = None

        # if there are no strong signals for the neighborhood then None is also assign to CI values for neighborhood and compartment
        else:
            n_ci = None # check that n unclass is c unclass always
            c_ci = None

        # adding all the info to the lists
        name_col.append(protein)
        neighborhood_col.append(neighborhood)
        compartment_col.append(compartment)
        n_ci_col.append(n_ci)
        c_ci_col.append(c_ci)
        #print(protein, "\t", neighborhood, "\t", compartment, "\t", n_ci, "\t", c_ci)

    # adding the lists together to create a df and and returning the df to the user
    df = pd.DataFrame({"Protein":name_col, "Neighborhood Class": neighborhood_col, "Compartment Class": compartment_col, "Neighborhood CI": n_ci_col, "Compartment CI": c_ci_col, "Cell line":cell_line, "Source":ref})
    return df

# getting subcell location info with locCSVtoDF function. Only one sheet is used in this case but others could be added too
location_data = locCSVtoDF("mmc4.xlsx", "MCF7", 0)
#print(location_data)

# creating a list of the proteins in the df
protein_list = location_data["Protein"].tolist()
#print(protein_list)
def nameToUniprotAndSp(protein_list):
    '''The function takes the protein list and changes the proteins to uniprot ids. It also collects the species ids. All this is done by connecting to uniprot
    API endpoint. The function returns a df'''
    # empty columns lists that will be populated
    name_col = []
    sp_col = []
    uniprot_col = []
    # going throuhg the list protein by protein
    for protein in protein_list:
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
                sp_id = int(results["organism"]["taxonId"])
                uniprot = results["primaryAccession"]
                #print(sp_id, "\t", uniprot)
                name_col.append(protein)
                sp_col.append(sp_id)
                uniprot_col.append(uniprot)
        # adding None values if no info found
        except:
            name_col.append(protein)
            sp_col.append(None)
            uniprot_col.append(None)
    
    # adding the lists together to create the df that is returned
    df = pd.DataFrame({"Name":name_col, "Uniprot": uniprot_col, "Species": sp_col})
    return df

# using the nameToUniprotAndSp to change the protein names to uniprot ids to match the sql table. Also acquiring the species ids
ids_df = nameToUniprotAndSp(protein_list)
#print(ids_df)

# merging the location and ids df together based on the protein names
loc_table = pd.merge(location_data, ids_df, left_on="Protein", right_on="Name", how="outer",)
# dropping the protein name column since the uniprot ids are used as identifiers
loc_table = loc_table.drop(columns="Name")
#print(loc_table.head(10))

def addLocToSQL(loc_table):
    '''The function takes the df for Subcell_location sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session  
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in loc_table.iterrows():
        location = Subcell_location(
            name = row["Protein"],
            protein = row["Uniprot"],
            localisation = row["Neighborhood Class"],
            compartment = row["Compartment Class"],
            confidence_l = row["Neighborhood CI"],
            confidence_c = row["Compartment CI"],
            cell_line = row["Cell line"],
            species = row["Species"],
            references = row["Source"])
        # the info will be added to the batch
        batch.append(location)
        
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

# adding the Subcell_location info to the db
addLocToSQL(loc_table)