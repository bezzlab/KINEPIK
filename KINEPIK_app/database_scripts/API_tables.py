import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Interaction, Phosphosite, Effect


def fromAPItoPandas(url):
    '''The function takes the url of the API endpoint and changes the format of it from json to df'''
    # to make it look like a normal internet search
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
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
    

# the url to omnipath
asso_url = "https://omnipathdb.org/annotations?format=json&databases=SignaLink_pathway,NetPath,SIGNOR"
inter_url = "https://omnipathdb.org/interactions/?datasets=kinaseextra&organisms=9606&fields=sources&format=json"
enz_url = "https://omnipathdb.org/enz_sub?genesymbols=1&fields=sources&organisms=9606&format=json"


def speciesID (url,length):
    '''The function takes the url of the API endpoint and the length of the df. It returns the species column
    that is populated with the species id based on the url and is repeated as many times as the length of the given df'''
    # spliting the url from where the taxonomy id is found and selecting the id from the splits
    organism_part = url.split("organisms=")
    organism = organism_part[1].split("&")[0]
    # then based on all the info above the species_col is filled with the for loop
    species_col = []
    for i in range(length):
        species_col.append(organism)
    return species_col

# creating the pandas dataframes with fromAPItoPandas function
interactions = fromAPItoPandas(inter_url)
enzyme_sub = fromAPItoPandas(enz_url)

# selecting the wanted columns that are fine as they are from the API
interaction1 = enzyme_sub[["enzyme","substrate","modification"]]

def phosphositeIDs(df):
    '''The function takes the df and then retuns a list with the phosphosite ids that are based on the info found from
    the other columns in the df'''
    phosphosite_list = []
    # the for loop is run based on the length of the df
    for index in range(len(df)):
        # selecting the specific gene and residue info from the df and changing it to str so that it can be use in sql table
        gene = str(df["substrate_genesymbol"].iloc[index])
        residue_type = str(df["residue_type"].iloc[index])
        residue_loc = str(df["residue_offset"].iloc[index])
        # combining the variables together so that they are in right format for the phosphosite id
        phosphosite_id = gene + "(" + residue_type + residue_loc + ")"
        # the id is then added to the list and once the for loop is done the list is retuned
        phosphosite_list.append(phosphosite_id)
        #print(phosphosite_id)
    return phosphosite_list

#print(len(enzyme_sub))
#print(len(phosphosite_list))
#print(len(species_list))

# species id list created with the speciesID function and the phosphosite id list with phosphositeIDs function
species_inter = speciesID(enz_url, len(enzyme_sub))
phosphosites_inter = phosphositeIDs(enzyme_sub)

# these two list are combined into a df
interaction2 = pd.DataFrame({
    "Phosphosite_ID":phosphosites_inter,
    "Species":species_inter}
)

# the column names of previous df are renamed in order to recognise them better later
interaction1 = interaction1.rename(columns = {"enzyme":"Source_Uniprot","substrate":"Target_Uniprot","modification":"Modification"})

def sourcesToString(sources):
    '''The function takes the source column from a df and changes the source lists from list format in the rows to 
    strings and returns it as a df'''
    source_list = []
    # goes through the column row by row and changes the content of it into a string for sql
    for row in sources:
        source = str(row)
        source_list.append(source)

    # creating a new df from the new list and returning that
    source_df = pd.DataFrame({"Sources":source_list})
    return source_df

# changing the source lists from list format to string
inter_sources = sourcesToString(enzyme_sub["sources"])
#print(interaction1.head(10))
#print(interaction2.head(10))

# combining all the dfs together to form the interaction table with all the needed info. The index is used as guidance when combining everything
interaction_table = pd.concat([interaction1, interaction2,inter_sources], axis=1)
#print(interaction_table.columns)
#print(interaction_table.tail(10))



##### Effect table
# species id list created with the speciesID function and creating a df out of it
effect_species = speciesID(inter_url,len(interactions))
effect1 = pd.DataFrame({"Species":effect_species})

# selecting the wanted columns that are fine as they are from the API
effect2 = interactions.rename(columns = {"source":"Source_Uniprot", "target":"Target_Uniprot", "is_directed":"Direction", "is_stimulation":"Stimulate", "is_inhibition":"Inhibit", "consensus_direction":"Consensus_direction", "consensus_stimulation":"Consensus_stimulate", "consensus_inhibition":"Consensus_inhibit", "sources":"Sources"})
# changing the source lists from list format to string
effect_sources = sourcesToString(effect2["Sources"])
# dropping the old Sources column from the df to not have duplicates
effect2 = effect2.drop(columns=["Sources"])
#print(effect2.columns)
# combining the dfs together in order to have the effect table with all the needed info for the sql table
effect_table = pd.concat([effect2,effect1,effect_sources], axis=1)
#print(effect_table.columns)

##### Phosphosite table
# selecting the wanted columns that are fine as they are from the API
phosphosites = enzyme_sub[["substrate", "substrate_genesymbol", "residue_type", "residue_offset"]]

# species id list created with the speciesID function and the phosphosite id list with phosphositeIDs function
phosphosite_ids = phosphositeIDs(enzyme_sub)
species_phosphosites = speciesID(enz_url, len(enzyme_sub))

# adding the above lists together into a df
phosphosite1 = pd.DataFrame({
    "Phosphosite_ID":phosphosite_ids,
    "Species":species_phosphosites}
)

# the column names of previous df are renamed in order to recognise them better later
phosphosites2 = phosphosites.rename(columns = {"substrate":"Uniprot_ID","substrate_genesymbol":"Gene_name","residue_type":"Residue","residue_offset":"Location"})

# combining the dfs together in order to have the phosphosite table with the needed info for the sql table
phosphosite_table = pd.concat([phosphosite1,phosphosites2], axis=1) 
#print(len(phosphosite_table))
# removing all the duplicate rows from the df
phosphosite_table = phosphosite_table.drop_duplicates()
##print(phosphosite_table)
#print(len(phosphosite_table))

#print(interaction_table.columns)


def addInteractionToSQL(interaction_table):
    '''The function takes the df for the Interaction sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size 
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in interaction_table.iterrows():
        interaction = Interaction(
            source = row["Source_Uniprot"],
            target = row["Target_Uniprot"],
            phosphosite_id = row["Phosphosite_ID"],
            modification = row["Modification"],
            species = row["Species"],
            references = row["Sources"]
        )
        # the info will be added to the batch
        batch.append(interaction)
        
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

# adding the interaction info to the db
addInteractionToSQL(interaction_table)



def addEffectToSQL(effect_table):
    '''The function takes the df for Effect sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session    
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in effect_table.iterrows():
        effect = Effect(
            source = row["Source_Uniprot"],
            target = row["Target_Uniprot"],
            direction = row["Direction"],
            stimulate = row["Stimulate"],
            inhibit = row["Inhibit"],
            con_direction = row["Consensus_direction"],
            con_stimulate = row["Consensus_stimulate"],
            con_inhibit = row["Consensus_inhibit"],
            species = row["Species"],
            references = row["Sources"]
        )
        # the info will be added to the batch
        batch.append(effect)
        
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

# adding the effect info to the db
addEffectToSQL(effect_table)





def addPhosphositeToSQL(phosphosite_table):
    '''The function takes the df for Phosphosite sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session  
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in phosphosite_table.iterrows():
        phosphosite = Phosphosite(
            phosphosite_id = row["Phosphosite_ID"],
            uniprot_id = row["Uniprot_ID"],
            residue = row["Residue"],
            location = row["Location"],
            gene = row["Gene_name"],
            species = row["Species"]  
        )
        # the info will be added to the batch
        batch.append(phosphosite)
        
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

# adding the Phosphosite info to the db
addPhosphositeToSQL(phosphosite_table)