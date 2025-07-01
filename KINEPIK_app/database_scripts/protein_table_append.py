import sqlalchemy
import pandas as pd
import requests
from sqlalchemy import create_engine
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Interaction, Phosphosite, Effect, Perturbation_interaction, Subcell_location, Tissue_location, Protein

# creating the connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

def collectProteins():
    '''The function collects all the protein infos from protein table and compares it to the list of proteins from Uniprot to add all the proteins 
    to the database. And then returns them as a list to the user'''
    # selecting all the proteins from the protein table
    old_proteins = session.query(Protein.id).distinct().all()
    #print(old_proteins)
    
    # selecting the uniprot ids from the file that has all the human proteins and creating a list from them
    all_proteins_tsv = pd.read_csv("uniprotkb_proteome_UP000005640_AND_revi_2025_06_12.tsv",sep = "\t")
    all_proteins = all_proteins_tsv["Entry"].tolist()
    #print(all_proteins)

    old_proteins_list = []
    for old in old_proteins:
        old_proteins_list.append(old[0])
    #print(old_proteins_list)

    protein_list = []
    duplicates = []
    # go through all the proteins in the list
    for protein in all_proteins:
        # if the protein is not already in the protein_list, it will be added to it
        if protein not in old_proteins_list and protein != None:
            protein_list.append(protein)
        else:
            duplicates.append(protein)

    #print(protein_list)
    return protein_list, duplicates

# save protein list from the function to the variable
proteins, duplicates = collectProteins()

#print(len(proteins))
#print(duplicates)
#print(len(proteins))
#print(len(duplicates))
print(1)

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
        # save the response as json
        x = response.json()
        # save it as a df
        df = pd.DataFrame(x)
        return df
    
    # if the connection does not work the status code returned
    else:
        return f"Error: {response.status_code}"

# to collect the kinase data connect to omnipath and use the fromAPItoPandas function to get the df
kinases_url = "https://omnipathdb.org/annotations?format=json&databases=kinase.com"
kinases_df = fromAPItoPandas(kinases_url)

# using pivot_table to align the rows with group/family/subfamily info for specific protein and creating a new clean df
pivoted_kinases_df = kinases_df.pivot_table(
    index=['uniprot', 'genesymbol', 'entity_type', 'source', 'record_id'],
    columns='label',
    values='value',
    aggfunc='first'
).reset_index()

print(2)

#print(pivoted_kinases_df.columns)
# empty columns for the protein df (for kinase infos)
kinase_column = []
family = []
group = []
subfamily = []
# using for loop to go through all the proteins in proteins list
for protein in proteins:
    # selecting the row from the pivoted_kinases_df based on the protein uniprot id
    row = pivoted_kinases_df[pivoted_kinases_df['uniprot'] == protein]
    #print(row)
    # if the length of the row is 0, pivoted_kinases_df does not have that protein and thus it is not kinase
    if len(row) == 0:
        # 0 for this column since the column uses boolean to indicate if the protein is kinase or not
        kinase_column.append(0)
        family.append(None)
        group.append(None)
        subfamily.append(None)
    # if the lenght is longer than 0 then the protein is kinase and the values are assigned to each columns
    else:
        #row = pivoted_kinases_df[pivoted_kinases_df['uniprot'] == protein]
        #print(row)
        # 1 for this column since the column uses boolean to indicate if the protein is kinase or not
        kinase_column.append(1)
        # some proteins have multiple families/groups/subfamilies due to contradictory sources and thus this info is presented as lists
        if len(row) > 1:
            family.append(str(row["family"].dropna().unique().tolist()))
            group.append(str(row["group"].dropna().unique().tolist()))
            subfamily.append(str(row["subfamily"].dropna().unique().tolist()))
        # if just one row of info then those values are used
        else:
            #print(row)
            #print(row.at[0,"family"])
            family.append(row["family"].values[0])
            group.append(row["group"].values[0])
            subfamily.append(row["subfamily"].values[0])

# creating a df with the generated columns
protein_df1 = pd.DataFrame({
    "Uniprot_ID":proteins,
    "Kinase":kinase_column,
    "Kinase_group":group,
    "Kinase_family":family,
    "Kinase_subfamily":subfamily
})

#print("df1 len:"+ str(len(protein_df1)))
    
#test = pivoted_kinases_df[pivoted_kinases_df['uniprot'] == "P31751"]
#print(len(kinase_column))
#print(len(family))
#print(len(group))
#print(len(subfamily))

print(3)

# collecting other info about the proteins
name_col = []
species_col = []
length_col = []
sequence_col = []
description_col = []
synonyms_col = []
gene_col = []

# again going through the list of proteins one by one
for protein in proteins:
    print(protein)
    # to make it look like a normal internet search
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    # connecting to uniprot to collect the info about the protein, using the protein id in the url to get specific info
    uniprot_url = "https://rest.uniprot.org/uniprotkb/"+ protein +"?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym%2Corganism_name%2Clength%2Csequence"

    # requesting the infomation from the url
    response = requests.get(uniprot_url, headers=headers)

    # if the connection to GWAS works
    if response.status_code == 200:
        # save the response as json
        protein_data = response.json()
       
       # collecting all the info from the json and saving them to different variables
        try:
            uni_name = protein_data["uniProtkbId"]
            #print(uni_name)
            species = int(protein_data["organism"]["taxonId"])
            #print(species)
            length = int(protein_data["sequence"]["length"])
            #print(length)
            seq = protein_data["sequence"]["value"]
            #print(seq)
            desc = protein_data["proteinDescription"]["recommendedName"]["fullName"]["value"]
            #print(desc)
            gene = protein_data["genes"][0]["geneName"]["value"]
            #print(gene)
            genes = protein_data["genes"]

            # since there are multiple synonyms this collection is done for loop to collect all the synonymns and not just the first one
            syn_list = []
            for gen in genes:
                for name in gen.get("synonyms",[]):
                    syn = name.get("value")
                    syn_list.append(syn)
            str_synonyms = str(syn_list)
        #print(str_synonyms)

        # if there are errors None is assigned to the protein
        except (KeyError, IndexError):
            uni_name =None
            species = None
            length = None
            seq = None
            desc = None
            gene = None
            str_synonyms = None

        # adding all the collected info to the columns
        name_col.append(uni_name)
        species_col.append(species)
        length_col.append(length)
        sequence_col.append(seq)
        description_col.append(desc)
        synonyms_col.append(str_synonyms)
        gene_col.append(gene)
        
    # if the connection does not work the status code returned
    else:
        print(f"Error: {response.status_code} with protein {protein}")

# print(len(name_col))
# print(len(species_col))
# print(len(length_col))
# print(len(sequence_col))
# print(len(description_col))
# print(len(synonyms_col))
# print(len(gene_col))

# combining other info to a new df
protein_df2 = pd.DataFrame({
    "Uniprot_name":name_col,
    "Species":species_col,
    "Length":length_col,
    "Sequence":sequence_col,
    "Description":description_col,
    "Gene_synonyms":synonyms_col,
    "Mapped_gene":gene_col
})

# combining kinase info and other info into one protein table
protein_table = pd.concat([protein_df1,protein_df2], axis=1)
print(protein_table)

print(4)

def addProteinToSQL(protein_table):
    '''The function takes the df for Protein sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session  
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in protein_table.iterrows():
        protein = Protein(
            id = row["Uniprot_ID"],
            name = row["Uniprot_name"],
            kinase = row["Kinase"],
            kinase_group = row["Kinase_group"],
            kinase_family = row["Kinase_family"],
            kinase_subfamily = row["Kinase_subfamily"],
            species = row["Species"],
            length = row["Length"],
            sequence = row["Sequence"],
            description = row["Description"],
            gene_synonyms = row["Gene_synonyms"],
            gene = row["Mapped_gene"])
        # the info will be added to the batch
        batch.append(protein)
        
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

# adding the Protein info to the db
addProteinToSQL(protein_table)