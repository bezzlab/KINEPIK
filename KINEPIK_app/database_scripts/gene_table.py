import sqlalchemy
import pandas as pd
import requests
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import gget
import time
from database_structure import Phosphosite, Protein, Gene

# creating a connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

# collecting all the genes from Protein and Phosphosite tables and combining them together
genes_prot = session.query(Protein.gene,Protein.species).distinct().all()
genes_phos = session.query(Phosphosite.gene,Phosphosite.species).distinct().all()
all_genes = genes_phos + genes_prot
#print(all_genes)
# removing the duplicates from the list
no_dup_genes = list(set(all_genes))
#print(len(all_genes))
#print(len(no_dup_genes))
# species id dictionary
species_ids = {9606:"homo_sapiens"} # create more efficient way of reconding all the species used
#print(species_ids.get)


# next the ensembl accession codes are collected from ensembl based on the symbol of the gene
# going through the gene list gene by gene
ensembl_list = []
for gene in no_dup_genes:
    # since the list is full of tuples gene ids are recorded with [0]
    gene_id = gene[0]
    #print(type(gene_id))
    # and species with [1]
    species_id = gene[1]
    # species name is required for the url so the species id is changed to name based on the above species dictionary
    species = species_ids.get(species_id)
    #print(type(species))

    # if the gene and species are not None then connection to Ensembl will be created
    if gene_id != None and species != None:
        # the url to Ensembl API endpoint is created by using the species name and gene id
        gene_url = "https://rest.ensembl.org/xrefs/symbol/"+species+"/"+gene_id+"?content-type=application/json"
        # to make it look like a normal internet search
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }

        # requesting the infomation from the url
        response = requests.get(gene_url, headers=headers)

        # if the connection works
        if response.status_code == 200:
            # save the response as json
            x = response.json()

            # going through the json until "ENSG" found as part of id result
            for i in range(len(x)):
                if "ENSG" in x[i]["id"]:
                    ensembl_id = x[i]["id"]
                    ensembl_list.append(ensembl_id)
                    # since some of the jsons can have multiple ids only the first one is recorded and then the for loop is broken
                    break

            
    
#print(ensembl_list)
#print(len(all_genes))
#print(len(ensembl_list))


# creating the empty columns for gene table
name_col = []
species_col = []
ensemlb_col = []
uniprot_col = []
description_col = []
full_col = []
ncbi_col = []

# handling the data in batches
batch_size = 50

# going through the list of ensembl id
for x in range(0, len(ensembl_list), batch_size):
    # defining which batch is used
    batch = ensembl_list[x:x + batch_size]
    try:
        # getting df with gene info for the genes in the batch from gget.info()
        df = gget.info(batch)

        # going through the df and recording the info from it to column lists
        for j in range(len(df)):
            # if records not found the try-except will handle it
            try:
                # parse the needed information from the gget df
                #one_id = gget.info([j[0]])
                ensembl = str(df.iloc[j,0])
                symbol = str(df.iloc[j,6])
                ncbi_id = str(df.iloc[j,3])
                uniprot_id = str(df.iloc[j,1])
                uniprot_desc = str(df.iloc[j,12])
                full_name = str(df.iloc[j,11]) 
                species_e = str(df.iloc[j,4]) 

                name_col.append(symbol)
                species_col.append(species_e)
                ensemlb_col.append(ensembl)
                uniprot_col.append(uniprot_id)
                description_col.append(uniprot_desc)
                full_col.append(full_name)
                ncbi_col.append(ncbi_id)
                    
            except Exception as inner_e:
                print(f"Skipped: {inner_e}")
                continue

    except Exception as e:
        print(f"{i}:{e} batch failed")
        continue
    
    #print(x)
    # time.sleep makes sure that the gget wont get banned from the APIs due to too many requests
    time.sleep(1)

# changing the species names in species column to species ids
sp_col = []
sp_ids = {"homo_sapiens":9606} # add other species to this dictionary when adding other species, also check what format gget uses for the species name
for species in species_col:
    sp = sp_ids[species]
    sp_col.append(sp)

# combing all the columns to a df
gene_table_with_dup = pd.DataFrame({"Symbol":name_col,"Full_name":full_col,"Description":description_col,"Uniprot_ID":uniprot_col,"Ensembl_ID":ensemlb_col,"NCBI_ID":ncbi_col,"Species":sp_col})

#print(gene_table_with_dup.head(10))

# if there are any duplicate genes in the df (based on symbol) then those are dropped to avoid error with sql
gene_table = gene_table_with_dup.drop_duplicates(subset=["Symbol"])


def addGeneToSQL(gene_table):
    '''The function takes the df for Gene sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session  
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in gene_table.iterrows():
        gene = Gene(
            symbol = row["Symbol"],
            full_name = row["Full_name"],
            description = row["Description"],
            uniprot = row["Uniprot_ID"],
            ensembl = row["Ensembl_ID"],
            ncbi = row["NCBI_ID"],
            species = row["Species"])
        # the info will be added to the batch
        batch.append(gene)
        
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

# adding the Gene info to the db
addGeneToSQL(gene_table)
