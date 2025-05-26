import sqlalchemy
import pandas as pd
import requests
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from database_structure import Base, Phosphosite

# reading the excel file to df
data = pd.read_excel("pr2c00131_si_010.xlsx", sheet_name="Table S10. Gold standard set")
# getting the list of proteins and sites from the df
protein_list = data["Protein"].tolist()
residue_list = data["Site"].tolist()
# acquiring the column that has the info about location
loc_list = data["Position in target protein"]
print(len(protein_list))
print(len(residue_list))
print(len(loc_list))

def uniprotToGene(protein_list):
    '''The function takes the list of protein ids and by using uniprot API collects gene names for the proteins. Gene names are added to a list
    and returned to the user'''
    gene_col = []
    # going through the protein list protein by protein
    for protein in protein_list:
        # url to the uniprot API endpoint with the given protein id
        id_url = "https://rest.uniprot.org/uniprotkb/"+protein+"?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym"
        # to make it look like a normal internet search
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }

        # requesting the infomation from the url
        response = requests.get(id_url, headers=headers)
        #print(protein)
        # if the connection works reading the json to get the gene name
        try:
            if response.status_code == 200: 
                protein_json = response.json()
                gene = protein_json["genes"][0]["geneName"]["value"]
                gene_col.append(gene)
        except:
            gene_col.append(None)
    
    return gene_col

# getting the gene list with uniprotToGene function
gene_list = uniprotToGene(protein_list)
#print(gene_list)


def phosphositeIDs(gene_list, residue_list, loc_list):
    '''The function takes the gene, residue and residue location lists and returns  a list of phosphosite ids that are generated based on the given lists'''
    phosphosite_list = []
    #species_list = []
    # going through the whole list
    for i in range(len(gene_list)):
        # if the gene name is not None a phosphosite id will be generated for it
        if gene_list[i] != None:
            # accessing the specific variable based on the i
            # changing all of them to str so that they can be combined together later
            gene = str(gene_list[i])
            residue_type = str(residue_list[i])
            residue_loc = str(loc_list[i])
            # combining above info together in right format for phosphosite id
            phosphosite_id = gene + "(" + residue_type + residue_loc + ")"
            #species = speciesID(enz_url)
            # adding the id to the column list
            phosphosite_list.append(phosphosite_id)
            #species_list.append(species)
            #print(phosphosite_id)
    # returning the generated list
    return phosphosite_list

# using phosphositeIDs function to get the phosphosite ids
confidence_phosphosites = phosphositeIDs(gene_list,residue_list,loc_list)


# creating a connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

Session = sessionmaker(bind=engine)
session = Session()

# going through the phophosite ids and if the site is found from the sql table its confidence is updated to be high
for site in confidence_phosphosites:
    phosphosite = session.query(Phosphosite).filter_by(phosphosite_id=site).first()
    if phosphosite:
        phosphosite.confidence = "High"
        session.commit()