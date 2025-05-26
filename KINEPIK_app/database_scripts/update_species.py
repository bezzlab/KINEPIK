import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy import select, distinct
import json
from sqlalchemy.orm import sessionmaker
import gget
import time
from database_structure import Gene

# creating the connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

# collecting the species and gene symbols from the gene table
species_col = session.query(Gene.species).all()
gene_list = session.query(Gene.symbol).all()
#print(len(species_col))

# creating a new list for the species column
sp_col = []
# dictionary that contains the species names and ids. Append this list when other species are added
sp_ids = {"homo_sapiens":9606}
# going through the species list from the gene table and changing the species name to the id
for species in species_col:
    sp = sp_ids[species[0]]
    sp_col.append(sp)

# genes are also collected in order to match the right species id to the gene
genes = []
#print(sp_col)
for g in gene_list:
    genes.append(g[0])

# combining these two together to create a df
df = pd.DataFrame({"Gene":genes, "Species":sp_col})
print(df)

#sp_df = pd.DataFrame({"Species":sp_col})

# updating the Gene table based on the new df
for idx, row in df.iterrows():
    gene_symbol = row["Gene"]
    species_id = row["Species"]
    gene_table = session.query(Gene).filter_by(symbol=gene_symbol).first()
    if gene_table:
        # only updating the species column
        gene_table.species = species_id
        session.commit()