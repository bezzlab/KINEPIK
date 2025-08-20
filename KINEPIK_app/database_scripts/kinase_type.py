import sqlalchemy
import pandas as pd
import requests
from sqlalchemy import create_engine
import json
import ast
from sqlalchemy.orm import sessionmaker
from database_structure import Interaction, Phosphosite, Effect, Perturbation_interaction, Subcell_location, Tissue_location, Protein

# creating the connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

protein_table = session.query(Protein).filter(Protein.kinase == 1).all()


def classifyKinaseTypes(kinase_group):
    '''Function returns the kinase type based on the given kinase group'''
    # linking the kinase groups to the amino acid residue types
    kinase_groups = {
        "TK": "Tyrosine",
        "TKL": "Tyrosine", 
        "AGC": "Ser/Thr",
        "CAMK": "Ser/Thr",
        "CMGC": "Ser/Thr",
        "CK1": "Ser/Thr",
        "STE": "Ser/Thr",
        "RGC": "Atypical",
        "Atypical": "Atypical",
        "Other": "Mixed/Uncategorized"
    }

    # handling special cases, if None (then Unknown type) or if the group is a list (then the first one chosen)
    if kinase_group is None:
        kinase_type = "Unknown"
    elif "[" in kinase_group:
        kinase_group = ast.literal_eval(kinase_group)
        kinase_group = kinase_group[0]
    
    # matching the group to the type with the dictionary
    if kinase_group is not None:
        kinase_type = kinase_groups.get(kinase_group)
    
    # returning the type
    return kinase_type

# going through all the kinases in the Protein table and getting the type and updating the db
for row in protein_table:
    kinase_type = classifyKinaseTypes(row.kinase_group)
    uniprot = row.id
    protein_sql = session.query(Protein).filter_by(id=uniprot).first()
    if protein_sql:
        protein_sql.kinase_type = kinase_type
        session.commit()
