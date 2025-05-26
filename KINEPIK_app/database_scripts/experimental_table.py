import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy import select, distinct
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Base, Study, Experimental

# connecting to the database and creating a session  
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

Session = sessionmaker(bind=engine)
session = Session()

# use these if need to create the tables again. The lines drops the table and then creates them again 
# Experimental.__table__.drop(engine)
# Base.metadata.create_all(engine)

def expDF(file_name):
    '''The function takes the file name and collects all the needed info for the Experimental sql 
    table and retuns it as df'''
    # reading the csv file into a df
    data = pd.read_csv(file_name, sep="\t")
    # adding all the needed columns from the df to a new df
    exp_df1 = data[["pst","perturbagen","fc","pval_eb","sid_score"]]

    # creating list for reference and cell line columns
    ref = []
    cell_line = []
    # using for loop to make sure the length of the lists match the lenght of the df
    for i in range(len(exp_df1)):
        # use this line to get the study ids from the sql table  (note slow method if thousands of lines of data)
        #study_id = session.query(Study.id).filter_by(name=file_name).first()

        # since file names known the references are created based on them (note faster method but requires manual editing of the code)
        if file_name == "ctamdb_dpoa_HL60.tsv":
            ref.append(7)
        elif file_name == "ctamdb_dpoa_MCF7.tsv":
            ref.append(8)
        elif file_name == "ctamdb_dpoa_NTERA2.tsv":
            ref.append(9)
        #ref.append(study_id[0])
        # collecting the cell line info from the file name
        cl = file_name.split("_")
        cl = cl[2].split(".")[0]
        cell_line.append(cl)

    # combining the lists together into a df
    exp_df2 = pd.DataFrame({"Source":ref,"Cell_line":cell_line})

    # combining the dfs together and returning the df with all the info for the sql table
    experimental_table = pd.concat([exp_df1,exp_df2], axis=1)
    return experimental_table

# creating dfs for all the experimental files
exp_HL60 = expDF("ctamdb_dpoa_HL60.tsv")
exp_MCF7 = expDF("ctamdb_dpoa_MCF7.tsv")
exp_NTERA2 = expDF("ctamdb_dpoa_NTERA2.tsv")


def addExperimentalToSQL(experimental_table):
    '''The function takes the df for Experimental sql table. It creates a connection the the db and populates it 
    based on the give df'''
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in experimental_table.iterrows():
        experimental = Experimental(
            perturbation = row["perturbagen"],
            phosphosite = row["pst"],
            cell_line = row["Cell_line"],
            fc = row["fc"],
            p_value = row["pval_eb"],
            sid = row["sid_score"],
            references = row["Source"],)
        # the info will be added to the batch
        batch.append(experimental)
    
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

# adding all the Experimental info to the db
addExperimentalToSQL(exp_HL60)
addExperimentalToSQL(exp_MCF7)
addExperimentalToSQL(exp_NTERA2)