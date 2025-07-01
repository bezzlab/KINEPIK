import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
import json
from sqlalchemy.orm import sessionmaker
from database_structure import Experimental, Phosphosite, Subcell_location

# creating a connection and session to the database
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

Session = sessionmaker(bind=engine)
session = Session()


#stu_list = session.query(Subcell_location.references).all()
#print(stu_list)

subloc = pd.read_sql_table(Subcell_location.__tablename__,session.bind)
print(subloc)

# updating the Gene table based on the new df
for idx, row in subloc.iterrows():
    study = 0
    uniprot = row["Uniprot"]
    loc_table = session.query(Subcell_location).filter_by(protein=uniprot).first()
    if loc_table:
        # only updating the species column
        loc_table.references = study
        session.commit()
        