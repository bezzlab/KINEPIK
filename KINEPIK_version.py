import sqlalchemy
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import Table, Column, Integer, String
from sqlalchemy import ForeignKey
from sqlalchemy.orm import declarative_base
from sqlalchemy import Numeric

# creating a connection to the db
engine = create_engine("sqlite+pysqlite:///KINEPIK_version.db", echo=True)

# defining the base
Base = declarative_base()

# creating ORM classes for each sql table and specifying each column (name and type)

class Version(Base):
    __tablename__= "Version_history"
    id = Column("version_number", Integer, primary_key=True)
    date = Column("date", String(30))
    api = Column("api_version", Integer)
    db = Column("db_name", String(30))


# generating the metadata for the db based on the ORMs
Base.metadata.create_all(engine)
