from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker,scoped_session

# creating the connection to the database
database_url = "sqlite:///KINEPIK.db"
engine = create_engine(database_url, echo = True)

# scoped sessions creates more secure sessions thta prevents sharing sessiosnbetween queries and users
session_local = scoped_session(sessionmaker(bind=engine))