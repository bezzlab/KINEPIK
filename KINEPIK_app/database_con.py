from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker,scoped_session
from KINEPIK_version import Version

# connection to the version history db
version_engine = create_engine("sqlite+pysqlite:///KINEPIK_version.db", echo=True)
version_con = version_engine.connect()

version_Session = sessionmaker(bind=version_engine)
version_session = version_Session()

def getSessionForVersion(api_version:int):
    '''Return SQLAlchemy session based the requested version api/version'''
    version = version_session.query(Version).filter(Version.id==api_version).first()
    if not version:
        raise ValueError(f"No database found for version {api_version}")
    db = version.db
    # creating the connection to the database
    database_url = "sqlite:///"+ db +".db"
    engine = create_engine(database_url, echo = True)

    # scoped sessions creates more secure sessions thta prevents sharing sessiosnbetween queries and users
    session_local = scoped_session(sessionmaker(bind=engine))
    return session_local