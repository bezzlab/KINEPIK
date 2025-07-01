import sqlalchemy
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy import Table, Column, Integer, String
from sqlalchemy import ForeignKey
from sqlalchemy.orm import declarative_base
from sqlalchemy import Numeric

# creating a connection to the db
engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

# defining the base
Base = declarative_base()

# creating ORM classes for each sql table and specifying each column (name and type)

class Protein(Base):
    __tablename__= "Protein"
    id = Column("Uniprot_ID", String(30), primary_key=True)
    name = Column("Uniprot_name", String(30))
    kinase = Column("Kinase", Integer)
    kinase_type = Column("Kinase_type", String(30))
    kinase_group = Column("Kinase_group", String(50))
    kinase_family = Column("Kinase_family", String(50))
    kinase_subfamily = Column("Kinase_subfamily", String(50))
    species = Column("Species", Integer)
    length = Column("Length", Integer)
    sequence = Column("Sequence", String(1000))
    go = Column("GO-ID", String(1000))
    description = Column("Description", String(1000))
    gene_synonyms = Column("Gene_synonyms", String(100))
    gene = Column("Mapped_gene", String(30), ForeignKey("Gene.Symbol"))

class Interaction(Base):
    __tablename__ = "Interaction"
    id = Column(Integer, primary_key=True, autoincrement=True)
    source = Column("Source_Uniprot", String(30), ForeignKey("Protein.Uniprot_ID"))
    target = Column("Target_Uniprot", String(30), ForeignKey("Protein.Uniprot_ID")) # should this link to the protein table or the interaction or phosphosite?
    phosphosite_id = Column("Phosphosite_ID", String(30), ForeignKey("Phosphosite.Phosphosite_ID"))
    modification = Column("Modification", String(50))
    species = Column("Species", Integer)
    references = Column("Sources", String(200))

class Effect(Base):
    __tablename__ = "Effect"
    id = Column(Integer, primary_key=True, autoincrement=True)
    source = Column("Source_Uniprot", String(30), ForeignKey("Protein.Uniprot_ID"))
    target = Column("Target_Uniprot", String(30), ForeignKey("Protein.Uniprot_ID")) # should this link to the protein table or the interaction or phosphosite?
    direction = Column("Direction", Integer)
    stimulate = Column("Stimulate", Integer)
    inhibit = Column("Inhibit", Integer)
    con_direction = Column("Consensus_direction", Integer)
    con_stimulate = Column("Consensus_stimulate", Integer)
    con_inhibit = Column("Consensus_inhibit", Integer)
    species = Column("Species", Integer)
    references = Column("Sources", String(200))

class Phosphosite(Base):
    __tablename__ = "Phosphosite"
    id = Column(Integer, primary_key=True, autoincrement=True)
    phosphosite_id = Column("Phosphosite_ID", String(30))
    uniprot_id = Column("Uniprot_ID", String(30))
    residue = Column("Residue", String)
    location = Column("Location", Integer)
    confidence = Column("Confidence", String(30))
    gene = Column("Gene_name", String(30), ForeignKey("Gene.Symbol"))
    species = Column("Species", Integer)                                                            

class Gene(Base):
    __tablename__ = "Gene"
    symbol = Column("Symbol", String(30), primary_key=True)
    full_name = Column("Full_name", String(50))
    description = Column("Description", String(1000))
    uniprot = Column("Uniprot_ID", String(30))
    ensembl = Column("Ensembl_ID", String(30))
    ncbi = Column("NCBI_ID", String(30))
    species = Column("Species", Integer)

class Cell_line(Base):
    __tablename__ = "Cell_line"
    name = Column("Name", String(30), primary_key=True)
    tissue = Column("Organ", String(50))
    disease = Column("Disease", String(50))
    species = Column("Species", Integer)

class Perturbation(Base):
    __tablename__ = "Perturbation"
    name = Column("Name", String(50), primary_key=True)
    type = Column("Type", String(30))
    pubchem = Column("PubChem_CID", String(30))
    smiles = Column("SMILES", String(100))
    synonyms = Column("Synonymns", String(1000))
    gene = Column("Gene", String(30), ForeignKey("Gene.Symbol"))
    action = Column("Action", String(30))

class Perturbation_interaction(Base):
    id = Column(Integer, primary_key=True, autoincrement=True)
    __tablename__ = "Perturbation_interaction"
    perturbation = Column("Perturbation", String(30), ForeignKey("Perturbation.Name"))
    target = Column("Target_protein", String(30), ForeignKey("Protein.Uniprot_ID"))
    score = Column("Score", Numeric(20, 15))
    cell_line = Column("Cell_line", String(30), ForeignKey("Cell_line.Name"))
    references = Column("Source", String(100), ForeignKey("Study.id"))

class Subcell_location(Base):
    __tablename__ = "Subcell_location"
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column("Protein name", String(30))
    protein = Column("Uniprot", String(30), ForeignKey("Protein.Uniprot_ID"))
    localisation = Column("Localisation", String(30))
    compartment = Column("Compartment", String(30))
    confidence_l = Column("Confidence_score_L", Numeric(25, 20))
    confidence_c = Column("Confidence_score_C", Numeric(25, 20))
    cell_line = Column("Cell_line", String(30), ForeignKey("Cell_line.Name"))
    species = Column("Species", Integer)
    references = Column("Source", String(100), ForeignKey("Study.id"))

class Tissue_location(Base):
    __tablename__ = "Tissue_location"
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein = Column("Uniprot", String(30), ForeignKey("Protein.Uniprot_ID"))
    tissue = Column("Tissue", String(100))
    level = Column("Level", String(30))
    status = Column("Status", String(30))
    organ = Column("Organ", String(50))
    species = Column("Species", Integer)
    references = Column("Source", String(100))

class Experimental(Base):
    __tablename__ = "Experimental"
    id = Column(Integer, primary_key=True, autoincrement=True)
    perturbation = Column("Perturbation", String(30), ForeignKey("Perturbation.PubChem_CID"))
    phosphosite = Column("Phosphosite", String(30))
    cell_line = Column("Cell_line",String(30), ForeignKey("Cell_line.Name"))
    fc = Column("Fold_change", Numeric(30, 25))
    p_value = Column("p_value", Numeric(30, 25))
    sid = Column("SID_score", Numeric(30, 25))
    references = Column("Source", String(100), ForeignKey("Study.id"))

class Study(Base):
    __tablename__ = "Study"
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column("Name", String(50))
    date = Column("Date", String(10))
    author = Column("Author", String(100))
    desc = Column("Description", String(500))
    link = Column("Link", String(100))

# generating the metadata for the db based on the ORMs
Base.metadata.create_all(engine)
