import sqlalchemy
import pandas as pd
import requests
from sqlalchemy import create_engine
import json
import datetime
import ast
from sqlalchemy.orm import sessionmaker
import time
from KINEPIK_app.database_scripts.database_structure import Interaction, Phosphosite, Effect, Perturbation_interaction, Subcell_location, Tissue_location, Protein, Known_perturbations, Experimental,Perturbation,Cell_line,Study,Gene
from KINEPIK_version import Version
import gget


######CONNECTIONS TO THE DATABASES
# Before running this script run database_structure.py to create an empty database file (need to edit the name of the database, aka add corresponding version title eg KINEPIK_v3.db)


# creating the connection and session to the version database
version_engine = create_engine("sqlite+pysqlite:///KINEPIK_version.db", echo=True)
version_con = version_engine.connect()

version_Session = sessionmaker(bind=version_engine)
version_session = version_Session()

old_version = version_session.query(Version).order_by(Version.id.desc()).first()
print(old_version.id)
old_db = old_version.db
new_id = old_version.id+1
new_db = "KINEPIK_v"+str(new_id)
new_api = "v1" # change this if api version changed
new_date = datetime.datetime.now().strftime("%Y-%m-%d")


# creating the connection and session to the previous database
# old_engine = create_engine("sqlite+pysqlite:///"+old_db+".db", echo=True)
old_engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)
old_con = old_engine.connect()

old_Session = sessionmaker(bind=old_engine)
old_session = old_Session()

# creating a connection to the new database
engine = create_engine("sqlite+pysqlite:///"+new_db+".db", echo=True)
con = engine.connect()

Session = sessionmaker(bind=engine)
session = Session()

# add new version to the version db
new_record = Version(id=new_id, date=new_date, api=new_api, db=new_db)
session.add(new_record)
# commit the changes to the database
session.commit()


#### ADDING TABLES THAT DO NOT USE APIS DIRECTLY TO THE NEW DB
#adding the tables that are not updated through apis directly to the new db
pert_df = pd.read_sql_table("Perturbation",old_con)
pertinter_df = pd.read_sql_table("Perturbation_interaction",old_con)
exp_df = pd.read_sql_table("Experimental",old_con)
cell_df = pd.read_sql_table("Cell_line",old_con)
study_df = pd.read_sql_table("Study",old_con)
known_df = pd.read_sql_table("Known_perturbations",old_con)
subcell_df = pd.read_sql_table("Subcell_location",old_con)

def addPertToSQL(perturbation_table):
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in perturbation_table.iterrows():
        perturbation = Perturbation(
            name = row["Name"],
            type = row["Type"],
            pubchem = row["PubChem_CID"],
            smiles = row["SMILES"],
            synonyms = row["Synonyms"],
            gene = row["Gene"],
            action = row["Action"])
        # the info will be added to the batch
        batch.append(perturbation)
        
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
# adding the Perturbation info to the new db
addPertToSQL(pert_df)

def addPertInterToSQL(pert_inter_table):
    '''The function takes the df for Perturbation_interaction sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in pert_inter_table.iterrows():
        cell_line = row.get("Cell_line",None)
        pert_interaction = Perturbation_interaction(
            perturbation = row["Perturbation"],
            target = row["Target_protein"],
            score = row["Score"],
            cell_line = cell_line,
            references = row["Source"]
        )
        # the info will be added to the batch
        batch.append(pert_interaction)
        
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
# adding the Perturbation interaction info to the new db
addPertInterToSQL(pertinter_df)

def addExperimentalToSQL(experimental_table):
    '''The function takes the df for Experimental sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in experimental_table.iterrows():
        experimental = Experimental(
            perturbation = row["Perturbation"],
            phosphosite = row["Phosphosite"],
            cell_line = row["Cell_line"],
            fc = row["Fold_change"],
            p_value = row["p_value"],
            sid = row["SID_score"],
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
# adding all the Experimental info to the new db
addExperimentalToSQL(exp_df)

def addCellLineToSQL(cell_line_table):
    '''The function takes the df for Cell_line sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in cell_line_table.iterrows():
        cell_line = Cell_line(
            name = row["Name"],
            tissue = row["Organ"],
            disease = row["Disease"],
            species = row["Species"]
        )
        session.add(cell_line)
        session.commit()
    
    # closing the connection to the db
    session.close()
# adding the cell line info to the new db
addCellLineToSQL(cell_df)

def addStudyToSQL(study_table):
    '''The function takes the df for Cell_line sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in study_table.iterrows():
        study = Study(
            id = row["id"],
            name = row["Name"],
            date = row["Date"],
            author = row["Author"],
            desc = row["Description"],
            link = row["Link"]
        )
        session.add(study)
        session.commit()
    
    # closing the connection to the db
    session.close()
# adding the study info to the new db
addStudyToSQL(study_df)

def addKnownPertToSQL(known_perturbation_table):
    '''The function takes the df for Known_perturbations sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in known_perturbation_table.iterrows():
        known_perturbation = Known_perturbations(
            kinase = row["Kinase"],
            perturbation = row["Perturbation"],
            source = row["Source"],
            score = row["Score"])
        # the info will be added to the batch
        batch.append(known_perturbation)

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
# the known perturbations are added to the new sql table
addKnownPertToSQL(known_df)

def addLocToSQL(loc_table):
    '''The function takes the df for Subcell_location sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in loc_table.iterrows():
        location = Subcell_location(
            name = row["Protein name"],
            protein = row["Uniprot"],
            localisation = row["Localisation"],
            compartment = row["Compartment"],
            confidence_l = row["Confidence_score_L"],
            confidence_c = row["Confidence_score_C"],
            cell_line = row["Cell_line"],
            species = row["Species"],
            references = row["Source"])
        # the info will be added to the batch
        batch.append(location)
        
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
# adding the Subcell_location info to the new db
addLocToSQL(subcell_df)

##### STARTING THE UPDATE WITH PROTEINS
# upload a complete list of human proteins from uniprot

def collectProteins():
    '''The function collects all the protein infos from protein table and compares it to the list of proteins from Uniprot to add all the proteins 
    to the database. And then returns them as a list to the user'''
    # selecting all the proteins from the protein table
    old_proteins = old_session.query(Protein.id).distinct().all()
    #print(old_proteins)
    
    # selecting the uniprot ids from the file that has all the human proteins and creating a list from them
    all_proteins_tsv = pd.read_csv("uniprotkb_proteome_UP000005640_AND_revi_2025_08_17.tsv",sep = "\t")
    all_proteins = all_proteins_tsv["Entry"].tolist()
    #print(all_proteins)

    old_proteins_list = []
    for old in old_proteins:
        old_proteins_list.append(old[0])
    #print(old_proteins_list)

    protein_list = []
    duplicates = []
    # go through all the proteins in the list
    for protein in all_proteins:
        # if the protein is not already in the protein_list, it will be added to it
        if protein not in old_proteins_list and protein != None:
            protein_list.append(protein)
        else:
            duplicates.append(protein)

    #print(protein_list)
    return protein_list, duplicates

# save protein list from the function to the variable
protein_list, duplicates = collectProteins()

print(protein_list)
print(len(duplicates))

def fromAPItoPandas(url):
    '''The function takes the url of the API endpoint and changes the format of it from json to df'''
    # to make it look like a normal internet search
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }

    # requesting the infomation from the url
    response = requests.get(url, headers=headers)

    # if the connection works
    if response.status_code == 200:
        # save the response as json
        x = response.json()
        # save it as a df
        df = pd.DataFrame(x)
        return df
    
    # if the connection does not work the status code returned
    else:
        return f"Error: {response.status_code}"
    
def addProteinToSQL(protein_table, db_version):
    '''The function takes the df for Protein sql table. It creates a connection the the db and populates it 
    based on the give df'''

    db_path = "sqlite+pysqlite:///"+db_version+".db"
    # connecting to the new database and creating a session  
    engine = create_engine(db_path, echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in protein_table.iterrows():
        protein = Protein(
            id = row["Uniprot_ID"],
            name = row["Uniprot_name"],
            kinase = row["Kinase"],
            kinase_group = row["Kinase_group"],
            kinase_family = row["Kinase_family"],
            kinase_subfamily = row["Kinase_subfamily"],
            species = row["Species"],
            length = row["Length"],
            sequence = row["Sequence"],
            description = row["Description"],
            gene_synonyms = row["Gene_synonyms"],
            gene = row["Mapped_gene"])
        # the info will be added to the batch
        batch.append(protein)
        
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

if len(protein_list) > 0:
    # to collect the kinase data connect to omnipath and use the fromAPItoPandas function to get the df
    kinases_url = "https://omnipathdb.org/annotations?format=json&databases=kinase.com"
    kinases_df = fromAPItoPandas(kinases_url)

    # using pivot_table to align the rows with group/family/subfamily info for specific protein and creating a new clean df
    pivoted_kinases_df = kinases_df.pivot_table(
        index=['uniprot', 'genesymbol', 'entity_type', 'source', 'record_id'],
        columns='label',
        values='value',
        aggfunc='first'
    ).reset_index()

    print(2)

    #print(pivoted_kinases_df.columns)
    # empty columns for the protein df (for kinase infos)
    kinase_column = []
    family = []
    group = []
    subfamily = []
    # using for loop to go through all the proteins in proteins list
    for protein in protein_list:
        # selecting the row from the pivoted_kinases_df based on the protein uniprot id
        row = pivoted_kinases_df[pivoted_kinases_df['uniprot'] == protein]
        #print(row)
        # if the length of the row is 0, pivoted_kinases_df does not have that protein and thus it is not kinase
        if len(row) == 0:
            # 0 for this column since the column uses boolean to indicate if the protein is kinase or not
            kinase_column.append(0)
            family.append(None)
            group.append(None)
            subfamily.append(None)
        # if the lenght is longer than 0 then the protein is kinase and the values are assigned to each columns
        else:
            #row = pivoted_kinases_df[pivoted_kinases_df['uniprot'] == protein]
            #print(row)
            # 1 for this column since the column uses boolean to indicate if the protein is kinase or not
            kinase_column.append(1)
            # some proteins have multiple families/groups/subfamilies due to contradictory sources and thus this info is presented as lists
            if len(row) > 1:
                family.append(str(row["family"].dropna().unique().tolist()))
                group.append(str(row["group"].dropna().unique().tolist()))
                subfamily.append(str(row["subfamily"].dropna().unique().tolist()))
            # if just one row of info then those values are used
            else:
                #print(row)
                #print(row.at[0,"family"])
                family.append(row["family"].values[0])
                group.append(row["group"].values[0])
                subfamily.append(row["subfamily"].values[0])

    # creating a df with the generated columns
    protein_df1 = pd.DataFrame({
        "Uniprot_ID":protein_list,
        "Kinase":kinase_column,
        "Kinase_group":group,
        "Kinase_family":family,
        "Kinase_subfamily":subfamily
    })

    #print("df1 len:"+ str(len(protein_df1)))
        
    #test = pivoted_kinases_df[pivoted_kinases_df['uniprot'] == "P31751"]
    #print(len(kinase_column))
    #print(len(family))
    #print(len(group))
    #print(len(subfamily))

    print(3)

    # collecting other info about the proteins
    name_col = []
    species_col = []
    length_col = []
    sequence_col = []
    description_col = []
    synonyms_col = []
    gene_col = []

    # again going through the list of proteins one by one
    for protein in protein_list:
        print(protein)
        # to make it look like a normal internet search
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
        # connecting to uniprot to collect the info about the protein, using the protein id in the url to get specific info
        uniprot_url = "https://rest.uniprot.org/uniprotkb/"+ protein +"?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym%2Corganism_name%2Clength%2Csequence"

        # requesting the infomation from the url
        response = requests.get(uniprot_url, headers=headers)

        # if the connection to GWAS works
        if response.status_code == 200:
            # save the response as json
            protein_data = response.json()
        
        # collecting all the info from the json and saving them to different variables
            try:
                uni_name = protein_data["uniProtkbId"]
                #print(uni_name)
                species = int(protein_data["organism"]["taxonId"])
                #print(species)
                length = int(protein_data["sequence"]["length"])
                #print(length)
                seq = protein_data["sequence"]["value"]
                #print(seq)
                desc = protein_data["proteinDescription"]["recommendedName"]["fullName"]["value"]
                #print(desc)
                gene = protein_data["genes"][0]["geneName"]["value"]
                #print(gene)
                genes = protein_data["genes"]

                # since there are multiple synonyms this collection is done for loop to collect all the synonymns and not just the first one
                syn_list = []
                for gen in genes:
                    for name in gen.get("synonyms",[]):
                        syn = name.get("value")
                        syn_list.append(syn)
                str_synonyms = str(syn_list)
            #print(str_synonyms)

            # if there are errors None is assigned to the protein
            except (KeyError, IndexError):
                uni_name =None
                species = None
                length = None
                seq = None
                desc = None
                gene = None
                str_synonyms = None

            # adding all the collected info to the columns
            name_col.append(uni_name)
            species_col.append(species)
            length_col.append(length)
            sequence_col.append(seq)
            description_col.append(desc)
            synonyms_col.append(str_synonyms)
            gene_col.append(gene)
            
        # if the connection does not work the status code returned
        else:
            print(f"Error: {response.status_code} with protein {protein}")

    # print(len(name_col))
    # print(len(species_col))
    # print(len(length_col))
    # print(len(sequence_col))
    # print(len(description_col))
    # print(len(synonyms_col))
    # print(len(gene_col))

    # combining other info to a new df
    protein_df2 = pd.DataFrame({
        "Uniprot_name":name_col,
        "Species":species_col,
        "Length":length_col,
        "Sequence":sequence_col,
        "Description":description_col,
        "Gene_synonyms":synonyms_col,
        "Mapped_gene":gene_col
    })

    # combining kinase info and other info into one protein table
    protein_table = pd.concat([protein_df1,protein_df2], axis=1)
    print(protein_table)

    # adding the Protein info to the db
    addProteinToSQL(protein_table,new_db)

elif len(protein_list) == 0:
    df = pd.read_sql_table("Protein",old_con)
    addProteinToSQL(df, new_db)


# adding the kinase type in the db
protein_table_kinases = session.query(Protein).filter(Protein.kinase == 1).all()


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
for row in protein_table_kinases:
    kinase_type = classifyKinaseTypes(row.kinase_group)
    uniprot = row.id
    protein_sql = session.query(Protein).filter_by(id=uniprot).first()
    if protein_sql:
        protein_sql.kinase_type = kinase_type
        session.commit()

# select all the proteins from the database
sql_proteins = session.query(Protein.id).distinct().all()


# create empty lists for the columns that will be populated
go_col = []
name_col = []

# the for loop will go through all the proteins that we collected from the database
for protein in sql_proteins:
    # to make it look like a normal internet search
    headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, kuten Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
    
    # the url for the api endpoint that will give the needed data. Since the protein is in tuple [0] is added
    url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=UniProtKB:"+protein[0]+"&aspect=biological_process"
        
    # requesting the infomation from the url
    response = requests.get(url, headers=headers)

    # if the connection works
    if response.status_code == 200:
        
        # saving the response in json format
        protein_data = response.json()
        
        # collecting the list of dictionaries that is found under "results"-key
        go_data = protein_data["results"]

        # empty list to be populated with the go terms. This will also clear the list from the previous protein
        go_id_list = []
        # the for loop goes through all the go terms in the api endpoint json
        for result in go_data:
            # only collecting the go ids that are manually curated and not IEA(inferred electronic annotation) which are created automaticcally by computational methods
            if result["goEvidence"] != "IEA":
                go_id = result["goId"]
                go_id_list.append(go_id)
        
        # adding both the go id and the uniprot id to the lists
        name_col.append(protein[0])
        go_col.append(go_id_list)
        
            
    # if the connection does not work the status code returned
    else:
        print(f"Error: {response.status_code} with protein {protein}")
    
    # short breaks to make sure that uniprot api does not close the connection due to hundreds of queries
    time.sleep(0.2) 

#print(len(name_col))
#print(len(go_col))

# combining the lists into a df
go_df = pd.DataFrame({
    "Uniprot":name_col,
    "GO-ids":go_col
})

#print(go_df)

# going through the new dataframe row by row and collecting the results and updating the go id column based on the uniprot id
for idx, row in go_df.iterrows():
    uniprot = row["Uniprot"]
    go_term = str(row["GO-ids"])
    protein_sql = session.query(Protein).filter_by(id=uniprot).first()
    if protein_sql:
        protein_sql.go = go_term
        session.commit()

##### RETREIVING DATA FROM OMNIPATH FOR INTERACTION, EFFECT AND PHOSPHOSITE TABLES
# the url to omnipath
asso_url = "https://omnipathdb.org/annotations?format=json&databases=SignaLink_pathway,NetPath,SIGNOR"
inter_url = "https://omnipathdb.org/interactions/?datasets=kinaseextra&organisms=9606&fields=sources&format=json"
enz_url = "https://omnipathdb.org/enz_sub?genesymbols=1&fields=sources&organisms=9606&format=json"


def speciesID (url,length):
    '''The function takes the url of the API endpoint and the length of the df. It returns the species column
    that is populated with the species id based on the url and is repeated as many times as the length of the given df'''
    # spliting the url from where the taxonomy id is found and selecting the id from the splits
    organism_part = url.split("organisms=")
    organism = organism_part[1].split("&")[0]
    # then based on all the info above the species_col is filled with the for loop
    species_col = []
    for i in range(length):
        species_col.append(organism)
    return species_col

# creating the pandas dataframes with fromAPItoPandas function
interactions = fromAPItoPandas(inter_url)
enzyme_sub = fromAPItoPandas(enz_url)

# selecting the wanted columns that are fine as they are from the API
interaction1 = enzyme_sub[["enzyme","substrate","modification"]]

def phosphositeIDs(df):
    '''The function takes the df and then retuns a list with the phosphosite ids that are based on the info found from
    the other columns in the df'''
    phosphosite_list = []
    # the for loop is run based on the length of the df
    for index in range(len(df)):
        # selecting the specific gene and residue info from the df and changing it to str so that it can be use in sql table
        gene = str(df["substrate_genesymbol"].iloc[index])
        residue_type = str(df["residue_type"].iloc[index])
        residue_loc = str(df["residue_offset"].iloc[index])
        # combining the variables together so that they are in right format for the phosphosite id
        phosphosite_id = gene + "(" + residue_type + residue_loc + ")"
        # the id is then added to the list and once the for loop is done the list is retuned
        phosphosite_list.append(phosphosite_id)
        #print(phosphosite_id)
    return phosphosite_list

#print(len(enzyme_sub))
#print(len(phosphosite_list))
#print(len(species_list))

# species id list created with the speciesID function and the phosphosite id list with phosphositeIDs function
species_inter = speciesID(enz_url, len(enzyme_sub))
phosphosites_inter = phosphositeIDs(enzyme_sub)

# these two list are combined into a df
interaction2 = pd.DataFrame({
    "Phosphosite_ID":phosphosites_inter,
    "Species":species_inter}
)

# the column names of previous df are renamed in order to recognise them better later
interaction1 = interaction1.rename(columns = {"enzyme":"Source_Uniprot","substrate":"Target_Uniprot","modification":"Modification"})

def sourcesToString(sources):
    '''The function takes the source column from a df and changes the source lists from list format in the rows to 
    strings and returns it as a df'''
    source_list = []
    # goes through the column row by row and changes the content of it into a string for sql
    for row in sources:
        source = str(row)
        source_list.append(source)

    # creating a new df from the new list and returning that
    source_df = pd.DataFrame({"Sources":source_list})
    return source_df

# changing the source lists from list format to string
inter_sources = sourcesToString(enzyme_sub["sources"])
#print(interaction1.head(10))
#print(interaction2.head(10))

# combining all the dfs together to form the interaction table with all the needed info. The index is used as guidance when combining everything
interaction_table = pd.concat([interaction1, interaction2,inter_sources], axis=1)
#print(interaction_table.columns)
#print(interaction_table.tail(10))

interaction_table = interaction_table.loc[interaction_table['Modification'].isin(["phosphorylation", "dephosphorylation"])]

##### Effect table
# species id list created with the speciesID function and creating a df out of it
effect_species = speciesID(inter_url,len(interactions))
effect1 = pd.DataFrame({"Species":effect_species})

# selecting the wanted columns that are fine as they are from the API
effect2 = interactions.rename(columns = {"source":"Source_Uniprot", "target":"Target_Uniprot", "is_directed":"Direction", "is_stimulation":"Stimulate", "is_inhibition":"Inhibit", "consensus_direction":"Consensus_direction", "consensus_stimulation":"Consensus_stimulate", "consensus_inhibition":"Consensus_inhibit", "sources":"Sources"})
# changing the source lists from list format to string
effect_sources = sourcesToString(effect2["Sources"])
# dropping the old Sources column from the df to not have duplicates
effect2 = effect2.drop(columns=["Sources"])
#print(effect2.columns)
# combining the dfs together in order to have the effect table with all the needed info for the sql table
effect_table = pd.concat([effect2,effect1,effect_sources], axis=1)
#print(effect_table.columns)

##### Phosphosite table
# selecting the wanted columns that are fine as they are from the API
phosphosites = enzyme_sub[["substrate", "substrate_genesymbol", "residue_type", "residue_offset", "modification"]]
print("enzyme_sub:", len(enzyme_sub))
print("phosphosites before filter:", len(phosphosites))
print("phosphosites after filter:", len(phosphosites.loc[phosphosites['modification'].isin(['phosphorylation','dephosphorylation'])]))

phosphosites = phosphosites.loc[phosphosites['modification'].isin(["phosphorylation", "dephosphorylation"])]

# species id list created with the speciesID function and the phosphosite id list with phosphositeIDs function
phosphosite_ids = phosphositeIDs(phosphosites)
species_phosphosites = speciesID(enz_url, len(phosphosites))

# adding the above lists together into a df
phosphosite1 = pd.DataFrame({
    "Phosphosite_ID":phosphosite_ids,
    "Species":species_phosphosites}
)

# the column names of previous df are renamed in order to recognise them better later
phosphosites2 = phosphosites.rename(columns = {"substrate":"Uniprot_ID","substrate_genesymbol":"Gene_name","residue_type":"Residue","residue_offset":"Location"})

# combining the dfs together in order to have the phosphosite table with the needed info for the sql table
phosphosite_table1 = pd.concat([phosphosite1,phosphosites2], axis=1) 
#print(len(phosphosite_table))
# removing all the duplicate rows from the df
phosphosite_table1 = phosphosite_table1.drop_duplicates(subset=["Phosphosite_ID"])

print(len(phosphosite_table1))
print(phosphosite_table1['Phosphosite_ID'].nunique())
phosphosite_table1 = phosphosite_table1.drop(columns=['modification'])

### ADDING THE EXPERIMENTAL PHOSPHOSITES TO THE PHOSPHOSITE DF
def phosphositeDF(phos_list):
    '''The function takes the list of phosphosites and then retuns a df for the Phosphosite table'''
    phosphosite_col = []
    uniprot_col = []
    residue_col = []
    loc_col = []
    gene_col = []
    sp_col = []
    problem =[]
    # the for loop is run based on the length of the df
    for pho in phos_list:
        # selecting the specific gene and residue info from the df and changing it to str so that it can be use in sql table
        gene = str(pho[0].split("(")[0])
        residue_type = pho[0].split("(")[1]
        residue_type = str(residue_type[0])
        sep = "(" + residue_type
        residue_loc = pho[0].split(sep)[1]
        residue_loc = str(residue_loc.replace(")",""))
        # combining the variables together so that they are in right format for the phosphosite id
        phosphosite_id = pho[0]
        # the id is then added to the list and once the for loop is done the list is retuned
        try:
            protein_row = session.query(Protein).filter_by(gene=gene).first()
            uniprot = protein_row.id
            sp = protein_row.species
        except:
            # default if no genes found from the db
            found = False
            # list all the synonyms that are not None
            syn_list = session.query(Protein).filter(Protein.gene_synonyms != None).all()
            # go through the list row by row
            for row in syn_list:
                # try except keeps the loop going even if results not found
                try:
                    # changing the synonym list from string to an actual list
                    synonyms = ast.literal_eval(row.gene_synonyms)
                    # if the gene that was not found earlier is found from the synonym list the values are collect from the row
                    if gene in synonyms:
                        uniprot = row.id
                        sp = row.species
                        gene = row.gene
                        found = True
                except (ValueError, SyntaxError):
                    continue
            # if the gene name is not found from the synonyms then it is recorded as none and also listed as a problem gene for further info
            if not found:
                uniprot = None
                sp = 9606
                problem.append(gene)
        print(phosphosite_id)

        # adding all the variables to the column lists
        phosphosite_col.append(phosphosite_id)
        uniprot_col.append(uniprot)
        residue_col.append(residue_type)
        loc_col.append(residue_loc)
        gene_col.append(gene)
        sp_col.append(sp)

        # combining the column lists togther into a df and returning that to the user
    df = pd.DataFrame({"Phosphosite_ID":phosphosite_col, "Uniprot_ID":uniprot_col, "Residue":residue_col, "Location":loc_col, "Gene_name":gene_col, "Species":sp_col})
    return df, problem
phos_list = session.query(Experimental.phosphosite).distinct().all()

# to collect only the phosphosites from the Experimental table that are not found from Phosphosite table
#old_phos = session.query(Phosphosite.phosphosite_id).distinct().all()
old_phos = list(phosphosite_table1["Phosphosite_ID"])
new_phos = []
for phos in phos_list:
    if phos not in old_phos and phos not in new_phos:
        new_phos.append(phos)

phosphosite_table2, problem = phosphositeDF(new_phos)
print(phosphosite_table2)
print(set(problem))

phosphosite_table = pd.concat([phosphosite_table1,phosphosite_table2])

def addInteractionToSQL(interaction_table):
    '''The function takes the df for the Interaction sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size 
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in interaction_table.iterrows():
        interaction = Interaction(
            source = row["Source_Uniprot"],
            target = row["Target_Uniprot"],
            phosphosite_id = row["Phosphosite_ID"],
            modification = row["Modification"],
            species = row["Species"],
            references = row["Sources"]
        )
        # the info will be added to the batch
        batch.append(interaction)
        
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

# adding the interaction info to the db
addInteractionToSQL(interaction_table)


def addEffectToSQL(effect_table):
    '''The function takes the df for Effect sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in effect_table.iterrows():
        effect = Effect(
            source = row["Source_Uniprot"],
            target = row["Target_Uniprot"],
            direction = row["Direction"],
            stimulate = row["Stimulate"],
            inhibit = row["Inhibit"],
            con_direction = row["Consensus_direction"],
            con_stimulate = row["Consensus_stimulate"],
            con_inhibit = row["Consensus_inhibit"],
            species = row["Species"],
            references = row["Sources"]
        )
        # the info will be added to the batch
        batch.append(effect)
        
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

# adding the effect info to the db
addEffectToSQL(effect_table)


def addPhosphositeToSQL(phosphosite_table):
    '''The function takes the df for Phosphosite sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in phosphosite_table.iterrows():
        phosphosite = Phosphosite(
            phosphosite_id = row["Phosphosite_ID"],
            uniprot_id = row["Uniprot_ID"],
            residue = row["Residue"],
            location = row["Location"],
            gene = row["Gene_name"],
            species = row["Species"]  
        )
        # the info will be added to the batch
        batch.append(phosphosite)
        
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

# adding the Phosphosite info to the db
addPhosphositeToSQL(phosphosite_table)

##### adding the phosphosite confidence to the table
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


def phosphositeIDsconf(gene_list, residue_list, loc_list):
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
confidence_phosphosites = phosphositeIDsconf(gene_list,residue_list,loc_list)

# going through the phophosite ids and if the site is found from the sql table its confidence is updated to be high
for site in confidence_phosphosites:
    phosphosite = session.query(Phosphosite).filter_by(phosphosite_id=site).first()
    if phosphosite:
        phosphosite.confidence = "High"
        session.commit()


#### ADDING TISSUE DATA TO THE DB FROM OMNIPATH
def combineTissueData(protein_list):
    '''Function takes a list of uniprot protein ids and collects tissue information from Omnipath regarding the tissue the proteins is linked to. 
    The function returns a df with the info about the tissues'''
    # empty df that will be populated later
    all_data = pd.DataFrame()
    # going through all the proteins in the list
    for protein in protein_list:
        # the protein is added to the query to access specific info. Since the protein is tuple [0] is added to access the protein id
        url = "https://omnipathdb.org/annotations?format=json&databases=HPA_tissue&proteins="+protein[0]
        print(protein)
        # using fromAPItoPandas function with the above url to get the API respunse to df format
        one_data = fromAPItoPandas(url)
        #print(HPA_loc_data)
        # if the connection was sunccessful the given df is pivoted since in Omnipath the info is in multiple rows. Theses rows are combined to one row woth pivot_table()
        if isinstance(one_data, pd.DataFrame) and not one_data.empty:
            one_data = one_data.pivot_table(
                index=['uniprot', 'genesymbol', 'entity_type', 'source', 'record_id'],
                columns='label',
                values='value',
                aggfunc='first'
            ).reset_index()

            #print(one_data)
            #print(one_data["level"].unique())
            #print(one_data["status"].unique())

            # the new df with one protein is added to the all_data df which will eventually have all the info
            all_data = pd.concat([all_data,one_data])
            # since the two dfs can have same protein multiple times they are dropped and also the index is reset
            all_data = all_data.reset_index(drop=True).drop_duplicates()
        
        # if the connection to API was not successful warning message is printed
        else:
            print(f"Failed with protein: {protein}")
    
    # the df with all the proteins is returned
    return all_data

# collecting all the protein id from the Protein table
proteins_list = session.query(Protein.id).distinct().all()
#print(protein_list)
# running the combineTissueData function with the above protein list
tissue_df = combineTissueData(proteins_list)
#print("here")
#print(tissue_df["uniprot"].unique())

def speciesIDtissue (length,sp):
    '''The function takes the lenght of the wanted df as int and the species id as string. It will create a new df that has a species column that is as 
    long as the original df'''
    species_col = []
    for i in range(length):
        species_col.append(sp)
    df = pd.DataFrame({"species":species_col})
    return df

# creating the species df with the speciesID function 
sp_df = speciesIDtissue(len(tissue_df),"9606")
# collecting the wanted columns from the tissue_df and combining them with the species column
tissue_table = tissue_df[["uniprot","source","level","organ","status","tissue"]]
tissue_table = pd.concat([tissue_table,sp_df],axis=1)
#print(tissue_table)

def addTissueToSQL(tissue_table):
    '''The function takes the df for Tissue_location sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in tissue_table.iterrows():
        tissue = Tissue_location(
            protein = row["uniprot"],
            tissue = row["tissue"],
            level = row["level"],
            status = row["status"],
            organ = row["organ"],
            species = row["species"],
            references = row["source"]
        )
        # the info will be added to the batch
        batch.append(tissue)

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

# adding the Tissue_location info to the db
addTissueToSQL(tissue_table)


#### ADDING THE GENES

def collectGenes():
    '''The function collects all the gene infos from the new db and compares it to the list of genes from the old db.'''
    # selecting all the proteins from the protein table
    new_genes_prot = session.query(Protein.gene).distinct().all()
    new_genes_phos = session.query(Phosphosite.gene).distinct().all()
    new_genes_pert = session.query(Perturbation.gene).distinct().all()
    new_all_genes = new_genes_phos + new_genes_prot + new_genes_pert

    new_no_dup_genes = list(set(new_all_genes))

    old_genes = old_session.query(Gene.symbol).distinct().all()

    old_genes_list = []
    for old in old_genes:
        old_genes_list.append(old[0])
    #print(old_proteins_list)

    new_gene_list = []
    for new in new_no_dup_genes:
        new_gene_list.append(new[0])

    gene_list = []
    duplicates = []
    # go through all the proteins in the list
    for gene in new_gene_list:
        # if the protein is not already in the protein_list, it will be added to it
        if gene not in old_genes_list and gene != None:
            gene_list.append(gene)
        else:
            duplicates.append(gene)

    #print(protein_list)
    return gene_list, duplicates

# save protein list from the function to the variable
gene_list, duplicates = collectGenes()

species_ids = {9606:"homo_sapiens"} # create more efficient way of reconding all the species used
#print(species_ids.get)


# next the ensembl accession codes are collected from ensembl based on the symbol of the gene
# going through the gene list gene by gene
ensembl_list = []
for gene in gene_list:
    # since the list is full of tuples gene ids are recorded with [0]
    gene_id = gene[0]
    #print(type(gene_id))
    # and species with [1]
    #species_id = gene[1]
    # species name is required for the url so the species id is changed to name based on the above species dictionary
    species = species_ids.get(9606)
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

print("through")
