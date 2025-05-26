import sqlalchemy
import pandas as pd
import requests
import os
from sqlalchemy import create_engine
from sqlalchemy import text
import json
from sqlalchemy.orm import sessionmaker
import numpy as np
from database_structure import Base, Cell_line, Perturbation_interaction

# use this part if new tables needs to be created, otherwise will append the existing tables
# engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

# Session = sessionmaker(bind=engine)
# session = Session()

# # note by deleting cell line table you will also delete the manually added cell lines
# Cell_line.__table__.drop(engine)
# Perturbation_interaction.__table__.drop(engine)
# Base.metadata.create_all(engine)

def removeNA(list,na):
    '''The function take a list and the na in string form. na in the list will be removed and list clear of them will be returned'''
    new_list = []
    # the for loop saves the values that are not na
    for x in list:
        if x != na:
            new_list.append(x)
    
    return new_list

def gctFileToPandas(file_name, sheet_name, reference):
    '''The function takes the excel file name, the sheet name and the reference id. The excel file is in gct file format and it will parsed through to format 
    df that contains all the needed info for the sql table. The function returns the cleared df'''
    # opening the excel file and the sheet that contains the needed info
    raw_data = pd.read_excel(file_name, sheet_name=sheet_name)
    #print(raw_data)
    # removing the first empty row
    raw_data = raw_data.drop(0, axis=0)
    #print(raw_data)

    # selecting the new first row as the headers and creating a new df with the new headers
    headers = raw_data.iloc[0]
    clean_data = pd.DataFrame(raw_data.values[1:], columns=headers)
    #print(clean_data.columns)

    # transposing the df to collect the info on the rows of the gct file for new headers and creating a new df out of them
    transpose_clean = clean_data.transpose()
    new_headers = transpose_clean.iloc[0]
    #print(new_headers)
    transpose_clean = pd.DataFrame(transpose_clean.values[1:],columns=new_headers)
    #print(transpose_clean)
    #print(transpose_clean.columns)

    # selecting the common names of the perturbations and the cell lines from the new df
    pert_id = transpose_clean.loc[:,"pert_iname"]
    cell_lines = transpose_clean.loc[:,"cell_id"]
    #print("pert ids")
    #print(pert_id)

    # collecting the uniprot ids from the not transposed df
    uniprot_ids = clean_data.loc[:,"pr_uniprot_id"]
    #print(uniprot_ids)

    # creating lists out of all the info collected above
    pert_id_list = pert_id.to_list()
    uniprot_ids_list = uniprot_ids.to_list()
    cell_lines_list = cell_lines.to_list()

    # removing nas from the lists since columns that have the uniprot ids, common names of the perturbations and cell lines have approx 20 nas at the beginning before the data
    uniprot_ids_list = removeNA(uniprot_ids_list,"na")
    pert_id_list = removeNA(pert_id_list,"na")
    cell_lines_list = removeNA(cell_lines_list,"na")

    # to remove the excess rows and only have the data and the ids the lenghts and shapes of the data is needed
    uni_len = len(uniprot_ids_list)
    pert_len = len(pert_id_list)
    clean_uni_len = clean_data.shape[0]
    clean_pert_len = clean_data.shape[1]

    # then the rows that will be excluded at the beginning (for rows and columns) are calculated based on the shape and length of the data
    remove_rows = clean_uni_len - uni_len
    remove_columns = clean_pert_len - pert_len

    # removing the unnecessary rows
    uniprot_ids = uniprot_ids.iloc[remove_rows:]
    cell_lines = cell_lines.iloc[remove_rows:]
    #print(cell_lines)
    #print(uniprot_ids)
    #print(remove_columns)
    #print(remove_rows)

    # removing the unnecessary data from the rows and columns from the not transposed df
    only_data = clean_data.iloc[remove_rows:,remove_columns:]

    #uniprot_df = pd.DataFrame({"Uniprot_ID":uniprot_ids_list})
    #print(uniprot_df)
    # collecting the current column names
    data_ids = only_data.columns
    done_df = only_data
    #print(data_ids)
    i = 0
    # replace the test ids with pubchem ids that where acquired earlier. Also adding the cell line if so that it is not lost 
    for name in data_ids:
        new_column = str(pert_id_list[i])+"-"+str(cell_lines_list[i])
        done_df = done_df.rename(columns = {name:new_column})
        i += 1


    # combining the pertubation data and the uniprot ids to one table
    done_df = pd.concat([uniprot_ids,done_df], axis=1)

    # making sure that all the column names are in string format
    done_df.columns = done_df.columns.astype(str)

    # changing the format of the df from wide format to long format where uniprot id, perturbation name and the related value are in the same row one at the time
    done_df = done_df.melt(id_vars=["pr_uniprot_id"])
    #print(done_df)

    # from the new long format df selecting the perturbation names and the linked cell line info and creating a list out of that column
    variable_col = done_df.loc[:,"variable"].to_list()
    #print(variable_col)

    # going through all the variable names and separating the cell line and perturbation name info to their respective columns
    id_list = []
    cell_list = []
    for name in variable_col:
        id = name.split("-")[0]
        # if the name if from the small molecule data the names are capitalised to have constant form for the names, also if any extra info (eg where the perturbation came from) added to the name (after _) will be removed
        if "chem" in file_name:
            id = id.split("_")[0]
            id = id.capitalize()
        #print(id)
        cell = name.split("-")[1]
        #print(cell)
        id_list.append(id)
        cell_list.append(cell)


    # based on the length of the other columns a reference column list will be created with the ref id that was given at the beginning
    length = len(id_list)
    ref = []
    for z in range(length):
        ref.append(reference)

    #print(length)
    #print(len(ref))
    
    # all the info is added to a new complete df and it is returned to the user
    perturbation_interaction_table = pd.DataFrame({"Perturbation":id_list, "Target_protein": done_df["pr_uniprot_id"], "Score": done_df["value"], "Cell_line":cell_list, "Source":ref})

    #perturbation_interaction_table

    #print(last_data)

    #perturbation_interaction_table = pd.concat([done_df,last_data], axis=1)

    #print(perturbation_interaction_table)
    return perturbation_interaction_table

# creating dfs for both gct files
pert_chem = gctFileToPandas("Data/Per_chem_P100.xlsx", "Chem data", "3")
#print(pert_chem)
pert_crispr = gctFileToPandas("Data/Per_CRISPR_P100.xlsx", "CRISPR data", "2")
#print(pert_crispr)

def pertCSV(dataset, sep, protein_meta, reference):
    '''The function takes the dataset name, the character that the dataset uses as separator, the metadata file name and the reference id, 
    and returns all the data as a df that is formated for the sql table'''
    # reading the files into dfs
    pert_data = pd.read_csv(dataset, sep=sep)
    pert_protein_meta = pd.read_csv(protein_meta,sep="\t")

    # the length of the df used to generate other columns
    length = len(pert_data)

    # creating a dictionary of the protein names and the uniprot ids
    name_to_id = dict(zip(pert_protein_meta["pr_name"],pert_protein_meta["verf_pr_uniprot_id"]))
    #print(id_to_name)

    # generating the reference column with the ref id
    ref = []
    for z in range(length):
        ref.append(reference)

    # getting the uniprot ids in to a column list based on the names in the file and using the previous dictionary as guide
    ids = []
    for name in pert_data["Protein Name"]:
        #print(name)
        id = name_to_id.get(name)
        ids.append(id)
        
    #print(ids)
    # creating a df out of the id list and then adding all the data together in one df that is returned to the user
    ids_df = pd.DataFrame({"Target_protein":ids})
    df = pd.DataFrame({"Perturbation":pert_data["Small Molecule Name"], "Target_protein": ids_df["Target_protein"], "Score": pert_data["% Control"], "Source":ref})
    return df

# creating dfs out of the other files
pert_data1 = pertCSV("Data/LDS-1507/Data/20342.csv",",","Data/LDS-1507/Metadata/Protein_Metadata.txt","6")
pert_data2 = pertCSV("Data/LDS-1505/Data/20340.csv",",","Data/LDS-1505/Metadata/Protein_Metadata.txt","5")
pert_data3 = pertCSV("Data/LDS-1110/Data/20124.txt","\t","Data/LDS-1110/Metadata/Protein_Metadata.txt","4")

def normaliseScore(scores, type):
    '''The function takes the scores column from a df and the information regarding the type of the scores (either logfc or %control). 
    It changes the control values to logfc so that the scores in the table would be harmonised. The function returns the scores as a new df'''
    # change the format to list from df
    score_list = scores.to_list()

    # change the scores to logfc
    transformed_score = []
    # if the scores are already logfc then they will be just added to the score list
    if type == "logfc":
        for score in score_list:
            transformed_score.append(score)
    # if the score type is control, the fc value will be calculated and then the new logfc values are added to the list
    elif type == "control":
        for score in score_list:
            fc = score/100
            # if the fc value is not 0 it is passed to log2 and then added to the list
            if fc != 0:
                logfc = np.log2(fc)
                transformed_score.append(logfc)
            # since log2 cannot handle 0s, if the fc value is 0, -10 is appointed to it because that signs that inhibition is drastic
            else:
                transformed_score.append(-10)
    # creating a df out of the list of the scores and returning it to the user
    df = pd.DataFrame({"Score":transformed_score})
    return df

# pert_data1_inhibit = normaliseScore(pert_data1["Score"],"control")

def changeScore(oldDf,type):
    '''Function uses the normaliseScore function to normalise the scores and replacing the old results with the new. 
    The function takes the df and the type of the scores and returns a new df'''
    scoreDf = normaliseScore(oldDf["Score"],type)
    newDf = oldDf.drop(columns=["Score"])
    df = pd.concat([newDf,scoreDf], axis=1)
    return df

#changing the scores of the dfs
pert_data1_done = changeScore(pert_data1,"control")
#print(pert_data1_done)
pert_data2_done = changeScore(pert_data2,"control")
#print(pert_data2_done)
pert_data3_done = changeScore(pert_data3,"control")
#print(pert_data3_done)
pert_crispr_done = changeScore(pert_crispr,"logfc")
#print(pert_crispr_done)
pert_chem_done = changeScore(pert_chem,"logfc")
#print(pert_chem_done)

def cellLineToPandas(dataset, sheet_name):
    '''The function creates a df with the cell line information based on the excel dataset and data sheet the user provided and returns the df to the user'''
    # collecting the data from the excel sheet
    raw_cell_line_data = pd.read_excel(dataset, sheet_name=sheet_name)
    # selecting the species info from the df and changing the format to list
    species_col = raw_cell_line_data["cl_organism"].tolist()
    species_list = []
    # dictionary to change the species latin name to taxonomy id. Append this manually when extending to other species
    sp_dict = {"homo sapiens":9606}
    # creating the list of the species ids
    for species in species_col:
        species = species.lower()
        species_list.append(sp_dict[species])
    # adding all the info into a new df and returning it to the user
    cell_line_table = pd.DataFrame({"Name":raw_cell_line_data["cl_center_canonical_id"],"Organ":raw_cell_line_data["cl_organ"],"Disease":raw_cell_line_data["cl_disease"],"Species":species_list})
    return cell_line_table

# collecting the celline info
pert_chem_cell = cellLineToPandas("Data/Per_chem_P100.xlsx", "Chem meta_cell_line")
pert_crispr_cell =cellLineToPandas("Data/Per_CRISPR_P100.xlsx", "CRISPR meta_cell_line")
#print(pert_chem_cell)
#print(pert_crispr_cell)

def addPertInterToSQL(pert_inter_table):
    '''The function takes the df for Perturbation_interaction sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session  
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

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

# adding the Perturbation interaction info to the db
addPertInterToSQL(pert_data1_done)
addPertInterToSQL(pert_data2_done)
addPertInterToSQL(pert_data3_done)
addPertInterToSQL(pert_crispr_done)
addPertInterToSQL(pert_chem_done)

def addCellLineToSQL(cell_line_table):
    '''The function takes the df for Cell_line sql table. It creates a connection the the db and populates it 
    based on the give df'''
    # connecting to the database and creating a session 
    engine = create_engine("sqlite+pysqlite:///KINEPIK.db", echo=True)

    Session = sessionmaker(bind=engine)
    session = Session()

    # defining the batch size
    batch_size = 500
    batch = []

    # using a for loop and iterrows to go through all the rows in the df and then add that info for the variable
    for i, row in cell_line_table.iterrows():
        cell_line = Cell_line(
            name = row["Name"],
            tissue = row["Organ"],
            disease = row["Disease"],
            species = row["Species"]
        )
        # the info will be added to the batch
        batch.append(cell_line)

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

# adding the cell line info to the db
addCellLineToSQL(pert_chem_cell)
addCellLineToSQL(pert_crispr_cell)