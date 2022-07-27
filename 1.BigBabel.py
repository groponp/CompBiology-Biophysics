#########################################################
# BigBabel:                                             #
# Is an object-orient script to perform several task in #
# large database of ligand from SDF files.              #
#                                                       #
# 	@autor: Ropón-Palacios G.                           # 
# 	date: 14, Jul 2022.                                 # 
#	e-mail: groponp@gmail.com                           #
#########################################################

#! Import all library 
#========================
import sys 
import glob
from openbabel import pybel 
import os 
import time 
import pandas as pd 
import numpy as np 


def search_abspath(directory): 
	'''
	Set parameters:
		- pass Current/Here directory name.
		- return absolute path to a give directory
	'''
	path = directory 
	abspath = os.path.abspath(path) 
	return abspath
	print("The Path is: {}".format(abspath))
	time.sleep(5) 

def read_and_split_sdf(path, molregex="*.sdf"):
	'''
	Read SDF files, split and convet to Smiles
	   - path : Is absolute path to directory containing files 
	   - molregex : Is an regular expresión to match with desired file. i.e "*.sdf". 
	'''
	db_list = glob.glob(path +"/"+ molregex)
	count = 0
	for database in db_list:
	    base_name1 = "tmp_" + database.split(".")[0]
	    base_name2 = base_name1.split("/")[6]
	    dirname = path+"/"+base_name2
	    if not os.path.exists(dirname):
	        os.mkdir(dirname)
	    for sdf in pybel.readfile("sdf", database): 
	        save_path = os.path.join(dirname, "L{}.smi".format(count)) 
	        sdf.write("smi", save_path, overwrite=True)
	        count += 1
	        print("Mol#{} convert to {}".format(count,save_path))
        
def make_dataframe_and_save(path, subdir_regex="*", remove_smi_db=True):
	''' Make Dataframe with smiles
		- path : absolute path to directory containing subdirs with smi files.
		- subdir_regex : Regular expresion to match with subfolder, i.e "*"
		- remove_smir_db : Defaul True to delete subdir of SMILES to free space. 
	'''
	names = []     #! Name of smiles
	smiles = []    #! Smiles character
	fromdb = []    #! Data base source 

	subdirs = glob.glob(path + "/"+ subdir_regex)
	subdirs
	for subdir in subdirs:
	    if os.path.isdir(subdir) and subdir != PATH+"/"+"ligand_pdb":
	        os.chdir(subdir)
	        smiles_list = glob.glob("*")
	        for smile in smiles_list:
	            name = smile.split(".")[0]
	            #print("Name is {}".format(name)) 
	            names.append(str(name))
	            with open(smile) as f:
	                lines = f.read()
	                smi = lines.split()
	                if type(smi) is list:
	                	smi= smi[0]
	                else:
	                	smi = smi 
	                #print("Smi string to add is {}".format(smi))
	                smiles.append(str(smi)) 
	                fromdb.append(str(subdir.split("/")[-1]))
	        os.chdir("../")
	        if remove_smi_db == True:
	        	os.system("rm -rf {}".format(subdir))  
	        else:
	        	print("Warning you're storing SMILES database.")

	#! Save BidData to smi and csv files
	#======================================
	df0= pd.DataFrame( { "Mol_Name" : names, "Smiles" : smiles, "DataBase" : fromdb})
	df1 = pd.DataFrame({ "Smiles" : smiles, "Mol_Name" : names})

	df0.to_csv('BigDatabase0.csv', index=False)
	i=500
	nrows=len(df1.index)     #! Total number of rows in dataframe. 
	for j in range(round((nrows/500)+1)):
	    df1[i-500:i].to_csv('BigDatabase1_' + str(i) +'.smi', index=False, sep="\t", header=False) 
	    i=i+500

def convert_smi_to_pdb(path, db, odir, minpH=7.4, maxpH=7.4, pka=0, maxconf=5, numproc=8):
	'''
	Convert BigDatabase to PDB
		- values taken from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7656723/ 
		- db : Name of database into SMI format. 
		- path : Absolute path where you are working.
		- odir : Output directory where will be stored ligand PDB. 
		- minpH and maxpH : Is range of pH to ionizing ligands.
		- pka : Is range of pka to an given pH.
		- maxconf : Is max number of conformomers to generate by ligand.
		- numproc : Is number of cores that you want use to make it task. 
	'''
	
	if str(os.getcwd()) == str(path):
		if not os.path.exists(odir):    
			os.mkdir(path+"/"+odir)
		if os.path.getsize(db) == 0:
			print("File has no smile data")
		else: 
	    		command = "python $GYPSUM/run_gypsum_dl.py --source  {}  --min_ph {} --max_ph {} ".format(db, minpH, maxpH) +  \
	    		"--pka_precision {} --max_variants_per_compound {} --output_folder {} --add_pdb_output ".format(pka, maxconf, odir+"/") + \
	    		"--separate_output_files --use_durrant_lab_filters --job_manager multiprocessing " + \
	    		"--num_processors {}".format(numproc)
	    		os.system(command) 
	    		os.system("rm -rf {}".format(path+"/"+odir+"/*.sdf")) 


def convet_ligand_from_pdb_to_pdbqt(directory, fpythosh, fpreplig, ligregex="*.pdb"):  
	'''
		Convert PDB to PDBQT
			- directory : Name of directory where're stored ligand convert into pdb 
			- fpythonsh : Absolute path to pythosh exec from AutodockTools.
			- fprelig : Absolute paht to prepare_ligand4.py from AutodokcTools
	'''
	workdir = os.getcwd()
	pdbs = glob.glob(workdir +"/" + directory + "/" + ligregex)

	path_pythonsh = fpythosh       #! Add absolute PATH to pythonsh from ADT
	path_prepare_lig = fpreplig        #! same to prepare_ligand4.py from ADT 

	os.mkdir(workdir+"/"+"ligand_pdbqt")
	for pdb in pdbs:
		if os.path.getsize(pdb) == 0:
			print("No data into PDB")

		else:
			name = pdb.split(".")[0]
			os.system("{exec} {script} -v -l {ligand} -o {dir}/{output}.pdbqt".format(exec=path_pythonsh,
				script=path_prepare_lig, ligand=pdb, dir="ligand_pdbqt", output=name)) 


#! Calling desired routines/functions
#==================================================
PATH = search_abspath("Nectin4") 
read_and_split_sdf(path=PATH, molregex="*.sdf")
make_dataframe_and_save(path=PATH, subdir_regex="*", remove_smi_db=True)

smis_list = glob.glob(PATH + "/*.smi")
for i in smis_list:
	convert_smi_to_pdb(path=PATH, db=i, odir="ligand_pdb", minpH=7.4, maxpH=7.4, pka=0, maxconf=5, numproc=6)
	os.remove(i)









