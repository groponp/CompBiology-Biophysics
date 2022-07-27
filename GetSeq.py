#! Get Sequence from PDB file and from NCBI after perform an BLASTp against Any phyllum.
#! using as query T. spiralis 5HT protein.
#! @autor: Ropón-Palacios G.
#! date: 19, Jul, 2022.
#! e-mail: groponp@gmail.com 

from Bio import SeqIO
import os
import warnings
warnings.filterwarnings("ignore")


command = """
grep -v "POPC" system.pdb | grep -v "TIP3" | grep -v "CLA" | grep -v "SOD" | grep -v "TER" | grep -v "HEADER" > system_only_protein.pdb

"""
os.system(command) 
target = "system_only_protein.pdb"


def pdb_to_seq(target, seqfile, header):
    for record in SeqIO.parse(target, "pdb-atom"):
        with open(seqfile, "w") as f:
            f.write(">"+header+"\n")
            print(record.seq)
            f.write(str(record.seq)) 

def get_pdb(pdbID, format):
    import Bio
    from Bio.PDB import PDBList 
    pdbl = PDBList()
    list = [pdbID] 
    for pdb in list: 
        pdbl.retrieve_pdb_file(pdb, pdir=".", file_format=format, overwrite=True, obsolete=False) 

    os.rename("pdb"+pdbID+".ent", "pdb"+pdbID+".pdb")


def make_remote_blastp(query, ofile, top, entrez=None, db="nr", e_val=10E-24):
   print("Performing Remote BLASTp ...") 

   if entrez != None:
       command = """
       blastp -out remote_blastp_{ofile}.tab -outfmt 7 \
       -query {query} -db {db} -evalue {e_val} -max_target_seqs {top} \
       -remote -entrez_query {entrez} -qcov_hsp_perc 90
       """.format(ofile=ofile, query=query, db=db, e_val=e_val, entrez=entrez, top=top)
       os.system(command) 
   else: 
        command = """
        blastp -out remote_blastp_{ofile}.tab -outfmt 7 \
        -query {query} -db {db} -evalue {e_val} -max_target_seqs {top} \
        -remote -qcov_hsp_perc 90
        """.format(ofile=ofile, query=query, db=db, e_val=e_val, top=top) 
        os.system(command) 


def filter_blastp_data(file):
    import pandas as pd 
    import numpy as np 

    grep = """
    grep -v '#' {} > {}
    """.format(file, file.split(".")[0]+"_fix.tab")  
    os.system(grep)  
    df = pd.read_table(file.split(".")[0]+"_fix.tab", usecols=[0,1,2,10])
    df.columns = ["QueryNameOgrn","SeqID","%Identity","E-value"]
    df1 = df[df["%Identity"] >= 80] 
    df1.to_csv(file.split(".")[0]+"_fix.csv", sep=",", index=False)

    print("Retrive #{} SeqIDs".format(len(df1.loc[:,"SeqID"].values))) 
    return(np.array(df1.loc[:,"SeqID"].values))


def get_seqs(email, SeqIDList, odir):
    from Bio import Entrez
    from Bio import SeqIO 
    Entrez.email = email

    count = 1 
    incr  = 50
    for SeqID in SeqIDList:
        print("Downloading SeqID #{} with ID {}".format(count, SeqID)) 
        if not os.path.exists(odir):  
            os.mkdir(odir) 

        if count <= incr:
            handle = Entrez.efetch(db="protein", id=SeqID, rettype="fasta")
            record = SeqIO.read(handle, "fasta")
            SeqIO.write(record, odir+"/"+SeqID+".fasta","fasta")
            count +=1

        else:
            import time
            time.sleep(60)
            incr += 50  
      
    #! Concat all files. 
    os.system("cp 5HT_Hsapiens.fasta {}".format(odir+"/"))
    os.system("cat {}*.fasta > 5HT_MSA.fasta".format(odir+"/"))
 

#! Call routines/functions
#==============================
#pdb_to_seq(target=target, seqfile="5HT_Tspiralis.fasta", header="Trichinella spiralis |5HT protein")
#get_pdb("5V54", format="pdb")  

make_remote_blastp(query="5HT_Tspiralis.fasta", ofile="nematodes", top=20, entrez="nematoda[orgn]", db="nr", e_val=10E-25)
#make_remote_blastp(query="5HT_Tspiralis.fasta", ofile="platyhelminthes", entrez="patyhelminthes[orgn]",db="nr", e_val=10E-25)
#make_remote_blastp(query="5HT_Tspiralis.fasta", ofile="all.tab", entrez=None, db="nr", e_val=10E-25)

SeqID_array1=filter_blastp_data("remote_blastp_nematodes.tab")
#SeqID_array2=filter_blastp_data("remote_blastp_platyhelminthes.tab")
#SeqID_array3=filter_blastp_data("remote_blastp_all.tab")

get_seqs(email="georcki.ropon@unmsm.edu.pe", SeqIDList=SeqID_array1, odir="FastasFiles_Nematodes") 
#get_seqs(email="georcki.ropon@unmsm.edu.pe", SeqIDList=SeqID_array2, odir="FastasFiles_Platyhelminthes") 
#get_seqs(email="georcki.ropon@unmsm.edu.pe", SeqIDList=SeqID_array3, odir="FastasFiles_all")  


#! Remove duplicate into MSA files. 
os.system("seqkit rmdup --by-seq < 5HT_MSA.fasta > 5HT_MSA_fix.fasta")  




