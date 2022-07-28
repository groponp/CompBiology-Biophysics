#! __author__ = 'Ropón-Palacios G.' 
#! __date__   = 'November 19, 2021'

import MDAnalysis as mda 
import optparse 
import os 
import warnings
warnings.filterwarnings("ignore")

info= """\
+++ E M G I: EASY MAKE GROMACS INDEX +++
"""

#! CLI 
parser = optparse.OptionParser(description=info) 
parser.add_option("--coord", help="coord [PDB]", type=str)
parser.add_option("--select", help="selection based into MDAnalysis [\"segid PROA and name CA\"]", type=str, action='store') 
parser.add_option("--name_index", help="group name of index based [\"PROA_CA\"]", type=str, action='store') 
parser.add_option("--dummy_index", help="true only if run for first time [ true or t]", type=str) 



options, args = parser.parse_args() 

#! INPUTS 
pdb = options.coord
sel = options.select 
name = options.name_index 

#! MAKE INDEX gromacs
if options.dummy_index == "true" or options.dummy_index == "t": 
	os.system("echo \"q\\n\" | gmx make_ndx -f em.gro -o index_array.ndx") 


#! RUN 
u = mda.Universe(pdb) 
index = u.select_atoms(sel) 

with mda.selections.gromacs.SelectionWriter('index_array.ndx', mode='a') as ndx:
     ndx.write(index, name=name)
