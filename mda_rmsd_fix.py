#! Calculate Root Mean Square Deviation (RMSD). 
#! _author : Ropón-Palacios G. 
#! _date   : Feb 27, 2022. 
#! _e-mail : groponp@gamil.com 

import optparse 
import subprocess 
import MDAnalysis as mda 
import numpy as np 
import pandas as pd 
import warnings
warnings.filterwarnings("ignore")
from MDAnalysis.lib.log import ProgressBar

disclaimer="""<MDAnalysis RMSD Tool>"""
parser = optparse.OptionParser(description=disclaimer) 
#! INPUTS options
parser.add_option("--coord", help="coord [PDB/GRO, PSF, PARM7]", type=str)
parser.add_option("--traj", help="traj [XTC, DCD, NETCDF]", type=str)
parser.add_option("--sel", help="syntaxis-based in MDAnalysis [\"protein and name CA\"]", type=str,action='store')
#! OUTPUT options 
parser.add_option("--ofile", help="type output name [rmsd_prot_nameCA.dat]", type=str)

options, args = parser.parse_args() 

#! Load trajectory 
universe = mda.Universe(options.coord, options.traj)
ref = universe 

def rmsd(traj1, ref, seltext1, ofile): 
    from MDAnalysis.analysis import rms  
    data = rms.RMSD(traj1.select_atoms(seltext1),      #; universe to align
             ref.select_atoms(seltext1),               #; reference universe or atomgroup
             #select=seltext1,   #; group to superimpose and calculate RMSD
             ref_frame=0)      #; frame index of the reference
    data.run() 
    
    time = data.rmsd[:,1] / 1000 ## Columna 1 is time in ps , verificar esto al correr
    rmsd = data.rmsd[:,2]      ## Colimna 2 is data to rmsd 
    #larray = np.array((time,rmsd))
    
    df = pd.DataFrame({"time":time, "rmsd":rmsd})
    df.to_csv(ofile, index=False, sep="\t")
    print ("file RMSD write to: ", ofile) 

#! Call routine
selAtoms = options.sel
ofile = options.ofile 
rmsd(universe, ref, selAtoms, ofile)

print ("RMSD calculations finished!!!") 
 	       
