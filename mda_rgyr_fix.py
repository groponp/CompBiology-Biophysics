#! Calculate Radius of Gyration (Rgyr). 
#! _author : Ropón-Palacios G. 
#! _date   : Feb 27, 2022. 
#! _e-mail : groponp@gamil.com 

import optparse 
import subprocess 
import MDAnalysis as mda 
import numpy as np 
import pandas as pd 

 
disclaimer="""<MDAnalysis Rgyr Tool>"""
parser = optparse.OptionParser(description=disclaimer) 
#! INPUTS options
parser.add_option("--coord", help="coord [GRO, PSF, PARM7]", type=str)
parser.add_option("--traj", help="traj [XTC, DCD, NETCDF or DRT]", type=str)
parser.add_option("--sel", help="syntaxis-based in MDAnalysis [\"protein and name CA\"]", type=str,action='store')
#! OUTPUT options 
parser.add_option("--ofile", help="type output name [rgyr_prot_nameCA.dat]", type=str)

options, args = parser.parse_args() 

#! Load trajectory 
universe = mda.Universe(options.coord, options.traj)

#! Define routine
def rgyr(traj, seltext1, ofile):
    rgyr = [] 
    time = [] 
    u1 = traj 
    sel = u1.select_atoms(seltext1)  
    for ts in u1.trajectory:
        time.appen(u1.trajectory.time/1000)
        rgyr.append(sel.radius_of_gyration())
    
    #larray = np.array((time,rgyr))
    df = pd.DataFrame({"time":time, "rgyr":rgyr})
    df.to_csv(ofile, index=False, sep="\t")
    print("file Rgyr write to : ", ofile)

#! Call routine
selAtoms = options.sel
ofile = options.ofile 
rgyr(universe, selAtoms, ofile)

print ("Rgyr calculations finished!!!") 






 	       
