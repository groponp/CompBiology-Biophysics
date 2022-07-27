#! Calculate Root Mean Square Fluctuation (RMSF). 
#! _author : Ropón-Palacios G. 
#! _date   : Feb 27, 2022. 
#! _e-mail : groponp@gamil.com 

import optparse 
import subprocess 
import MDAnalysis as mda 
import numpy as np 
import pandas as pd 

disclaimer="""<MDAnalysis RMSD Tool>"""
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
def rmsf(traj, seltext1, ofile): 
  from MDAnalysis.analysis import rms
  u1 = traj 
  bb = u1.select_atoms(seltext1)
  r = rms.RMSF(bb).run() 
  resids = bb.resids 
  rmsf  = r.rmsf 
  #larray = np.array((resids,rmsf))
  df = pd.DataFrame({"resid":resids, "rmsf":rmsf})
  df.to_csv(ofile, index=False, sep="\t")
  #! rmsf to beta PDB 
  u1.add_TopologyAttr('tempfactors')    # add empty attribute for all atoms
  protein = u1.select_atoms('protein')  # select protein atoms
  for residue, r_value in zip(protein.residues, r.rmsf):
    residue.atoms.tempfactors = r_value
  
  protein.write("rmsf_to_beta.pdb")

  print ("file RMSF write to: ", ofile)     

#! Call routine
selAtoms = options.sel
ofile = options.ofile 
rmsf(universe, selAtoms, ofile)

print ("RMSF calculations finished!!!") 






 	       
