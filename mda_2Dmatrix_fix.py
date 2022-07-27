#! Calculate Root Mean Square Deviation (RMSD). 
#! _author : Ropón-Palacios G. 
#! _date   : March 9, 2022. 
#! _e-mail : groponp@gamil.com 


import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms
import numpy as np
import matplotlib.pyplot as plt
import optparse 
import subprocess   

import warnings
warnings.filterwarnings("ignore")

disclaimer="""<MDAnalysis 2D RMSD Tool>"""
parser = optparse.OptionParser(description=disclaimer) 
#! INPUTS options
parser.add_option("--coord", help="coord [PDB/GRO, PSF, PARM7]", type=str)
parser.add_option("--traj", help="traj [XTC, DCD, NETCDF]", type=str)
parser.add_option("--sel", help="syntaxis-based in MDAnalysis [\"resname LIG and not name H\"]", type=str,action='store')
#! OUTPUT options 
parser.add_option("--ofile", help="type output name [2Dmatrix_lig.svg]", type=str)


options, args = parser.parse_args() 


gro1 = options.coord
xtc1 = options.traj 
u1 = mda.Universe(gro1, xtc1) 

#! 2D matrix to itself traj 
aligner = align.AlignTraj(u1, u1, select=options.sel,
				in_memory=True).run() 
matrix = diffusionmap.DistanceMatrix(u1, select=options.sel).run() 

#! Plot.
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15

plt.imshow(matrix.dist_matrix, cmap="hsv") 
plt.xlabel("Molecular conformation [frame]")
plt.ylabel("Molecular conformation [frame]") 
plt.colorbar(label=r'RMSD [$\AA$]')
plt.tight_layout()
plt.savefig(options.ofile, dpi=300, format="svg", bbox_inches='tight') 

print("2Dmatrix itsefl done and file is save to", options.ofile) 

#! 2D matrix to two different traj 
#traj1 = mda.Universe(gro1, xtc1, in_memory=True)
#traj2 = mda.Universe(gro2, xtc2, in_memory=True)

#prmsd = np.zeros((len(traj1.trajectory),
#                 len(traj2.trajectory)))

#for i, frame in enumerate(traj1.trajectory[0:999]):
#        r = rms.RMSD(traj1, traj2, select='name CA',
#                        ref_frame=i).run()
#        prmsd[i] = r.rmsd[:, -1]
#        print(i,"\t",prmsd[i])

#plt.matshow(prmsd, cmap='hsv', interpolation='nearest')
#plt.xlabel('traj1 # frame')
#plt.ylabel('traj2 # frame')
#plt.colorbar(label=r'RMSD ($\AA$)', shrink=0.75)
#plt.rcParams["figure.figsize"] = (15,5)
#plt.rcParams.update({'font.size': 15})
#plt.show()
#plt.savefig("pairwese-rmsd.png") 


