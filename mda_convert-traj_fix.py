i## import library
import MDAnalysis
import optparse
import subprocess


disclaimer="""\
Parallel Molecular Dynamics Engine for Analysis of long trajectories from packages as NAMD, GROMACS, AMBER or DESMOND. 
Autors: Ropón-Palacios G. and Carmen-Sifuentes S.
Department of Physics, Universidade Federal de Alfenas, Minais Gerais, Brasil. 
(C) All rigth to Autors. Contact: groponp@gmail.com"""


#parser = argparse.ArgumentParser(description=disclaimer)
parser = optparse.OptionParser(description=disclaimer)
parser.add_option("--icoord", help="type [GRO, PARM7, CMS or PDB]", type=str)
parser.add_option("--ipsf", help="type [PSF only for NAMD]", type=str)
parser.add_option("--ocoord", help="type [GRO, PSF, PARM7, CMS or PDB]", type=str)
parser.add_option("--itraj", help="traj [XTC, DCD, NETCDF or DRT]", type=str)
parser.add_option("--otraj", help="type [XTC, DCD, NETCDF or DRT]", type=str)
parser.add_option("--sel", help="sele mdanalys-based syntaxis [\"protein or resname LIG\"]", type=str, action='store')    
options, args = parser.parse_args()

##; Inputs 
COORD = options.icoord
TRAJ = options.itraj
PSF = options.ipsf

##; Write coordinates 
u1 = MDAnalysis.Universe(COORD)
COORD1 = u1.select_atoms(options.sel)
COORD1.write(options.ocoord)

##; Write trajectory 
if PSF:
	u2 = MDAnalysis.Universe(PSF, TRAJ)
	system = u2.select_atoms(options.sel)
	with MDAnalysis.Writer(options.otraj, system.n_atoms) as W:
    		for ts in u2.trajectory:
        		W.write(system)
        	print ("Converted frame %d" % ts.frame)
if COORD: 
	u2 = MDAnalysis.Universe(COORD, TRAJ)
	system = u2.select_atoms(options.sel)
	with MDAnalysis.Writer(options.otraj, system.n_atoms) as W:
    		for ts in u2.trajectory:
        		W.write(system)
        	print ("Converted frame %d" % ts.frame)


print (" Coordinate convert is :" , options.ocoord)
print (" Trajectory convert is :" , options.otraj)
