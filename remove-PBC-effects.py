#! Remove PBC effects during MD 
#! __author__ = 'Ropón-Palacios G.' 
#! __date__   = 'November 28, 2021'

#! script correcto, para usar
#! modified March 8, 2022. 

import MDAnalysis as mda 
from MDAnalysis import transformations 
import warnings 
warnings.filterwarnings('ignore')
import optparse
import os


#---- CLI ------

info = """\
++++ Remove PBC Effects into MD Simulations ++++

"""

parser = optparse.OptionParser(description=info) 
parser.add_option("--top", help="topology [TPR or PSF] type=str)", type=str) 
parser.add_option("--itraj", help="traj [XTC, DCD]", type=str)

parser.add_option("--sel", help="selection based en VMD [\"protein\"]", type=str, action='store')


parser.add_option("--rewrap", help="rewrap type [nopbc, nojump]", type=str)
parser.add_option("--type", help="type fitting i.e [solution or membrane]", type=str)
parser.add_option("--prefix_otraj", help="prefix output name [md_0_300]", type=str)

					
options, args = parser.parse_args() 

#----- I/O -----#
itop = options.top
itraj = options.itraj
otraj = options.prefix_otraj 
sel_atoms = options.sel 


if options.rewrap == "nopbc": 
	##! Call gmx, it why gromacs it very fast than MDanalysis transformations wrap or unwrap!.  
	os.system("echo \"1\n0\n\" | gmx trjconv -s " + itop + " -f " + itraj + " -o " + otraj + "_noPBC.xtc " + "-pbc mol -center")


elif option.rewrap == "nojump": 
	os.system("echo \"0\n\" | gmx trjconv -s " + itop + " -f " + itraj + " -o " + otraj + "_nojump.xtc " + "-pbc nojump")
	os.system("echo \"1\n0\n\" | gmx trjconv -s " + itop + " -f " + otraj + "_nojump.xtc" + " -o " + otraj + "_noPBC.xtc " + "-pbc mol -center")

else: 

	print("You not have select rewrap method!") 

## +++++++ Fitting +++++++++++++++++++++++++++++++++ 

if options.type == "solution":
	
	u1 = mda.Universe(itop, otraj + "_noPBC.xtc")
	ref_u1 = u1.copy()
	reference = ref_u1.select_atoms(sel_atoms)
	prot = u1.select_atoms(sel_atoms)

	workflow = transformations.fit_rot_trans(prot, reference)
	u1.trajectory.add_transformations(workflow) 

	system = u1.select_atoms("all") 
	with mda.Writer(otraj + "_fit.xtc", system.n_atoms) as W:
		for ts in u1.trajectory:
			W.write(system)  
	
	print("Fitting rot+trans complete to biomolecule into solution!")


elif options.type == "membrane":
	
	u1 = mda.Universe(itop, otraj + "_noPBC.xtc")
	ref_u1 = u1.copy()
	reference = ref_u1.select_atoms(sel_atoms)
	prot = u1.select_atoms(sel_atoms)

	workflow = transformations.fit_rot_trans(prot, reference, plane='xy', weights="mass") 
	u1.trajectory.add_transformations(workflow) 

	system = u1.select_atoms("all") 
	with mda.Writer(otraj + "_fit.xtc", system.n_atoms) as W:
		for ts in u1.trajectory:
			W.write(system)  
		
	print("Fitting membrane rotxy+transxy complete!")



else: 
	print("Error, not have select any type of transformation!")


print("Done !!!")
