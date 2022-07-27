###; Calculate Spring constant K for applied harmonic potential into Umbrella sampling 
###; Write : Ropón-Palacios G. 
###; date : April 23, 2021. 
###; Source from : http://www.strodel.info/index_files/lecture/html/US-MD.html 
###; Update to MDAnalysis it script: use follows link --> https://www.mdanalysis.org/MDAnalysisTutorial/atomgroups.html 
###; update September 19, 2021.

import sys
import pandas as pd
import numpy as np
import MDAnalysis as mda 

####; Argparse functions 
####; ~~~~~~~~~~~~~~~~~~~~~~~~~
itop  = sys.argv[1]
itraj = sys.argv[2]

####; Load trajectory into mda universe and measure \
####: center of mass 4 glycerol mols!
####; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
u1 = mda.Universe(itop, itraj) 

glp = u1.select_atoms("segid GLP and noh")
glx = u1.select_atoms("segid GLX and noh")
gly = u1.select_atoms("segid GLY and noh")
glz = u1.select_atoms("segid GLZ and noh")


####; Calculate distance to each glycerol mol
####; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rows = len(u1.trajectory)
cols = 4  #; num molecules

dist_matrix_mols = np.zeros(rows,cols)  #; (rows, columns)

for frame in u1.trajectory:

"""

Here using only Z axis, for restraint collective variable bias, into other cases into
gromacs for example an np.lingal would be used to calcualte norm to (X,Y,Z) axis distance.
More rigth approach is applied harmonic force to reaction coordiante driven the bias.  

Run a MD simulation without restraint for 5 ns, for example!!!!. 

"""

	if frame == 0: 
		comp = glp.center_of_mass()    ##; make reference to glycerol id name, into all cases. not confuse with (x,y,z) axis.
		comx = glx.center_of_mass() 
		comy = gly.center_of_mass() 
		comz = glz.center_of_mass() 
		
		dist_matrix_mols[frame][0] += 0   ##; it it beacause it difference into firt frame is 0. wyt com1 - com1 = 0 ; (15-15)
		dist_matrix_mols[frame][1] += 0
		dist_matrix_mols[frame][2] += 0
		dist_matrix_mols[frame][3] += 0

 	else: 

	glp_com = glp.center_of_mass() 
	glx_com = glx.center_of_mass() 
	gly_com = gly.center_of_mass() 
	glz_com = glz.center_of_mass() 

	dist_matrix_mols[frame][0] += np.absolute(glp_com - comp) + distance_matrix_mols[frame-1][0]  ##; it lat sum it to accumulate distance.
	dist_matrix_mols[frame][1] += np.absolute(glx_com[2] - comx[2]) + distance_matrix_mols[frame-1][1]
	dist_matrix_mols[frame][2] += np.absolute(glx_com[2] - comy[2]) + distance_matrix_mols[frame-1][2]
	dist_matrix_mols[frame][3] += np.absolute(glz_com[2] - comz[2]) + distance_matrix_mols[frame-1][3]



#####; Calculate harmonic bias, to be applied each windows 
#####; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""" calculate K using following eq1 
     kus = (kb*T*Na)/variance^2
     kb --> boltzamn constat 
     T ---> temperature into K 
     Na ---> Avogadro number 

""" 

##; constants 
kb = 1.380649E-23   #; in  J/K
T  = 298            #; K  from simualtion temperature 
Na = 6.02214076E23  #; mol^-1 

##; Solving equation
com_std_z_axis = np.std(distance_matrix_mol, axis=0) #; in Angstrom, calculate STD by each colum. 
kus = (kb*T*Na)/np.power(com_std_z_axis,2)
kus_kJ = kus/1000 

kus_kcal_molAngstrom = (kus_kJ*0.239006)/1

print ( " Standard deviation: " + str(com_std_z_axis))
print ( "Harmonic potential is: " + str(kus_kJ) + " kJ/mol*A^2")
print ("Harmonic potential is: " + str(kus_kcal_molAngstrom) + " kcal/mol*A^2 ")
