import os
import time 
import numpy as np
import pandas as pd
import optparse
from vmd import molecule
from vmd import atomsel

#conda install -c conda-forge vmd-python

disclaimer="""\
This script create a PDB input to perform solute-based concentration Flooding MD.
@Authors: Ropón-Palacios G. & Olivos-Ramirez Gustavo E.
(C) All rigth to Authors. Contact: groponp@gmail.com / gustavo.olivos@upch.pe"""

#Here the inputs are declare
parser = optparse.OptionParser(description=disclaimer)
parser.add_option("--conc",  help="Concentration in mol/L", type=float)
parser.add_option("--boxvol",  help="Box volumen in A^3 [100 A^3 = 1000000]", type=float)
parser.add_option("--pdbi",  help="File name containing the pdb file", type=str)
parser.add_option("--lig",  help="File name containing the ligand pdb file", type=str)
parser.add_option("--maxz",  help="maximum distance in z to place the solute", type=float)

options, args = parser.parse_args()

"""
VMD equations
N_ions = N_A * ionConc(Mol/l) * volume(l)
N_A = 6.022e23
volume(l) = volume_per_water_molecule(l) * nWater
volume_per_water_molecule = 31.05A^3 = 3.105e-26l
N_ions = 6.022e23 * 3.105e-26 * ionConc * nWater
= 0.0187 * ionConc * nWater
"""

concentracion = options.conc
boxvol = options.boxvol
pdbi = options.pdbi
lig = options.lig
maxz = options.maxz

print("Parameters used:\n" "----------------------\n" + 
      "Volume of box is: " + str(boxvol) + "\n"  +
     "Concentration is: " + str(concentracion) + "\n" +
     "PDB name is: " + str(pdbi) + "\n" +
     "Ligand name is: " + str(lig))

time.sleep(3)

navo = 6.022e23 		#number of avogadro
wvol = 31.05 			#water volumen

def n_wat(box):			#to calculate the number of water
    return box/wvol

def nmol(nwat,conc):		#to calculate the number of mols
    total = 0.0187 * conc * nwat
    return int(np.round(total))

#Calling functions

number_water = n_wat(boxvol)
number_mol = nmol(number_water, concentracion)

#Getting distance with vmd-python
molecule.load("pdb", pdbi)

box = atomsel("all")
protein = atomsel("protein")

l_box= np.array(box.minmax())
lz = l_box[1,2]

l_prot = np.array(protein.minmax())

xmin = l_prot[0,0] - 10
xmax = l_prot[1,0] + 10
ymin = l_prot[0,1] - 10
ymax = l_prot[1,1] + 10

lzp = l_prot[1,2] + 5

#Setting packmol Protocol
packmol = """
tolerance 2.0
filetype pdb
output pack%s.pdb
seed -1
avoid_overlap yes
nloop 1000

structure step4_lipid.pdb
    centerofmass
    fixed 0. 0. 0. 0. 0. 0.
end structure

structure %s
    number %s
    resnumbers 3 ##; for numering ligands
    inside box %s %s %s %s %s %s
end structure
""" % (concentracion, lig, number_mol, xmin, ymin, lzp, xmax, ymax, maxz)

f1 = open("pack.inp", "w")
f1.write(packmol)
f1.close()

#Runing packmol
os.system("packmol < pack.inp")

print("FINISHED")
