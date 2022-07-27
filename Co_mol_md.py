import optparse
import numpy as np  
import pandas as pd 

disclaimer="""\
Script for calculate Ligand concentration for Flooding MD. @autor: Ropón-Palacios G.  
(C) All rigth to Autors. Contact: groponp@gmail.com"""

parser = optparse.OptionParser(description=disclaimer)
parser.add_option("-c",  help="Concentration in mol/L", type=float) 
parser.add_option("-b",  help="Box volumen in A^3 [100 A^3 = 1000000]", type=float)
parser.add_option("-o",  help="File name containing number of molecules to add", type=str)
      
options, args = parser.parse_args()

##; Constant values 
Litre = 10e27            ##; this is number of molecules H2O in 1 water litre 
avogadro = 6.02 * 10e23  ##; avogadro number 

##; solving equation 
num_mol = (options.b * options.c * avogadro)/Litre 

##; write outfile 
f = open(options.o, 'w')  
f.write("Number of molecules: " + "%s" % str(num_mol) + " molecules" + "\n")
f.write("Concentration used: " + "%s" % str(options.c) + " mol/L" + "\n")
f.write("Box volumen used: " + "%s" % str(options.b) + " A^3" + "\n")    
f.close() 

print("Done, thanks for use my script") 


