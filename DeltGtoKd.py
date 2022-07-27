#!/usr/env python 

"""
@autor: ROPON-PALACIOS G.
date: March 11, 2019.
Copyright 2019 Umbrella Bioinformatics.
e-mail: groponp@gmail.com
usage: python DeltGtoKd.py --binding_energy 12.7 --conc nM -out_name outname.txt 
"""
from __future__ import print_function
import argparse as arg
import os
import sys
import math

#--Defining arguments
if __name__ == "__main__": 
	ap = arg.ArgumentParser(description=__doc__) 
	io = ap.add_argument_group('Input options')
	io.add_argument('-be','--binding_energy', required=True,
                help='Binding energy from Docking in Kcal/mol')

	out = ap.add_argument_group('Output options')
	out.add_argument('--conc', type=str, required=True, 
                 help='Concentration of Kd in M, mM or nM')
	out.add_argument('-o','--out_name', type=str, 
                 help='Name of your output')
	cmd = ap.parse_args()  
#-- Math Algorithm 
delg = float(cmd.binding_energy) 
if cmd.conc == 'M': 
	Kd = math.exp((delg*1000)/(1.98*298.15)) 
	conc = 'M'
elif cmd.conc == 'mM':
	Kd = (math.exp((delg*1000)/(1.98*298.15)))*1000000 
	conc = 'mM'
elif cmd.conc == 'nM': 
	Kd = (math.exp((delg*1000)/(1.98*298.15)))*1000000000 
	conc = 'nM'
else:
	print ('unsopported operation')
#--Print Results 
#print("Kd =", Kd, conc)
#--Defining output 	
if cmd.out_name: 	
	filename = cmd.out_name 
	f = open(filename, "a")
	f.write(str(Kd) + str(conc) +  '\n')








