#!/bin/bash

#######################################################################################
#      –––––– This is a script for convert library ligands of sdf to pdb ––––––	      #
# AUTOR: ROPÓN-PALACIOS G., BSc, MSc(c),					      #
# 		 Bioinformatics & Theoretical Biophysics – Molecular Genetics Lab,    #
#		 Departament of Biology & Physics, Faculty of Sciences,               #
#		 Universidad Nacional de San Antonio Abad del Cusco		      # 
#                Av. La Cultura 733, Wanchaq, Cusco, Peru.		              #               
#		 E-mail: biodano.geo@gmail.com					      #
#######################################################################################

# Run script usign babel module of Openbabel: 
for i in example_*.sdf; do 
	b=`basename $i .sdf`
	echo Processing ligando $b
	mkdir -p $b
	babel -isdf $i -opdb ${b}/out.pdb --gen3D -p 7.4;
done 

# -isdf is a option  for input in sdf format of ligand
# -opdb is a option for output in pdb format 
# --gen3D is a option for generate 3D coordinates
# -p 7.4 is a option for add Hydrogens to pH 7.4

# REMEMBER: 
# If you are usign this script please, cited link of github
# For more information o problems please send me a e-mail.
