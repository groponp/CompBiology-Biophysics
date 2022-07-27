#!/bin/bash

#######################################################################################
#      –––––– This is a script for run automatically to virtual screening ––––––	  #
#   Written By: ROPÓN-PALACIOS G., BSc, MSc(c),						                  #
# 		 Bioinformatics & Thereotical Biophysics – Molecular Genetics Lab,            #
#		 Departament of Biology & Physics, Faculty of Sciences,                       #
#		 Universidad Nacional de San Antonio Abad del Cusco							  # 
#        Av. La Cultura 733, Wanchaq, Cusco, Peru.									  #               
#		 E-mail: biodano.geo@gmail.com												  #
#######################################################################################

# Run Script using vina software:
for filename in ligando_*.pdbqt; do  #ligando name for spanish users that can be changed to ligand for english users. 
    b=`basename $filename .pdbqt`
    echo Processing ligando $b
    mkdir -p $b
    ./vina --config config.text --ligand $filename --out ${b}/ligand_output.pdbqt --log ${b}/log.text
done

# ./vina is a command for run software 
# --config is a command for run "config.text", which have all paremeters for run vina
# --ligand is a command for enter ligand in pdbqt format
# --out is a command that for output results in pdbqt format
# --log is a command for do a file "log.text" which have all scoring in kcal/mol
# for filename in ligando_*.pdbqt is a loop in bash scripting that run docking for each ligand whre * is {1;2;3...n} ligands numbers
# b=`basename $filename .pdbqt` is a variable name for make working directory for each ligand
# echo Processing ligando $b is a function for print "Processing ligando" string when run the script for each ligand 
# mdkir -p $b is a function for make working directry for each ligand, where -p indicate make directory if this not exist

# REMENBER:
# This script is a very small modification and description of original script from autodock vina manual
# The original script can be find in autodock vina manual page: http://vina.scripps.edu/manual.html
# If you are using this script please, cited link of github
# For more information or problems please send me a e-mail.


