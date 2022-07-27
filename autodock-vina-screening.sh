#!/bin/bash

echo "#--@utor: ROPON-PALACIOS G."
echo "#--date: March 19, 2019"
echo "#--Copyrigth 2019 Umbrella Bioinformatics"
echo "#--E-mail: groponp@gmail.com"
echo "#--This script run a bash loops to perfomed virtual screening"
echo "#--and extract best resuslt of Autodock vina"
echo "#--usage: bash autodock-vina-screening.sh"
echo "¡Waiting ...!"
sleep 5
echo "¡waiting ...  bite more!"
sleep 2
echo "¡I <3 Molecular modeling!"
sleep 10
echo "¡Ready ... processing ligands!"
sleep 3
#--Part0 convert smi to pdbqt more geometry optimization
split -l 1 smi.txt 
for g in x*; do mv "$g" $(echo "$g" | sed 's/^xa/ligando/g'); done 
g=1; for smi in ligando*; do 
	mv "$smi" ligando_$g
	g=$((g+1)); done 
for i in ligando_*; do babel -ismi $i -osdf ${i%}.sdf --gen2D; done 
for j in ligando_*.sdf; do babel -isdf $j -opdb ${j%.sdf}.pdb --gen3D -p 7.4; done 
for l in *.pdb; do babel -ipdb $l -opdbqt ${l%.pdb}.pdbqt; done 
clear 
echo "processing ... geometry optimization using force field MMFF94"
sleep 5 
echo "¡Waiting ... !"
sleep 10
obminimize -ff MMFF94 -n 10000 -sd -c 1e-9 *.pdbqt  
mkdir -p ligands 
mv ligando_*.pdbqt ligands/ 
rm ligando_* 
clear 
echo "Ready to real Docking assay?"
echo "processing ... data for virtual screening"
sleep 5
echo "¡Wainting ... !"
sleep 10 
#--Part1 loop for run autodock-vina 
for file in ligands/*; do 
	tmp=${file%.pdbqt}
	name="${tmp##*/}"
	echo "processing ... $name"
	vina --config config.text --ligand $file --out $name.pdbqt --log $name.log --cpu 2 
	awk '/^[-+]+$/{getline;print FILENAME,$0}' $name.log >> summary; done  
#--Part2 sort best binding energy from autodock-vina 
sort summary -nk 3 > summary_sorted.txt  
rm summary
mkdir -p result
mv *.log ligando_*.pdbqt result/ 
#--Part3 ligand convert from pdbqt to pdb 
cd result/
for m in ligando*.pdbqt; do cut -c-66 $m > $m.pdb; done 
for n in *.pdb; do mv "$n" "${n%.pdbqt.pdb}.pdb"; done 
for p in ligando*.pdbqt; do rm $p; done 
clear
echo "¡Check-in summary_sorted.txt file for select best lingand "
echo "¡The Script a Finished, thanks!"
sleep 5
echo "you can openning ligands and receptor in pymol ...."
sleep 10 


