#! _author : Ropón-Palacios G. 
#! _date   : March 9, 2022. 
#! _e-mail : groponp@gamil.com 

#! Making Data analysis. 

for file in C*; do 
	cd $file
	echo "Setting Domain to $file"
	sleep 2
	clear
	vmd -dispdev text -e ../scripts/set_domain.tcl 
	cd ../
done; 


for file in C*; do 
	cd $file 
	echo "Running RMSD analysis to $file"
	sleep 2 
	clear
	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=rmsd_prot.dat 
	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=rmsd_kin.dat
	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=rmsd_roc.dat
	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=rmsd_cor.dat

	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname LIG and not name H" --ofile=rmsd_LIG.dat 
	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname ATP and not name H" --ofile=rmsd_ATP.dat 
	python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname GDP and not name H" --ofile=rmsd_GDP.dat 

#	if [[ $file == "Control" ]] 
#	then
#		python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=rmsd_prot.dat 
#		python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=rmsd_kin.dat
#		python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=rmsd_roc.dat
#		python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=rmsd_cor.dat
#		
#		python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname ATP and not name H" --ofile=rmsd_ATP.dat 
#		python ../scripts/mda_rmsd_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname GDP and not name H" --ofile=rmsd_GDP.dat 
#	fi 

	cd ../
done; 

for file in C*; do 
	cd $file 
	echo "Running RMSF analysis to $file"
	sleep 2
	clear  
	python ../scripts/mda_rmsf_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=rmsf_prot.dat 
	python ../scripts/mda_rmsf_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=rmsf_kin.dat 
	python ../scripts/mda_rmsf_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=rmsf_roc.dat 
	python ../scripts/mda_rmsf_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=rmsf_cor.dat 	
	cd ../
done; 

for file in C*; do 
	cd $file 
	echo "Running 2Dmatrix itself to $file"
	sleep 2
	clear  
	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=2Dmatrix_prot.svg 
	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=2Dmatrix_kin.svg 
	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=2Dmatrix_roc.svg 
	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=2Dmatrix_cor.svg 

	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname LIG and not name H" --ofile=2Dmatrix_LIG.svg 
	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname ATP and not name H" --ofile=2Dmatrix_ATP.svg 
	python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname GDP and not name H" --ofile=2Dmatrix_GDP.svg 

#	if [[ $file == "Control" ]]
#	then
#		
#		python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=2Dmatrix_prot.svg 
#		python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=2Dmatrix_kin.svg 
#		python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=2Dmatrix_roc.svg 
#		python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=2Dmatrix_cor.svg 
#	
#		python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname ATP and not name H" --ofile=2Dmatrix_ATP.svg 
#		python ../scripts/mda_2Dmatrix_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname GDP and not name H" --ofile=2Dmatrix_GDP.svg 	
#	fi 
	cd ../

done; 

for file in C*; do 
	cd $file 
	echo "Running Rgyr to $file"
	sleep 2
	clear  
	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=rgyr_prot.dat 
	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=rgyr_kin.dat 
	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=rgyr_roc.dat 
	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=rgyr_cor.dat 

	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname LIG and not name H" --ofile=rgyr_LIG.dat 
	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname ATP and not name H" --ofile=rgyr_ATP.dat 
	python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname GDP and not name H" --ofile=rgyr_GDP.dat 

#	if [[ $file == "Control" ]]
#	then 
#		
#		python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="protein and name CA" --ofile=rgyr_prot.dat 
#		python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid KIN and name CA" --ofile=rgyr_kin.dat 
#		python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid ROC and name CA" --ofile=rgyr_roc.dat 
#		python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="segid COR and name CA" --ofile=rgyr_cor.dat 
#
#		python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname ATP and not name H" --ofile=rgyr_ATP.dat 
#		python ../scripts/mda_rgyr_fix.py --coord=system_fix_domain.pdb --traj=md_0_200_noPBC.xtc --sel="resname GDP and not name H" --ofile=rgyr_GDP.dat 
#	fi 
	cd ../ 
done; 

for file in C*; do 
	cd $file 
	echo "Running SASA to $file"
	sleep 2
	clear  
	python ../scripts/vmd_sasa_fix.py --coord=system_fix_domain.pdb --typecoord=pdb --traj=md_0_200_noPBC.xtc --typetraj=xtc --sel1="protein and same residue as within 5 of resname LIG" --sel2="protein and same residue as within 5 of resname LIG" --first=0 --last=-1 --trajfreq=400 --dt=0.002 --ofile=sasa_prot_LIG.dat
	python ../scripts/vmd_sasa_fix.py --coord=system_fix_domain.pdb --typecoord=pdb --traj=md_0_200_noPBC.xtc --typetraj=xtc --sel1="protein and same residue as within 5 of resname ATP" --sel2="protein and same residue as within 5 of resname ATP" --first=0 --last=-1 --trajfreq=400 --dt=0.002 --ofile=sasa_prot_ATP.dat
	python ../scripts/vmd_sasa_fix.py --coord=system_fix_domain.pdb --typecoord=pdb --traj=md_0_200_noPBC.xtc --typetraj=xtc --sel1="protein and same residue as within 5 of resname GDP" --sel2="protein and same residue as within 5 of resname GDP" --first=0 --last=-1 --trajfreq=400 --dt=0.002 --ofile=sasa_prot_GDP.dat
	
#	if [[ $file == "Control" ]]
#	then 
#		python ../scripts/vmd_sasa_fix.py --coord=system_fix_domain.pdb --typecoord=pdb --traj=md_0_200_noPBC.xtc --typetraj=xtc --sel1="protein and same residue as within 5 of resname ATP" --sel2="protein and same residue as within 5 of resname ATP" --first=0 --last=-1 --trajfreq=400 --dt=0.002  --ofile=sasa_prot_ATP.dat
#		python ../scripts/vmd_sasa_fix.py --coord=system_fix_domain.pdb --typecoord=pdb --traj=md_0_200_noPBC.xtc --typetraj=xtc --sel1="protein and same residue as within 5 of resname GDP" --sel2="protein and same residue as within 5 of resname GDP" --first=0 --last=-1 --trajfreq=400 --dt=0.002  --ofile=sasa_prot_GDP.dat
#	fi 
	cd ../
done;

for file in C*; do 
	cd $file 
	echo "Running FEL to $file"
	sleep 2 
	clear 
	python ../scripts/FEL.py --file1=rgyr_LIG.dat --file2=rmsd_LIG.dat --temperature=310.0 --bin=25 --ofile=FEL_LIG.svg
	python ../scripts/FEL.py --file1=rgyr_ATP.dat --file2=rmsd_ATP.dat --temperature=310.0 --bin=25 --ofile=FEL_ATP.svg
	python ../scripts/FEL.py --file1=rgyr_GDP.dat --file2=rmsd_GDP.dat --temperature=310.0 --bin=25 --ofile=FEL_GDP.svg

#	if [[ $file == "Control" ]]
#	then
#		python ../scripts/FEL.py --file1=rgyr_ATP.dat --file2=rmsd_ATP.dat --temperature=310.0 --bin=25 --ofile=FEL_ATP.svg
#		python ../scripts/FEL.py --file1=rgyr_GDP.dat --file2=rmsd_GDP.dat --temperature=310.0 --bin=25 --ofile=FEL_GDP.svg
#	fi 
	cd ../
done; 

echo "Done!!"	
	
