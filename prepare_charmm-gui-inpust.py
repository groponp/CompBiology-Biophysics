###; An PyMOL script into python to create inputs para charmm-gui 
###; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###; @utor: Ropón-Palacios G.
###; date: September 21. 2021. 
###; e-mail: groponp@gamil.com / ggdrpalacios@uesc.br  

###; How run it?:  pymol -cq prepare_charmm-gui_inputs.py 
from pymol import cmd 
import os 

#; Load original PDB to fix
cmd.load('../CHEMBL685__input3__variant1r_complex_6UL9.pdb')
cmd.remove('hydro')  ##; remove hydrogens. 
cmd.select('het')
cmd.alter('sele', 'chain=\"Z\"')
cmd.h_add('chain Z')
cmd.delete('sele')

#; Generate inputs to Charmm-gui

lig_name = "lig.mol2"
prot_name = "complex.pdb"

cmd.save(lig_name, selection="chain Z")
cmd.save(prot_name, selection="all")

#; Fix name ligand+ and complex
#; Note: name of ligand are 4 last letter from chembl. 
#os.system('obabel tmp.pdb -O lig.mol2 -h')

#; merge lig/prot 
#cmd.delete("all")
#cmd.load('prot.pdb')
#cmd.load('lig.mol2')
#cmd.select('het')
#cmd.alter('sele', 'chain=\"Z\"')
#cmd.save('complex.pdb', selection="all")

#; Fix name into files
os.system('sed \'s/UNL/LIG/g\' complex.pdb > complex_fix.pdb')
os.system('sed \'s/UNL1/LIG/g\' lig.mol2 > tmp.mol2')
os.system('sed \'s/UNL/LIG/g\' tmp.mol2 > tmp1.mol2')
os.system('sed \'3 s/CH.*/LIG/g\' tmp1.mol2 > lig_fix.mol2') 
os.system('rm tmp*.mol2')

#; Generate PQR file into amber format to check protonation state. 
#; Note: check what protein not have hydrogens, beacause generate problems with pdb2pqr. 

os.system('pdb2pqr30 complex_fix.pdb complex_fix.pqr --ff=AMBER --keep-chain --ffout=AMBER --drop-water --include-header --titration-state-method=propka --with-ph 7.0 --pH 7.0 2>/dev/null')

print('Done!!')

