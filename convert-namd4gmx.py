#! It convert PSF to TOP and PDB to GRO coordinates into charmm format.
#! __author__ = Ropón-Palacios G.
#! __date__   = 27 Jul, 2022.
#! __e-mail__ = groponp@gmail.com 

import parmed as pmd 

## Nota: if you use it output to generate simulation into gromacs with 
## following command : gmx pdb2gmx -f namd.gro -p name.top 
## and you get following error: passed to fgets2 has size 4096. The line starts with: '                    '
## remove all files start with ._name into charmmff download. 

top = pmd.load_file('viro3a_sola_POPC_150mM.psf')
top.save('structure_gmx.top') 

gro = pmd.load_file('viro3a_sola_POPC_150mM.pdb')
gro.save('structure_gmx.gro', format='gro') 


