#! It convert PSF/PDB to CRD coordinates into charmm format.
#! __author__ = Ropón-Palacios G.
#! __date__   = 27 Jul, 2022.
#! __e-mail__ = groponp@gmail.com 

import parmed as pmd 

psf = pmd.load_file('wt_150mM.psf')
psf.save('structure_charmm.psf') 

crd = pmd.load_file('wt_150mM.pdb')
crd.save('structure_charmm.crd', format='charmmcrd') 
