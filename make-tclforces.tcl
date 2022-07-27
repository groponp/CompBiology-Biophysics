;## Generate, mark into B-colunm to process for TCLForces or TCLBC 
;## Write by Ropón-Palacios G. 
;## date: March 25, 2021. 
;## e-mail: groponp@gmail.com 

mol new input.pdb type pdb waitfor all 
# mol addfile input.psf type psf waitfor all 

set targetpdb    "marks_lig_helixs_ca.pdb" ;# PDB for calling into TCLForces into config file

##; make marks 
set ligmark   "resname TYL and noh"   ;# Ligand pulling no using hydrogen 
set hexlimark "protein and helix"     ;# fixing helix structure 

set mark1 "1.00"  ;# for ligand, pulling 
set mark2 "2.00"  ;# for helix, for fixing 

set all [atomselect top all] 
$all set beta 0 
$all set occupancy 0 

set lig [atomselect top $ligmark] 
set ligmass [$lig get mass] 
$lig set beta $mark1
$lig set occupancy $ligmass 

set helix [atomselect top $helixmark] 
set helixmass [$helix get mass] 
$helix set beta $mark2 
$helix set occupancy $helixmass 

$all writepdb $targetpdb 
quit 
