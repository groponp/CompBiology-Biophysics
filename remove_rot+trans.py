##; Script to remove rot+trans into long trajectory
##; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##; @autor: Ropón-Palacios G. 
##; date: September 26, 2021. 

import os
import optparse 

info = """
It is an in-house script to remove rot+translation into DCD traj
@autor: Ropón-Palacios G.
e-mail: ggdrpalacios@uesc.br"""
 
disclaimer= info 				

#parser = argparse.ArgumentParser(description=disclaimer)
parser = optparse.OptionParser(description=disclaimer) 
parser.add_option("--top", help="coord [PSF]", type=str)
parser.add_option("--traj", help="traj [DCD]", type=str)
parser.add_option("--sel", help="selection based en VMD [\"protein or (resname LIG)\"]", type=str,
					 action='store')

option, args = parser.parse_args() 

tcl_script="""
package require pbctools 

#! Load traj
mol new %s type pdb waitfor all
mol addfile %s type xtc waitfor all

#! selection 
set sel1 "%s"


#! PBC wrap trj; it's necessary for unwrapped traj
#! without wrapAll on, option on NAMD config. 
pbc wrap -all -compound res -center bb -centersel "protein" 


#! ++++++++++ Not use ++++++++++++++++
#set sel1 [atomselect top \"$sel1\"]
#set com_fix [measure center $sel1]

#! Remove translation
#set nf [ molinfo top get numframes ]

#for {set j 0} {$j < $nf} {incr j} { 
#$sel1 frame $j
#set com_drag [measure center $sel1]
#set diff_com [vecsub $com_drag $com_fix]
#set move [vecscale -1 $diff_com]
#$sel1 moveby $move 
#}
#! +++++++++ Not use ++++++++++++++++++++

#! Remove rotation
#! Note: It part is inspirate into work of Jim Phillips (jim@ks.uiuc.edu)
set nf [ molinfo top get numframes ]
set ref [atomselect top \"$sel1\" frame 0]
set sel2 [atomselect top \"$sel1\"]
set all [atomselect top "all"]

for {set i 2} {$i < $nf} {incr i} {
	$sel2 frame $i
	$all frame $i 
	$all move [measure fit $sel2 $ref]

}

#! Write file 
animate write dcd prod_w_rot+trans.dcd beg 0 end -1 skip 1 waitfor all top 

quit
""" % (option.top, option.traj, option.sel)

##; write procedures 
f = open('remove_rot+trans.tcl', 'w')
f.write(tcl_script) 
f.close() 

os.system("vmd4 -dispdev text -e remove_rot+trans.tcl 2>/dev/null")
 

				

