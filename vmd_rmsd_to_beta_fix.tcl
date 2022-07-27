# This script contains a procedure called rmsd_residue_over_time that calculates the average RMSD for each residue in a selection over all frames in a trajectory. The procedure is called as:
# rmsd_residue_over_time mol sel_resid
#where mol is the molecule in VMD and sel_resid is a list of the residue numbers in that selection.

#You can use the procedure for any residue or list of residues. Here, as an example, we will make a selection for all residues in the protein. Note that this will take a long time to calculate:
mol new md1_renamed.pdb  type pdb waitfor all
mol addfile traj_1000.xtc type xtc first 0 last -1 waitfor all


set sel_resid [[atomselect top "protein and name CA"] get resid]

#The procedure is presented below.  It also sets the B value to the value calculated, so you can color the protein by RMSD.  The call for the procedure is at the end of the file.

proc rmsd_residue_over_time {{mol top} res} {

    # use frame 0 for the reference
    set reference [atomselect $mol "protein" frame 0]
    # the frame being compared
    set compare [atomselect $mol "protein"]
    #make a selection with all atoms
    set all [atomselect top all]
    #get the number of frames
    set num_steps [molinfo $mol get numframes]
    
    foreach r $res {
	set rmsd($r) 0
    }
    
    #loop over all frames in the trajectory
    for {set frame 0} {$frame < $num_steps} {incr frame} {
	puts "Calculating rmsd for frame $frame ..."
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$all move $trans_mat
	# compute the RMSD
	
	#loop through all residues
	foreach r $res {
	    set ref [atomselect $mol "segid PROA and resid $r and noh" frame 0]
	    set comp [atomselect $mol "segid PROA and resid $r and noh" frame $frame]
	    set rmsd($r) [expr $rmsd($r) + [measure rmsd $comp $ref]]
	    $comp delete
	    $ref delete
	}
    }
    set ave 0
	foreach r $res {
	    set rmsd($r) [expr $rmsd($r)/$num_steps]
	    # print the RMSD
	    puts "RMSD of residue $r is $rmsd($r)"
	    set res_b [atomselect $mol "resid $r"] 
            $res_b set beta $rmsd($r)
            $res_b delete
	    set ave [expr $ave + $rmsd($r)]
	    
	}
    set ave [expr $ave/[llength $res]]
    puts " Average rmsd per residue:   $ave"
    
}

#Call the procedure

rmsd_residue_over_time top $sel_resid
