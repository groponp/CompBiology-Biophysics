;## TclForces script for run smd to calculate Jarzynski equality 
;## Write by Ropón-Palacios G. 
;## date: March 25, 2021
;## e-mail: groponp@gmail.com 
;## reference: Acta Mech. Sin. (2012) 28(3):891–903. DOI 10.1007/s10409-012-0112-9. 

##; setting marks into PDB file 
set mark1 1.0 
set mark2 2.0

##; looping through set marks 
set ligtargets   {}
set ligmasses    {}
set helixtargets {}
set helixmasses  {}

set inStream [open $targetAtomPdb r] 
foreach line [split [read $inStream] \n] {
	set type [string trim [string range $line 0 5]]
	set name [string trim [string range $line 12 15]] 
	set resid [string trim [string range $line 22 25]]
	set beta [string trim [string range $line 60 65]]
	set occupancy [string trim [string range $line 54 59]]
	set segname [string trim [string range $line 72 75]] 
	
	if { ($type eq "ATOM" || $type eq "HETATM") && $beta == $mark1 } { 
		lappend ligtargets "$segname $resid $name" 
		lappend ligmasses  $occupancy 
	   } 

	else { ($type eq "ATOM" || $type eq "HETATM") && $beta == $mark2 } { 
		lappend helixtargets "$segname $resid $name" 
		lappend helixmasses  $occupancy 
	  }
} 

close $inStream 

##; Getting atoms index for atoms used in tclforces calculation 
set ligatoms   {}
foreach target1 $ligtargets { 
	lassing $target1 segname resid atom 
	set atomindex [atomid $segname $resid $atom] 
	lappend ligatoms $atomindex 
	addtom $atomindex 
} 

set helixatoms {} 
foreach target2 $helixtargets { 
	lassing $target2 segname resid atom 
	set atomindex [atomid $segname $resid $atom] 
	lappend helixatoms $atomindex 
	addtom $atomindex 
} 

##; adding atomindex to groups 
set ligand [addgroup $ligatoms] 
set helix  [addgroup $helixatoms] 


##; set output frequency, initialize the time counter 
set Tclfreq 50 
set t 0 

##; Constraint point 
set c1x	0.0  ;## center of mass of molecule to fix. not add values. 
set c1y 0.0 
set c1z 0.0 

set c2x 0.0 ;## pulling molecule center of masss. 
set c2y 0.0 
set c2z 0.0 ;## add center of mass in Z, if you want pulling into Z. 

##; force constant (kcal/mol/A^2) 
set k 5.0 

##; pulling velocity (A/timestep) 
set v 0.00002  ;# pulling to an velocity of: 2 A/timespre whean time integrations is 0.002. 

set outfilename force_pulling.dat 
open $outfilename w 

##; procedure to tclforces sintaxys 
proc calcforces {} { 
	global Tclfreq t k v ligand helix c1x c1y c1z c2x c2y c2z outfilename 
	
	;## load coordinates 
	loadcoords coordinate 
	
	set r1 $coordinate($helix) 
	set r1x [lindex $r1 0] 
	set r1y [lindex $r1 1] 
	set r1z [lindex $r1 2] 

	set r2 $coorinate($ligand) 
	set r2x [lindex $r2 0] 
	set r2y [lindex $r2 1] 
	set r2z [lindex $r2 2] 

	;## calculate forces 
	
	set f1x [ expr $k*($c1x-$r1x)] 
	set f1y [ expr $k*($c1y-$r1y)] 
	set f1z [ expr $k*($c1z-$r1z)]
	lappend f1 $f1x $f1y $f1z  

	##; /// equation is F = k(vt - (r2z-c2z)  
	##; /// Where k : is spring constat; v: velocity pulling t: is time simulation
	##; /// x1 : position of inhibitor or molecule for pulling. 

	#set f2x [ expr $k*($c2x-$r2x)] 
	#set f2y [ expr $k*($c2x-$r2x)] 
	set f2z [ expr $k*($v*$t-($r2z-$c2z)] ; ##; Our applied force to Z axis. 
	#lappend f2 $f2x $f2y $f2z 
	lappend f2 $f2z 

	;## apply forces 	
	addforce $helix  $f1 
	addforce $ligand $f2

	;## output 
	
	set foo [expr $t % $tclfreq] 
	if { $foo == 0} {
		set outfile [open $outfilename a] 
		set time [open $t*2/1000.0] 
 		##; /// delta-Work is calculate using DelW = f2z*($rz2-$c2z) 
		##; Where deltaW : is change of work into Z axis. 
		#set Workz [expr (f2z*($r2z-$c2z))/69.479] ;## divide 69.479 pN/A^2 to conver a kcal/mol
		set dt 0.002 
		set Workz [expr ($f2z*$v*$dt)/69.479]  ;## divide 69.479 pN/A^2 to conver a kcal/mol 
		puts $outfile "$time\t$r2z\t$f2z\t$Workz"
		close $outfile 
	} 
	incr t 
	return 
}


  









