;# Script to build complete biomolecule in water, membrane, glycosilated. . 
;# written by : Ropón-Palacios G., Vega-Chozo K., Vilca-Quispe J. 
;# allright : Umbrella Bioinformatics Inc. 
;# date : Dic 27, 2020. 
;# modified : April 9, 2022. 

;# Funcion for orient
;# ========================================================
proc orientMol {molID typeMD pullatoms fixatoms ofile} {
;# * This part the code was taken and modified from source code's QWIKMD plugin v1.3 * 
;# Inputs: 
;# 	- molID     : Is mol identifier after it is loaded into VMD. 
;# 	- typeMD    : This is type of MD to orient protein. values acepted ares MD and SMD. 
;#	- pullatoms : Group's atoms to be fix during SMD. 
;#	- fixatoms  : Group's atoms to be pulling during SMD. 
;#      - ofile     : Output name of file oriented. 

    if {$typeMD == "MD"} {
    	set selall [atomselect $molID all]
    	#set selmove ""
    	# set zaxis {0 0 -1}
    	set center [measure center $selall]
    	set move_dist [transoffset [vecsub {0 0 0} $center]]
    
    	$selall move $move_dist 
    	$selall writepdb $ofile
    
    } elseif {$typeMD == "SMD"} {
        set selall [atomselect $molID all]
        set selanchor [atomselect $molID "$fixatoms"] 
        set selpulling [atomselect $molID "$pullatoms"]
        set anchor [measure center $selanchor]
        set pulling [measure center $selpulling]
        
        ## This align the system in two steps: transvecinv rotates the vector
        ## to be along the x axis, and then transaxis rotates about the y axis to
        ## align your vector with z. 
        ## By Peter Freddolino (http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/6725.html)
        set axis [vecsub $pulling $anchor]
        $selpulling delete
        $selanchor delete
        set M [transvecinv $axis] 
        $selall move $M 
        set M [transaxis y -90] 
        $selall move $M 
        $selall writepdb $ofile

    } else {
    	puts "Error you do not have select nothing typeMD"
    }
    $selall delete
    mol delete all

}

;# Llamando la funcion de orientación
set molID [mol new step1_pdbreader.pdb type pdb waitfor all]           ;# This return molID number, example 0,1,3 .. etc.  
set fixatoms "protein and (segname PROB or segname PROC) and name CA"  ;# This is atoms to be fix during pulling. 
set pullatoms "protein and segname PROA and name CA"                   ;# This is atoms to be pulling alonx Z-axis.  
set typeMD "MD"                                                        ;# When typeMD is set to MD, arguments above are pass but not reads. 
set ofile "proteinFix.pdb"                                             ;# This is output file name of orient protein. 

orientMol $molID $typeMD $pullatoms $fixatoms $ofile
puts "Done"

;# Script for generate topology
;# =================================================================
if {1} {
package require psfgen 

;# Fix topology 
topology ../../toppar/top_all36_prot.rtf
topology ../../toppar/top_all36_na.rtf
topology ../../toppar/toppar_water_ions.str
topology ../../toppar/top_all36_lipid.rtf  ; ## lipids 
topology ../../toppar/top_all36_carb.rtf 
topology ../../toppar/stream/carb/toppar_all36_carb_glycopeptide.str 
}

;# Protein select   
mol new $ofile type pdb waitfor all


set sel [atomselect top "protein"]
set segnames [lsort -unique [$sel get segname]]
foreach segname $segnames {
    puts "Adding protein segname $segname to psfgen"
    #pdbalias residue HIS HSD
    #pdbalias residue GLUP GLU 
    #pdbalias residue CYSP CYS 
    pdbalias atom ILE CD1 CD
    
    set seg ${segname} 
    set sel [atomselect top "protein and segname $segname"]
    $sel set segid $seg 
    
    $sel writepdb tmp.pdb 
    segment $seg {pdb tmp.pdb}
    if {$segname == "PROB"} {patch GLUP PROB:6}  
    
    coordpdb tmp.pdb $seg 
    guesscoord   ;## adiciona, H, y algunos atomos perdidos. 
    
}    

if {1} {
;# Glycan select   
set sel [atomselect top "hetero and not protein"]
set segnames [lsort -unique [$sel get segname]]
foreach segname $segnames {
    puts "Adding glycan segname $segname to psfgen"
    pdbalias residue ANE ANE5AC   ; ## carbohydrates
    pdbalias residue BGL  BGLCNA   ; ## carbohydrates
    #pdbalias residue BGL NGLA 
    pdbalias residue AMA  AMAN   
    pdbalias residue BMA  BMAN 
    pdbalias residue BGA  BGAL 
    pdbalias residue AFU  AFUC
    pdbalias residue AGA  AGALNA
    set seg ${segname}
    set sel [atomselect top "hetero and not protein and segname $segname"]
    $sel set segid $seg
    $sel writepdb tmp.pdb
    segment $seg {pdb tmp.pdb}
    coordpdb tmp.pdb $seg 
    guesscoord   ;## adiciona, H, y algunos atomos perdidos. 
}
}

if {1} { 
##; Patching Disulfide bonds 
patch DISU PROA:336 PROA:361
patch DISU PROA:379 PROA:432 
patch DISU PROA:391 PROA:525 
patch DISU PROA:480 PROA:488 
patch DISU PROB:22  PROB:95 
patch DISU PROB:142 PROB:198 
patch DISU PROC:23  PROC:88 
patch DISU PROC:135 PROC:195 
} 

if {1} {
##; Patching glycans
##; information was get from: step1_pdbreader_glycan.str charmm-gui  
##; only RBD. 
patch NGLB PROA:343 CARA:1 
patch 14BB CARA:1 CARA:2
patch 14BB CARA:2 CARA:3 
patch 13AB CARA:3 CARA:4
patch 12BA CARA:4 CARA:5 
patch 16AT CARA:3 CARA:6 
patch 12BA CARA:6 CARA:7 
patch 16BT CARA:1 CARA:8 
} 

regenerate angles dihedrals 

;# Escribiendo los output
set molofile "cal"     ;# Name of outputs. 
writepsf $molofile.psf ;# Output PSF.
writepdb $molofile.pdb ;# Output PDB.

mol delete all 

;#! Calculate box size, solvate and ionized. 
;# ===============================================

proc cytoplasm {molID boxpad ionconc typeMD smdboxpad dirpull ifile} {
;# * This part the code was taken and modified from source code's QWIKMD plugin v1.3 * 
;# Generate cytoplasm conditions. 
;#         - Calculate Box size 
;#	   - Add Waters to Box
;#	   - Add Ionic atoms to Box at 150 mM (NaCL)
;# Inputs: 
;#         - molID     : Is mol identifier after it is loaded into VMD.    
;#         - boxpad    : Is buffering space add to box 
;#         - ionconc   : Is ionic concentrantion of NaCL into mol/L units.
;#         - typeMD    : This is type of MD to orient protein. values acepted ares MD and SMD. 
;#         - smdboxpad : Is additional spacing into Z-axis to favored the pulling. 
;#         - dirpull   : Is direction of puling into Z-acis, where: 1 +z and -1 -z
;#         - ifile     : This is input name the used into topology build to save psf and pdb. 

    set sel [atomselect $molID "all"]
    set minmax [measure minmax $sel]
    $sel delete
    set xsp [lindex [lindex $minmax 0] 0]
    set ysp [lindex [lindex $minmax 0] 1]
    set zsp [lindex [lindex $minmax 0] 2]

    set xep [lindex [lindex $minmax 1] 0]
    set yep [lindex [lindex $minmax 1] 1]
    set zep [lindex [lindex $minmax 1] 2]
    
    set xp [expr abs($xep - $xsp)]
    set yp [expr abs($yep - $ysp)]
    set zp [expr abs($zep - $zsp)]

    set xsb  ""
    set ysb  ""
    set zsb  ""

    set xeb  ""
    set yeb  ""
    set zeb  ""

    if {$typeMD != "SMD"} {
        set dp [expr sqrt($xp*$xp+$yp*$yp+$zp*$zp)]
        set box_length [expr $dp + 2*$boxpad]
    
        set xsb  [expr $xsp - ($box_length-$xp)/2]
        set ysb  [expr $ysp - ($box_length-$yp)/2]
        set zsb  [expr $zsp - ($box_length-$zp)/2]

        set xeb  [expr $xep + ($box_length-$xp)/2]
        set yeb  [expr $yep + ($box_length-$yp)/2]
        set zeb  [expr $zep + ($box_length-$zp)/2]
    } else {
        set dp [expr sqrt($xp*$xp+$yp*$yp)]
        set box_length [expr $dp + 2*$boxpad]

        set xsb  [expr $xsp - ($box_length-$xp)/2]
        set ysb  [expr $ysp - ($box_length-$yp)/2]
        
	if {$dirpull == -1} {
		set zsb  [expr $zsp - $boxpad - $smdboxpad] 
	} else {
		set zsb  [expr $zsp - $boxpad]
	} 
  
        set xeb  [expr $xep + ($box_length-$xp)/2]
        set yeb  [expr $yep + ($box_length-$yp)/2]
        
	if {$dirpull == -1} {
		set zeb  [expr $zep + $boxpad]

	} else {
		set zeb  [expr $zep + $boxpad + $smdboxpad]

	}
    } 
    

    set boxmin [list $xsb $ysb $zsb]
    set boxmax [list $xeb $yeb $zeb]

    #! call solvate 
    package require solvate
    solvate ${ifile}.psf ${ifile}.pdb -minmax [list $boxmin $boxmax] -o ${ifile}_solvated -s W
    
    #! call Autoionize 
    package require autoionize 
    autoionize -psf ${ifile}_solvated.psf -pdb ${ifile}_solvated.pdb -cation SOD -anion CLA -sc $ionconc -o ${ifile}_150mM -seg ION 

    set xfile [open "cellbox.inp" w]
    
    puts $xfile "centerX : [expr [expr $xsb + $xeb] /2]" 
    puts $xfile "centerY : [expr [expr $ysb + $yeb] /2]" 
    puts $xfile "centerZ : [expr [expr $zsb + $zeb] /2]" 

    puts $xfile " cB1x :  [expr abs($xeb - $xsb)]"
    puts $xfile " cB2y :  [expr abs($yeb - $ysb)]"
    puts $xfile " cB3z :  [expr abs($zeb - $zsb)]" 

    close $xfile

    #set center [list [format %.2f $centerX] [format %.2f $centerY] [format %.2f $centerZ]]
    #set length [list [format %.2f $cB1] [format %.2f $cB2] [format %.2f $cB3]]
    #set QWIKMD::cellDim [list $boxmin $boxmax $center $length]

    mol delete $molID

} 
 
;# Call function cytoplasm 
set molID [mol new $molofile.psf type psf waitfor all]
mol addfile $molofile.pdb type pdb waitfor all

set boxpad 0       ;# If you want build minimal box, set this value to 0, and if you want add padding set to values between 10-15.  
set ionconc 0.150  ;# Ionic concentration ino mol/L units. 
set typeMD  "MD"   ;# In MD type, dirpull and smdboxpad not are used, but need be pass. 
set smdboxpad 35   ;# Aditional box spacing for perform smd simulation.
set dirpull -1      ;# This is direction of puling in Z-axis, where: 1 => +z and -1 => -z.
set ifile $molofile  

cytoplasm $molID $boxpad $ionconc $typeMD $smdboxpad $dirpull $ifile

;#! call Hyrogen Mass Repartition (HMR) if you need! 
;# ====================================================
set hmr 1   ;# 1 => yes and 0 => no 
set oS  "macOS" ;# MacOS or Linux

if {$hmr} {

if { $oS == "macOS"} {
	set vmd_exec "/Applications/VMD\ 1.9.4a48-Catalina-Rev7.app/Contents/MacOS/startup.command" 
	exec $vmd_exec -dispdev text -e do_hmr.tcl -args ${molofile}_150mM.psf ${molofile}_150mM.pdb 
	
} elseif {$oS == "Linux"} {
	exec vmd -dispdev text -e do_hmr.tcl -args ${molofile}_150mM.psf ${molofile}_150mM.pdb 
} else {
	puts "Error not operative system given"
}

} 

quit

