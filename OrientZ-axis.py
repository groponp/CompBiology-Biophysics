#! Orient to Z axis 
#! __author__  = 'Ropon-Palacios G.' 
#! __date__    = 'March 17, 2021.'
#! __version__ = '1.0'

import os
import optparse  

header_code = """+++ Orient to Z-axis by Ropón-Palacios G. +++ """

parser = optparse.OptionParser(description=header_code) 
parser.add_option("--pdb", help="coord [PDB]", type=str) 
parser.add_option("--ofile_orient", help="name to output file [proteinOrintZ]", type=str)

options, args = parser.parse_args()

#! inputs
pdb_charmmgui = options.pdb
ofile = options.ofile_orient

#! TCL scripts 
tcl_orient="""
package require Orient
namespace import Orient::orient

set pdb %s
set ofile %s
set seltext "all"

mol load pdb $pdb 

proc orient_z { seltext } { 
	;#! orient in z axis 
	
	set sel [atomselect top $seltext]
	set I [draw principalaxes $sel]

	set A [orient $sel [lindex $I 2] {0 0 1}]
	$sel move $A
	set I [draw principalaxes $sel]

	set A [orient $sel [lindex $I 2] {0 0 1}]
	$sel move $A
	set I [draw principalaxes $sel]

	#set all [atomselect top "all"]
	#$all writepdb ${ofile}.pdb
}

#! run routine
[orient_z $seltext]

#! add chains 
set all [atomselect top $seltext]

set segnames [lsort -unique [$all get segname]]
foreach seg $segnames {
        set selseg [atomselect top "segname $seg"]
        $selseg set chain [string index $seg 3]
}

$all writepdb $ofile.pdb


quit
""" % (pdb_charmmgui, ofile)

#! WRITE scripts 
f = open('orientZ-axis.tcl', 'w')
f.write(tcl_orient) 
f.close() 

#! RUN scripts 
import os
os.system("vmd -dispdev text -e orientAxis.tcl ")
print("Protein was Oriented to Z-axis") 


