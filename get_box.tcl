#! Getting center and box for docking with AD4 or Vina
#! __authors__ = 'Ropón-Palacios G and Gustavo Olivos Ramirez'
#! __date__    = 'December 1, 2021.'

#!Note: if x,y or z values are higher than 126, just write 126
mol new receptor.pdb
set sel [atomselect top "protein"]
set geom [measure center $sel]
set minmax [measure minmax $sel]
set boxsize_vina [vecsub [lindex $minmax 1] [lindex $minmax 0]]
set ad4_x [expr [lindex $boxsize_vina 0]/0.375]
set ad4_y [expr [lindex $boxsize_vina 1]/0.375]
set ad4_z [expr [lindex $boxsize_vina 2]/0.375]
set boxsize_autodock [list $ad4_x $ad4_y $ad4_z]

set out_geom [open "box.dat" w]
puts $out_geom "Center is: $geom\nBoxsize for VINA (1 spc): $boxsize_vina\nBoxsize for AUTODOCK (0.375 spc): $boxsize_autodock"
close $out_geom
exit

