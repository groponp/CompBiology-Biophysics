#! Calculate Surface Accesible Solvent Area (SASA). 
#! _author : Ropón-Palacios G. 
#! _date   : Feb 27, 2022. 
#! _e-mail : groponp@gamil.com 

import optparse 
import subprocess 
 
disclaimer="""<VMD SASA Tool>"""
parser = optparse.OptionParser(description=disclaimer) 
#! INPUTS options
parser.add_option("--coord", help="top/coord [GRO, PSF, PARM7]", type=str)
parser.add_option("--typecoord", help="type [GRO, PSF, PARM7]", type=str)
parser.add_option("--traj", help="traj [XTC, DCD, NETCDF]", type=str)
parser.add_option("--typetraj", help="type [XTC, DCD, NETCDF]", type=str)
parser.add_option("--sel1", help="syntaxis-based in VMD [\"protein and name CA\"]", type=str, action='store')
parser.add_option("--sel2", help="syntaxis-based en VMD [\"protein and name CA\"]", type=str, action='store')

#! Traj options 
parser.add_option("--first", help="firts frame [0]", type=int, default=0)				
parser.add_option("--last", help="last frame [-1]", type=int, default=-1)
parser.add_option("--trajfreq", help="frequency which was save traj ", type=int)
parser.add_option("--dt", help="timestep of integration", type=float, default=0.002)
parser.add_option("--frame1", help="frame 1 only for gromacs [default=0, suggest 1]", type=float, default=0)

#! OUTPUTS options 
parser.add_option("--ofile", help="type output name [prot_nameCA.dat]", type=str)
					
option, args = parser.parse_args() 

sasa="""  
##; Set variables 
set coord       {coor}  
set typecoord   {tcoor}
set traj        {traj}
set typetraj    {ttraj} 
set first       {f} 
set last        {l} 
set selection1  {s1}
set selection2  {s2} 
set outname     {ofile}
set trajfreq    {tf} 
set fs          {0:.3f}
set frame1      {fgro}

##; Load trajectory 
mol new $coord type $typecoord waitfor all 
mol addfile $traj type $typetraj first $first last $last waitfor all 	

##; ! Run SASA calculation

set outfile [ open $outname w ] 
set nf [ molinfo top get numframes ]
set srad 1 
set freq 1 
	
for {set i $frame1} {$i < $nf} {incr i} { 
	set freq $trajfreq 
	set step $fs 
	set t [ expr ($i*$freq*$step)/2 ] 
    set sasa [measure sasa $srad [atomselect top "$select1" frame $i] -restrict [atomselect top "$select2" frame $i]]
    puts $outfile "$t $sasa"
    set i [expr $i + $freq]
}
close $outfile 
exit""".format(coor=option.coord, tcoor=option.typecoord, traj=option.traj, ttraj=option.typetraj, f=option.first, l=option.last,
s1=option.sel1, s2=option.sel2, ofile=option.ofile, tf=option.trajfreq, option.dt, fgro=option.frame1)

##; write procedures 
f1 = open('sasa.tcl', 'w')
f1.write(sasa) 
f1.close() 

subprocess.call("vmd -dispdev text -e sasa.tcl", shell=True) 
print ("SASA analysis was finished !!")






 	       
