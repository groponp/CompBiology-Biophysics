##; make movie from traj 
##; Write by Ropón-Palacios G. 
##; date: June 4, 2021. 


###; Load VMD scene into tcl or vmd extention. 
source block_ion_pass_viro_lig.vmd
display resize 700 700
exec mkdir -p movie 

###; Get info from trajectory and render frame by frame using tachyon 
set nf [molinfo top get numframes]

for {set i 0} { $i < $nf } {incr i} {
	animate goto $i 
        display update on 
	puts  "Rendering Frame #: $i"
	render TachyonInternal movie/frames$i.ppm 

}


###; make the mpeg4
catch { exec ffmpeg -pattern_type glob -r 10 -s 700x700 -i "movie/*.ppm" -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf fps=25 movie/ionblock.mp4}

exit 
