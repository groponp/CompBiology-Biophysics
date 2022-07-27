! Create quality Movie for MD traj
! __author__  = Ropón-Palacios G.
! __date__    = 14 Jul, 2021. 
! __e-mail__  = groponp@gmail.com 

$ global fps=15 render=t draft=f name=movieHD restart=t
$ scene0 visualization=example/scene1.tcl resolution=600,600  ambient_occlusion=t
!$ scene1 visualization=example/scene1.tcl resolution=600,600 after=scene0
!$ scene2 visualization=example/scene2.tcl resolution=600,600 after=scene1
!$ scene3 visualization=example/scene3.tcl resolution=600,600 after=scene2


# scene0
do_nothing         t=1s
rotate             axis=z angle=180 t=4s sigmoid=sls  
zoom_in            scale=1.2


! Generate second scene with surface protein 
highlight  selection='protein' style=surf color=red t=2s
{rotate  axis=y angle=720 t=2s sigmoid=sls fraction=:0.25; zoom_in scale=1.5}
rotate  axis=y angle=720 t=2s sigmoid=sls fraction=0.25:0.5

highlight  selection='chain A' style=surf color=red t=2s mode=u 
rotate  axis=y angle=720 t=2s sigmoid=sls fraction=0.5:0.75
{rotate  axis=y angle=720 t=2s sigmoid=sls fraction=0.75:; zoom_out scale=1.5}





!# scene1
!zoom_in             scale=1.2

!# scene2
!{rotate  axis=y angle=720 t=2s sigmoid=sls fraction=:0.25; zoom_in scale=1.5}
!rotate  axis=y angle=720 t=2s sigmoid=sls fraction=0.25:0.5

!# scene3 
!rotate  axis=y angle=720 t=2s sigmoid=sls fraction=0.5:0.75
!{rotate  axis=y angle=720 t=2s sigmoid=sls fraction=0.75:; zoom_out scale=1.5}
