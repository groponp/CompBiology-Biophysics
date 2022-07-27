package require topotools
set pdblist [glob *pdb]
set midlist [list ]
foreach pdb $pdblist {
                set mid [mol new $pdb]
                lappend midlist $mid
}
set mol [::TopoTools::mergemols $midlist]
animate write psf merged.psf $mol
animate write pdb merged.pdb $mol
