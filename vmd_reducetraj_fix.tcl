mol new ../01_system/viro3a_lig3925_3_POPC_150mM.psf type psf waitfor all
mol addfile eq1-out.dcd type dcd first 0 last -1 

animate write dcd reduce_1000f.dcd beg 0 end -1 skip 100 waitfor all top
quit  
