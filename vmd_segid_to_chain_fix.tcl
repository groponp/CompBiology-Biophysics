##; script converto segname to chainID 
##; @Autor: Ropón-Palacios G. 
###; 14 Feb 2021. 

mol new step1_pdbreader.pdb type pdb

set all [atomselect top "all"]

set segnames [lsort -unique [$all get segname]]
foreach seg $segnames {
        set selseg [atomselect top "segname $seg"]
        $selseg set chain [string index $seg 3]
}

$all writepdb orient.pdb

quit


