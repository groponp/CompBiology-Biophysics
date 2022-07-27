###; Make segments name for easy analysis. 
###; By: Ropón-Palacios G. <groponp@gmail.com>.


mol new step3_input.pdb type pdb waitfor all 


##; Set Spike region for the ectodomain 
##; example to select into analysis: 
##; segid RDB RMB and name CA ---> to select all RDB 

set rbd [atomselect top "protein and resid 509 to 541"]  ;## for RBD region, RDB total : RDB + RBM 
set rbm [atomselect top "protein and resid 437 to 508"]  ;## for RBM region 
set ntd [atomselect top "protein and resid 13 to 303"]   ;## for NTD region, NTD total : NTD + S1 
set s1 [atomselect top "protein and resid 304 to 685"]   ;## for S1 region 
set s2 [atomselect top "protein and resid 686 to 1273"]  ;## for S2 region 

$rbd set segid "RBD"
$rbm set segid "RBM"
$ntd set segid "NTD"
$s1  set segid "S1" 
$s2  set segid "S2"


###; set glycosilation, ligand and other.  
set md1 1
set md2 0
set md3 0
set md4 0
set md5 0
set md6 0

if { $md1 } {
	set lig [atomselect top "resname UNL"]
	$lig set resname "LIG"
	set pocket [atomselect top "protein and same residue as within 4 of resname LIG"]
	set presid [$pocket get resid] 
	set ofile [open "pocket_resids.dat" w]
	puts $ofile "$presid" 
	$pocket set segid "BPO" ;## for binding pocket into RBD place. 
	
	set all [atomselect top "all"]
	$all writepdb md1_renamed.pdb 
	quit 
	
}  


if { $md2 } {
	set lig [atomselect top "resname UNL"]
	$lig set resname "LIG"
	set pocket [atomselect top "protein and same residue as within 4 of resname LIG"]
	$pocket set segid "BPO" ;## for binding pocket into RBD place. 
	
	##; select all residue that are glycosilated
	set r1 [atomselect top "protein and resid 17"]
	set g1 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL BGLC AMAN BGLC AFUC"]
	$r1 set segid "R1"
	$g1 set segid "G1"

	set r2 [atomselec top "protein and resid 616"]
	set g2 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"] 
	$r2 set segid "R2"
	$g1 set segid "G2"
	
	set r3 [atomselecct top "protein and resid 61"]
	set g3 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r3 set segid "R3"
	$g3 set segid "G3"

	set r4 [atomselect top "protein and resid 74"]
	set g4 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL ANE5 BGLC AMAN BGLC BGLC AFUC"]
	$r4 set segid "R4"
	$g4 set segid "G4"

	set r5 [atoselect top "protein and resid 122"]
	set g5 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r5 set segid "R5"
	$g5 set segid "G5"

	set r6 [atomselect top "protein and resid 149"]
	set g6 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r6 set segid "R6"
	$g6 set sedig "G6"

	set r7 [atomselect top "protein and resid 165"]
	set g7 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL"]
	$r7 set segid "R7"
	$g7 set segid "G7"

	set r8 [atomselect top "protein and resid 331"]
	set g8 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r8 set segid "R8"
	$g8 set segid "G8"

	set r9 [atomselect top "protein and resid 234"]
	set g9 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN AMAN AMAN AMAN"]
	$r9 set segid "R9"
	$r9 set segid "G9"

	set r10 [atomselect top "protein and resid 657"]
	set g10 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN BGLC BGAL"]
	$r10 set segid "R10"
	$g10 set segid "G10"

	set r11 [atomselect top "protein and resid 282"]
	set g11 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGLC AMAN BGLC AFUC"] 
	$r11 set segid "R11"
	$g11 set segid "G11"

	set r12 [atomselect top "protein and resid 323"]
	set g12 [atomselect top "resname AGAL BGAL ANE5"]
	$r12 set segid "R12"
	$g12 set segid "G12"

	set r13 [atomselect top "protein and resid 343"]
	set g13 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r13 set segid "R13"
	$g13 set segid "G13"

	set r14 [atomselect top "protein and resid 603"]
	set g14 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r14 set segid "R14"
	$g14 set segid "G14"

	set r15 [atomselect top "protein and resid 709"]
	set g15 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r15 set segid "R15"
	$g15 set segid "R15"

	set r16 [atomselect top "protein and resid 717"]
	set g16 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r16 set segid "R16"
	$g16 set segid "G16"

	set r17 [atomselect top "protein and resid 801"]
	set g17 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r17 set segid "R17"
	$g17 set segid "G17"

	set r18 [atomselect top "protein and resid 1074"]
	set g18 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r18 set segid "R18"
	$g18 set segid "G18"

	set r19 [atomselect top "protein and resid 1098"]
	set g19 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL ANE5 AMAN BGLC AFUC"]
	$r19 set segid "R19"
	$g19 set segid "G19"

	set r20 [atomselect top "protein and resid 1134"]
	set g20 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"] 
	$r20 set segid "R20"
	$g20 set segid "G20"
	
	set all [atomselect top "all"]
	$all writepdb md2_renamed.pdb 
	quit 
	
}  
 

if { $md3 } {
	set lig [atomselect top "resname UNL"]
	$lig set resname "LIG"
	set cryptic [atomselect top "protein and same residue as within 4 of resname LIG"]
	set cresid [$cryptic get resid] 
	set ofile [open "cryptic_resids.dat" w]
	puts $ofile "$cresid" 
	$cryptic set segid "CPO" ;## for cryptic pocket pocket into RBD place. 
	
	set all [atomselect top "all"]
	$all writepdb md3_renamed.pdb 
	quit 
	
}  


if { $md4 } {
	set lig [atomselect top "resname UNL"]
	$lig set resname "LIG"
	set cryptic [atomselect top "protein and same residue as within 4 of resname LIG"]
	$cryptic set segid "CPO" ;## for cryptic pocket into RBD place. 
	
	##; select all residue that are glycosilated
	set r1 [atomselect top "protein and resid 17"]
	set g1 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL BGLC AMAN BGLC AFUC"]
	$r1 set segid "R1"
	$g1 set segid "G1"

	set r2 [atomselec top "protein and resid 616"]
	set g2 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"] 
	$r2 set segid "R2"
	$g1 set segid "G2"
	
	set r3 [atomselecct top "protein and resid 61"]
	set g3 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r3 set segid "R3"
	$g3 set segid "G3"

	set r4 [atomselect top "protein and resid 74"]
	set g4 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL ANE5 BGLC AMAN BGLC BGLC AFUC"]
	$r4 set segid "R4"
	$g4 set segid "G4"

	set r5 [atoselect top "protein and resid 122"]
	set g5 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r5 set segid "R5"
	$g5 set segid "G5"

	set r6 [atomselect top "protein and resid 149"]
	set g6 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r6 set segid "R6"
	$g6 set sedig "G6"

	set r7 [atomselect top "protein and resid 165"]
	set g7 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL"]
	$r7 set segid "R7"
	$g7 set segid "G7"

	set r8 [atomselect top "protein and resid 331"]
	set g8 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r8 set segid "R8"
	$g8 set segid "G8"

	set r9 [atomselect top "protein and resid 234"]
	set g9 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN AMAN AMAN AMAN"]
	$r9 set segid "R9"
	$r9 set segid "G9"

	set r10 [atomselect top "protein and resid 657"]
	set g10 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN BGLC BGAL"]
	$r10 set segid "R10"
	$g10 set segid "G10"

	set r11 [atomselect top "protein and resid 282"]
	set g11 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGLC AMAN BGLC AFUC"] 
	$r11 set segid "R11"
	$g11 set segid "G11"

	set r12 [atomselect top "protein and resid 323"]
	set g12 [atomselect top "resname AGAL BGAL ANE5"]
	$r12 set segid "R12"
	$g12 set segid "G12"

	set r13 [atomselect top "protein and resid 343"]
	set g13 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r13 set segid "R13"
	$g13 set segid "G13"

	set r14 [atomselect top "protein and resid 603"]
	set g14 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r14 set segid "R14"
	$g14 set segid "G14"

	set r15 [atomselect top "protein and resid 709"]
	set g15 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r15 set segid "R15"
	$g15 set segid "R15"

	set r16 [atomselect top "protein and resid 717"]
	set g16 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r16 set segid "R16"
	$g16 set segid "G16"

	set r17 [atomselect top "protein and resid 801"]
	set g17 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r17 set segid "R17"
	$g17 set segid "G17"

	set r18 [atomselect top "protein and resid 1074"]
	set g18 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r18 set segid "R18"
	$g18 set segid "G18"

	set r19 [atomselect top "protein and resid 1098"]
	set g19 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL ANE5 AMAN BGLC AFUC"]
	$r19 set segid "R19"
	$g19 set segid "G19"

	set r20 [atomselect top "protein and resid 1134"]
	set g20 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"] 
	$r20 set segid "R20"
	$g20 set segid "G20"
	
	set all [atomselect top "all"]
	$all writepdb md4_renamed.pdb 
	quit 
	
}  


if { $md5 } {
	##; revisar los archivos generados en MD1 y MD3 para obtener los residuos 
	set cryptic [atomselect top "protein and resid"]
	set pocket  [atomselect top "protein and resid "]

	$pocket set segid "BPO" ;## for binding pocket pocket into RBD place. 
	$cryptic set segid "CPO" ;## for cryptic pocket pocket into RBD place. 
	
	set all [atomselect top "all"]
	$all writepdb md5_renamed.pdb 
	quit 
	
}  


if { $md6 } {
	##; revisar los archivos generados en MD1 y MD3 para obtener los residuos 
	set cryptic [atomselect top "protein and resid"] ;## adicionar los residues. 
	set pocket  [atomselect top "protein and resid "]

	$pocket set segid "BPO" ;## for binding pocket pocket into RBD place. 
	$cryptic set segid "CPO" ;## for cryptic pocket pocket into RBD place.

		##; select all residue that are glycosilated
	set r1 [atomselect top "protein and resid 17"]
	set g1 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL BGLC AMAN BGLC AFUC"]
	$r1 set segid "R1"
	$g1 set segid "G1"

	set r2 [atomselec top "protein and resid 616"]
	set g2 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"] 
	$r2 set segid "R2"
	$g1 set segid "G2"
	
	set r3 [atomselecct top "protein and resid 61"]
	set g3 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r3 set segid "R3"
	$g3 set segid "G3"

	set r4 [atomselect top "protein and resid 74"]
	set g4 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL ANE5 BGLC AMAN BGLC BGLC AFUC"]
	$r4 set segid "R4"
	$g4 set segid "G4"

	set r5 [atoselect top "protein and resid 122"]
	set g5 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r5 set segid "R5"
	$g5 set segid "G5"

	set r6 [atomselect top "protein and resid 149"]
	set g6 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r6 set segid "R6"
	$g6 set sedig "G6"

	set r7 [atomselect top "protein and resid 165"]
	set g7 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL"]
	$r7 set segid "R7"
	$g7 set segid "G7"

	set r8 [atomselect top "protein and resid 331"]
	set g8 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r8 set segid "R8"
	$g8 set segid "G8"

	set r9 [atomselect top "protein and resid 234"]
	set g9 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN AMAN AMAN AMAN"]
	$r9 set segid "R9"
	$r9 set segid "G9"

	set r10 [atomselect top "protein and resid 657"]
	set g10 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN BGLC BGAL"]
	$r10 set segid "R10"
	$g10 set segid "G10"

	set r11 [atomselect top "protein and resid 282"]
	set g11 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGLC AMAN BGLC AFUC"] 
	$r11 set segid "R11"
	$g11 set segid "G11"

	set r12 [atomselect top "protein and resid 323"]
	set g12 [atomselect top "resname AGAL BGAL ANE5"]
	$r12 set segid "R12"
	$g12 set segid "G12"

	set r13 [atomselect top "protein and resid 343"]
	set g13 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]
	$r13 set segid "R13"
	$g13 set segid "G13"

	set r14 [atomselect top "protein and resid 603"]
	set g14 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r14 set segid "R14"
	$g14 set segid "G14"

	set r15 [atomselect top "protein and resid 709"]
	set g15 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r15 set segid "R15"
	$g15 set segid "R15"

	set r16 [atomselect top "protein and resid 717"]
	set g16 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r16 set segid "R16"
	$g16 set segid "G16"

	set r17 [atomselect top "protein and resid 801"]
	set g17 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r17 set segid "R17"
	$g17 set segid "G17"

	set r18 [atomselect top "protein and resid 1074"]
	set g18 [atomselect top "resname BGLC BGLC BMAN AMAN AMAN AMAN AMAN"]
	$r18 set segid "R18"
	$g18 set segid "G18"

	set r19 [atomselect top "protein and resid 1098"]
	set g19 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC BGAL ANE5 AMAN BGLC AFUC"]
	$r19 set segid "R19"
	$g19 set segid "G19"

	set r20 [atomselect top "protein and resid 1134"]
	set g20 [atomselect top "resname BGLC BGLC BMAN AMAN BGLC AMAN BGLC AFUC"]  
	$r20 set segid "R20"
	$g20 set segid "G20"	
	
	set all [atomselect top "all"]
	$all writepdb md5_renamed.pdb 
	quit 
	
}  

puts "thanks for use."
