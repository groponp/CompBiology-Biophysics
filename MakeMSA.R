#! Perform MSA analysis for 5HT protein from Trichinella spiralis 
#! @utor: Ropón-Palacios G. 
#! date: 20, Jul, 2022. 
#! e-mail: groponp@gmail.com 

setwd("/Volumes/Kepler/Research/5HT_Tspiralis")

#! Define function Tool-kits.
#!===================================================
install_msa <- function() {
  #! Function to install It package if is necessary. 
  if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
  BiocManager::install("msa")
  
}

run_msa <- function(seqfile) {
  seqFile <- readAAStringSet(seqfile)
  x
  seqAln <- msa(seqFile, method="ClustalW", verbose=TRUE)
  print(seqAln)
  msaPrettyPrint(seqAln, output="tex", showNames="left", showLogo="top",
                 askForOverwrite=FALSE, shadingColors="black",
                  consensusColors="HotCold", logoColors="rasmol", 
                 shadingMode="similar")
}

#! use only if you haven't installed it package.
#======================================================
# install_msa() 

#! Importing library and setting parameters of library. 
#======================================================
library("msa")
system.file("tex", "texshade.sty", package="msa")
file.copy("/Library/Frameworks/R.framework/Versions/4.2/Resources/library/msa/tex/texshade.sty", 
          ".", overwrite = TRUE) 

#! Call function. 
#=================================
run_msa(seqfile="5HT_MSA_fix.fasta") 
tools::texi2pdf("seqAln.tex", clean=TRUE)


