# Computational Biology and Biophysics Repository
<div style="text-align: justify">
It's an general repository containing multiple script for several task into computational biophysics or biology, which should be update frequently for me.
If you've any questions, please write me to groponp@gmail.com. Programming lenguage used into this repo are: C++, Fortran95, Tcl, Python, Bash, AWK. 

If you want get more detail about this script can join to our discord channel: https://discord.gg/ZGpAWtvW 
</div>

## Programming language 
<img src="figures/python.png" width="40" height="40"> <img src="figures/gnu-bash.png" width="40" height="40"> <img src="figures/tcl.png" width="40" height="40"> <img src="figures/r.png" width="40" height="40"> <img src="figures/c-.png" width="40" height="40">

## Clone repo
```bash
git clone git@github.com:groponp/CompBiology-Biophysics.git
```

## Citations
If you're using any code, please cite our papers: 
1. Ropón-Palacios et al.  (2020) Potential novel inhibitors against emerging zoonotic pathogen Nipah virus: a virtual screening and molecular dynamics approach, Journal of Biomolecular Structure and Dynamics, 38:11, 3225-3234, DOI: 10.1080/07391102.2019.1655480
2. Ropón Palacios et al. (2019) Novel multi-epitope protein containing conserved epitopes from different Leishmania species as potential vaccine candidate: Integrated immunoinformatics and molecular dynamics approach, Computational Biology and Chemistry, DOI: https://doi.org/10.1016/j.compbiolchem.2019.107157.
3. Otazu et al. (2020) Targeting Receptor Binding Domain and Cryptic Pocket of Spike glycoprotein from SARS-CoV-2 by biomolecular modeling, ARXIV Quantitative Biology, Biomolecules, DOI: https://doi.org/10.48550/arXiv.2006.06452
4. Ropón-Palacios e tal. (2022) Glycosylation is key for enhancing drug recognition into spike glycoprotein of SARS-CoV-2, Computational Biology and Chemistry, DOI: https://doi.org/10.1016/j.compbiolchem.2022.107668
5. Osorio-Mogollón et al. (2022) Attacking the SARS-CoV-2 Replication Machinery with the Pathogen Box’s Molecules, Letters in Drug Desing & Discovery, DOI: 10.2174/1570180819666220622085659 
6. Atanda et al. (2023). In silico study revealed the inhibitory activity of selected phytomolecules of C. rotundus against VacA implicated in gastric ulcer, Journal of Biomolecular Structure and Dynamics, DOI: https://doi.org/10.1080/07391102.2022.2160814 


## Information of the Tools
It secction given a brief overview of all tools here: 
1. <img src="figures/python.png" width="15" height="15"> ***BigBabel.py:*** It is a script that allows, to convert a large base (> 1 million molecules) of data in 1D-SDF format for PBQT useful to run massive molecular docking.
2. <img src="figures/tcl.png" width="15" height="15"> ***2_movie_render_fix.tcl:*** This script generates a video from a molecular dynamics trajectory in full HD, for this you have to have a VMD scene with all the parameters set (amb. occl, light, material, color, representation), that are of your interest.
3. <img src="figures/python.png" width="15" height="15"> ***Co_mol_md.py:*** This is a simple script, which calculates the number of molecules that can be placed in a simulation box, based on their molar concentration, i.e. Imagine that you want to place 10 mM of Urea in your simulation box, this script will be able to calculate how many molecules you need to place to reach that concentration.
4. <img src="figures/python.png" width="15" height="15"> <img src="figures/c-.png" width="15" height="15"> ***DeltGtoKd.py/cpp:*** These two scripts do the same task, use the deltaG value of autodokc4/vina to calculate the dissociation constant (Kd) for a given molecule. In case of the python script it is useful for data < 1 Million , while the one written in C++ is more useful for data > 1 Million.
5. <img src="figures/python.png" width="15" height="15"> ***FEL.py:*** It is a physical algorithm, which allows the use of two coordinate reactions (commonly RMSD, Rgyr) for a given biological system in order to find its Free Energy lanscape, and determine the sampled metastases. It generate un 2D plot too. 
6. <img src="figures/python.png" width="15" height="15"> ***Getseq.py:*** This is a script, which allows to obtain all the possible homologous sequences for a given sequence query, and for specific taxa, either to model or perform multiple alignments.
7. <img src="figures/python.png" width="15" height="15"> ***MDEquiStats.py:*** This is a script that allows to extract the information of pressure, temperature, energy during the equilibration process in molecular dynamics using NAMD.
8. <img src="figures/r.png" width="15" height="15"> ***MakeMSA.R:*** This is a simple code in R that allows quality alignment and customization, in LaTex format.
9. <img src="figures/tcl.png" width="15" height="15"> ***MolPack.tcl:*** It is a long script, which allows to prepare systems in solution (such as protein in water) and systems to perform stretching molecular dynamics.
10. <img src="figures/python.png" width="15" height="15"> ***OrientZ-axis.py:*** This script takes a molecule and orients it on the Z axis, for later use, such as in molecular dynamics stretching.
11. <img src="figures/gnu-bash.png" width="15" height="15"> ***autodock-vina-screening.sh:*** It is a script that allows you to run virtual screening for a large database, from the conversion of the files to the processing of the results.
12. <img src="figures/python.png" width="15" height="15"> ***contact_map_prot.py:*** This script allows you to find the number of contacts between two interacting molecules. Mainly between protein-protein.
13. <img src="figures/python.png" width="15" height="15"> ***convert-namd2charmm.py & convert-namd4gmx.py:*** These scripts convert the topology and coordinates from NAMD/GRO to charmm. Be careful with this, you might get some errors if the input files are not correct.
14. <img src="figures/gnu-bash.png" width="15" height="15"> ***essential-dynamics.sh:*** This script allows processing molecular dynamics trajectories to obtain the normal modes of motion of a given biomolecule.
15. <img src="figures/tcl.png" width="15" height="15"> ***get_box.tcl:*** This script allows to obtain the dimensions and the docking grid, this allows to use molecules with more than 1M atoms, overcoming the limitation of MGLTools.
16. <img src="figures/python.png" width="15" height="15"> ***harm-potential-us.py:*** It is a script that allows estimating the constant k in the bias potential 1/2k(w0-w1), to carry out umbrella sampling.
17. <img src="figures/tcl.png" width="15" height="15"> ***jazynski-tclforces.tcl:*** This is a script based on NAMD's TCLForces, to pull non-equilibrium dynamics.
18. <img src="figures/tcl.png" width="15" height="15"> ***make_tclforces.tcl:*** This script generates the necessary files to be used in TCLForces.
19. <img src="figures/python.png" width="15" height="15"> ***make_flooding.py:*** This script allows to generate a flooding system, i.e. a concentration of N ligand in the extracellular part of an aquaporin inserted in a membrane.
20. <img src="figures/tcl.png" width="15" height="15"> ***make_segname.tcl:*** This is a script that allows you to add segname to the PDB files of a given system, to identify parts of a system, to facilitate analysis.
21. <img src="figures/python.png" width="15" height="15"> ****mda_2Dmatrix_fix.py:*** This generates a 2D matrix of the conformational arrangement of a biomolecule.
22. <img src="figures/python.png" width="15" height="15">  ***mda_convert-traj_fix.py:*** It interconverts NAMD/GROMACS path files.
23. <img src="figures/python.png" width="15" height="15"> ***mda_(rgyr,rmsd,rmsf)_fix.py:*** These scripts perform Rgyr, RMSD, RMSF analysis on a NAMD/GROMACS trajectory.
24. <img src="figures/tcl.png" width="15" height="15"> ***molywood_movie.tcl:*** This is a script based on the molywood library, which allows you to customize videos from simulation data in a more intuitive way.
25. <img src="figures/python.png" width="15" height="15"> ***prepare_charmm_gui-inputs.py:*** This generates input fix for use in charmm-gui, mainly for preparing protein-ligand systems. 
26. <img src="figures/python.png" width="15" height="15"> ***remove-PBC-effects.py:***  This script removes the effects of PBC in a dynamic and makes further analysis suitable.
27. <img src="figures/python.png" width="15" height="15"> ***remove_rot+trans.py:*** This script removes the rotation and translation movements to smooth the molecular dynamics trajectory and be able to visualize it without jump effects.
28. <img src="figures/gnu-bash.png" width="15" height="15"> ***running_analysis.sh:*** This is a meta-script that gives an example of running different scripts in python or another language in bash, to automate analysis on big data.
29. <img src="figures/gnu-bash.png" width="15" height="15"> ***sf2pdb_script.sh:*** This script uses openbabel to convert files from SDF to PDB.
30. <img src="figures/gnu-bash.png" width="15" height="15"> ***vitscreen_script.sh:*** This script uses autodockvina to run a virtual screening.
31. <img src="figures/tcl.png" width="15" height="15"> ***vmd_merge_pdb_fix.tcl:*** This script merges multiple pdb files into a single file using the TopoTools library.
32. <img src="figures/tcl.png" width="15" height="15"> ***vmd_reducetraj_fix.tcl:*** This script allows reducing the size of the trajectories of the molecular dynamics simulation, to have a better handling of them, on a workstation with little computational power.
33. <img src="figures/tcl.png" width="15" height="15"> ***vmd_rmsd_to_beta_fix.tcl:*** This script adds the average RMSD values for the beta column of the PDB file and allows coloring based on this value.
34. <img src="figures/python.png" width="15" height="15"> ***vmd_sasa_fix.py:*** This script allows to calculate the SASA value for a given selection.
35. <img src="figures/tcl.png" width="15" height="15"> ***vmd_segid_to_chain_fix.tcl:*** This script adds the string name for a PBD from the SEGID.
36. <img src="figures/gnu-bash.png" width="15" height="15"> ***RunHPCMD.sh:*** This is a long script that allows you to run NAMD/GROMACS molecular dynamics in a segmented manner. Example: Imagine that you have to run 100 ns of simulation, but you are running it on an HPC, which sets limited usage times per job (5 days maximum), and you cannot estimate the performance, so send the 100 ns in 5 days, it could result in an error, since it could not be completed, this script divides the simulation into discrete times and automatically resends them, until the entire simulation is completed.
37. <img src="figures/python.png" width="15" height="15"> ***make_gmx_atom_index.py:*** It is a simple script, which allows to create atom index to gromacs, with a very simple selection, which is very similar to VMD, but based on MDAnalysis. For example in gromacs it is very difficult to create an index to restrict atoms to within 5 A of a given atom, but here this can be done in a simple way, "name CA around 5 of rename LIG", for example.


## License 
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Disclaimer
Icons in this repo were taken from [flaticon](https://www.flaticon.com/free-icons/programming-language) 
