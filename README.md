# CompBiology-Biophysics
It's an general repository containing multiple script for several task into computational biophysics or biology, which should be update frequently for me.
If you've any questions, please write me to groponp@gmail.com. Programming lenguage used into this repo are: C++, Fortran95, Tcl, Python, Bash, AWK. 

If you're using any code, please cite our papers: 
1. Ropón-Palacios et al.  (2020) Potential novel inhibitors against emerging zoonotic pathogen Nipah virus: a virtual screening and molecular dynamics approach, Journal of Biomolecular Structure and Dynamics, 38:11, 3225-3234, DOI: 10.1080/07391102.2019.1655480
2. Ropón Palacios et al. (2019) Novel multi-epitope protein containing conserved epitopes from different Leishmania species as potential vaccine candidate: Integrated immunoinformatics and molecular dynamics approach, Computational Biology and Chemistry, DOI: https://doi.org/10.1016/j.compbiolchem.2019.107157.
3. Otazu et al. (2020) Targeting Receptor Binding Domain and Cryptic Pocket of Spike glycoprotein from SARS-CoV-2 by biomolecular modeling, ARXIV Quantitative Biology, Biomolecules, DOI: https://doi.org/10.48550/arXiv.2006.06452
4. Ropón-Palacios e tal. (2022) Glycosylation is key for enhancing drug recognition into spike glycoprotein of SARS-CoV-2, Computational Biology and Chemistry, DOI: https://doi.org/10.1016/j.compbiolchem.2022.107668
5. Osorio-Mogollón et al. (2022) Attacking the SARS-CoV-2 Replication Machinery with the Pathogen Box’s Molecules, Letters in Drug Desing & Discovery, DOI: 10.2174/1570180819666220622085659 

====================================================================================
```
It secction given a brief overview of all tools here: 
1. BigBabel.py : It is a script that allows, to convert a large base (> 1 million molecules) of data in 1D-SDF format for PBQT useful to run massive molecular docking.
2. 2_movie_render_fix.tcl: This script generates a video from a molecular dynamics trajectory in full HD, for this you have to have a VMD scene with all the parameters set (amb. occl, light, material, color, representation), that are of your interest.
3. Co_mol_md.py: This is a simple script, which calculates the number of molecules that can be placed in a simulation box, based on their molar concentration, i.e. Imagine that you want to place 10 mM of Urea in your simulation box, this script will be able to calculate how many molecules you need to place to reach that concentration.
4. DeltGtoKd.py/cpp: These two scripts do the same task, use the deltaG value of autodokc4/vina to calculate the dissociation constant (Kd) for a given molecule. In case of the python script it is useful for data < 1 Million , while the one written in C++ is more useful for data > 1 Million.
5. FEL.py: It is a physical algorithm, which allows the use of two coordinate reactions (commonly RMSD, Rgyr) for a given biological system in order to find its Free Energy lanscape, and determine the sampled metastases. It generate un 2D plot too. 

6. Getseq.py: Este es un script, que permite obtener todas las posibles secuencias homólogas para un dado query de secuencias, y para taxones específicos, ya sea para modelar o realizar alineamientos múltiples. 
7. MDEquiStats.py: Este es un script que permite extraer la información de la presión, temperatura, energía durante el proceso de equilibración en la dinámica molecular usando NAMD. 
8. MakeMSA.R: Este es un simple código en R que permite hacer alineamiento de calidad y perfonalizarlos, en formato LaTex. 
9. MolPack.tcl: Es largo script, que permite preparar sistemas den solución (como proteina en agua) y sistemás para realizar dinámica molecular de estiramiento. 
10. OrientZ-axis.py: Este escript tomado una molécula y la orienta en Z axis, para su posterior uso, como en dinámica molecular de estiramiento. 
11. autodock-vina-screening.sh: Es un script que permite correr virtual screening para una base de datos larga. 
12. contact_map_prot.py: Este script permite encontrar el número de contactos que se da entre dos moléculas que interaccionan. Principalmente entre proteína-proteína. 
13. convert-namd2charmm.py & convert-namd4gmx.py: Estos script convierten la topología y las coordenadas de NAMD/GRO para charmm. Tenga cuidao con esto, podría tener algunos errores si los archivos de input no son los correctos.
14. essential-dynamics.sh: Este script permite procesar trayectorias de dinámica molecular obtener los modos normales del movimiento de una dada biomolécula. 
15. get_box.tcl: Este script permite obtener las dimenciones y el del grid de docking, esto permite usar moléculas con más de 1M de átomos, sobrellevando la limitación de MGLTools. 
16. harm-potential-us.py: Es un script que permite estimar la constante k en el bias potential 1/2k(w0-w1), para realizar umbrella sampling. 
17. jazynski-tclforces.tcl: Este es un script basado en TCLForces 
```
