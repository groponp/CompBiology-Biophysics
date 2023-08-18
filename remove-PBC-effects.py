#########-#########-#########-#########-#########-#########-#########-#########
#! Remove PBC effects during MD 
#! __author__ = 'Ropón-Palacios G.' 
#! __date__   = 'November 28, 2021'
#! check: pass
#! functional : pass  
#########-#########-#########-#########-#########-#########-#########-#########

#########-#########-#########-#########-#########-#########-#########-######### 
# modified Mon 15 Aug 2023.  
# Change logs:
# - It now version is util for NAMD and GROMACS into membrane or solutions.
# - Not use in_memory speed transformations, but write is slow.
# - Add secction into pbFix method to enable pbctools from VMD (more fast)
#########-#########-#########-#########-#########-#########-#########-#########

import MDAnalysis as mda 
from MDAnalysis import transformations  as trans 
import warnings 
warnings.filterwarnings('ignore')
import optparse
import os

#########-#########-#########-#########-#########-#########-#########-#########
# Class definitions
#########-#########-#########-#########-#########-#########-#########-#########
class MOLTransformation:

	def __init__(self, top:str, traj:str, otraj:str, sel1:str, sel2:str, typeFix:str, method:str,
	 stride=int(1)): 
		
		self.stride = stride                        # Stride to pass frame 
		self.top = top                              # Topology 
		self.traj = traj                            # Trajectory
		self.otraj = otraj 
		if method != "vmd":
			self.u = mda.Universe(self.top, self.traj)#,
		#in_memory=True, in_memory_step=1)  # Universe 
		#self.otraj = otraj                          # Name of oFile i.e 
													# `oFile.xtc/dcd`
		self.selPBC = sel1                          # Select Atoms to fix PBC
		self.selOFile = sel2                        # Select Atoms to save Traj
		self.type = typeFix  

	def pbcFix(self, method:str): 
		"""
		It method use VMD or MDAnalysis transformation for remove PBC from membrane 
		or solution. 
		Notes: 
		- If you system have one protein, then select only `protein`. 
		- If you system have two proteins, then select the protein chain 
			`segname PROA`. 
		
		Examples:
			>> molt = MOLTransformation(top, traj, sel1, sel2, typeFix, 
						stride) 
			>> molt.pbcFix(method=vmd) 
		""" 

		# It's VMD method , more fast that MDA
		#########-#########-#########-#########-#########-#########-#########-#########
		if method == "vmd":
			ext_top  = self.top.split(".")[1]
			ext_traj = self.traj.split(".")[1]
			outname  = self.otraj.split(".")[0] + "_noPBC." + self.otraj.split(".")[1]
		
			f = open("pbc.tcl", "w+")
			f.write("package require pbctools\n")
			f.write("mol new {} type {}\n".format(self.top, ext_top))
			f.write("mol addfile {} type {} step {} waitfor all\n".format(self.traj, ext_traj, self.stride))
				
			if self.type == "solution":
				print("[INFO   ] Starting PBC fix from solution.")
				f.write("pbc wrap -center com -centersel {} -all -compound res\n".format(self.selPBC))
				f.write("set sel [atomselect top \"{}\"]\n".format(self.selOFile))
				f.write("animate write {} {} beg 0 end -1 waitfor all sel $sel\n".format(ext_traj, outname, 
					self.selOFile))
				f.write("quit")
				f.close()
				os.system("vmd -dispdev text -e pbc.tcl") 

			elif self.type == "membrane":
				print("[INFO   ] Starting PBC fix from membrane.")
				f.write("pbc wrap -center com -centersel {} -all -compound res\n".format(self.selPBC))
				f.write("set sel [atomselect top \"{}\"]\n".format(self.selOFile))
				f.write("animate write {} {} beg 0 end -1 waitfor all sel $sel\n".format(ext_traj, outname, 
					self.selOFile))
				f.write("quit")
				f.close()
				os.system("vmd -dispdev text -e pbc.tcl") 

		# It's MDA method , more slow that VMD
		#########-#########-#########-#########-#########-#########-#########-#########
		else:
			if self.type == "solution":
				print("[INFO   ] Starting PBC fix from solution.")
				protein = self.u.select_atoms(selPBC) 
				waters = self.u.select_atoms("not {}".format(selPBC)) 
				workflow = [trans.unwrap(protein, max_threads=4),
			 	trans.center_in_box(protein, center="geometry", max_threads=4, 
			 	), trans.wrap(waters, compound="residues", max_threads=4, 
			 	)]
				self.u.trajectory.add_transformations(*workflow)
				print("[INFO   ] Finish PBC fix from solution.")

			elif self.type == "membrane":
				print("[INFO   ] Starting PBC fix from membrane.")
				protein = self.u.select_atoms(selPBC) 
				ag = self.u.atoms 
				workflow = (trans.unwrap(ag), trans.center_in_box(protein, center="mass"),
					trans.wrap(ag, compound="fragment"))
				self.u.trajectory.add_transformations(*workflow)
				print("[INFO   ] Finish PBC fix from solution.")

	def writeTRAJ(self, otraj, u):
		"""
		It method use MDAnalysis to write traj output (so slow now). 
		Notes: 
		- Atoms into TRAJ are only save from `self.OFile`.
		Examples:
			>> molt = MOLTransformation(top, traj, sel1, sel2, typeFix, 
						stride) 
			>> molt.writeTRAJ(otraj=name_fit, u=universe) 
		"""
		
		selAtoms = u.select_atoms(self.selOFile) 

		# Write output traj
		#########-#########-#########-#########-#########-#########-#########-#########
		with  mda.Writer(otraj, selAtoms.n_atoms) as W:
			counter = 1 
			total   = len(u.trajectory)
			#u.transfer_to_memory(verbose=True)
			for ts in u.trajectory:
				frac = (counter/total)
				per  = frac * 100 
				W.write(selAtoms)
				print("Write Frame [{}/{} - {:.2f}%]".format(counter,total, per))
				counter +=1 

	def fit(self, noPBC_traj):
		"""
		It method use MDAnalysis to remove rot+trans in the TRAJ. 
		Notes: 
		- It use an TRAJ upload in memory to speed the write process.
		Examples:
			>> trajnoPBC = mda.Universe(itop, name_noPBC, n_memory=True, in_memory_step=1) 
			>> u1 = molt.fit(noPBC_traj=trajnoPBC) 
			>> molt.writeTRAJ(otraj=name_fit, u=u1)
		"""

		if self.type == "solution":
			u1 = noPBC_traj
			ref_u1 = u1.copy() 
			reference = ref_u1.select_atoms(self.selPBC)
			prot = u1.select_atoms(self.selPBC)

			workflow = trans.fit_rot_trans(prot, reference)
			print("[INFO   ] Fitting rot+trans complete to biomolecule into solution.")
			u1.trajectory.add_transformations(workflow) 
			return u1 
			
		elif options.type == "membrane":
			u1 = noPBC_traj
			ref_u1 = u1.copy()
			reference = ref_u1.select_atoms(self.selPBC)
			prot = u1.select_atoms(self.selPBC)
			workflow = trans.fit_rot_trans(prot, reference, plane='xy'
				, weights="mass") 
			print("[INFO   ] Fitting membrane rotxy+transxy complete.")
			u1.trajectory.add_transformations(workflow) 
			return u1 
			


#########-#########-#########-#########-#########-#########-#########-#########
# CLI 
#########-#########-#########-#########-#########-#########-#########-#########

info = '''Remove and Fit PBC Effect into GMX and NAMD TRAJ'''
parser = optparse.OptionParser(description=info) 
parser.add_option("--top", help="Topology [TPR or PSF] type=str)", type=str) 
parser.add_option("--itraj", help="Input Trajectory [XTC, DCD]", type=str)
parser.add_option("--selPBC", help="Selection PBC atoms to centred, MDAnalysis/VMD-based i.e \
	[\"protein\" or \"segname PROA\"]", type=str, action='store')
parser.add_option("--selOFile", help="Selection atoms to output, MDAnalysis/VMD-based i.e\
	[\"protein\" or \"all\"]", type=str, action='store')


parser.add_option("--method", help="Select method [vmd or mda]", type=str)
parser.add_option("--steps", help="Slicing frame, useful only in VMD", type=int)
parser.add_option("--type", help="Type system [solution or membrane]", type=str)
parser.add_option("--otraj", help="output name of TRAJ [md_0_300.xtc]", type=str)
parser.add_option("--fittraj", help="Enable fit rot+trans", type=str)
parser.add_option("--usage", help="Print an example the as use it", action="store_true")

					
options, args = parser.parse_args() 

#########-#########-#########-#########-#########-#########-#########-#########
# IO
#########-#########-#########-#########-#########-#########-#########-#########

itop     = options.top
itraj    = options.itraj
otraj    = options.otraj
selPBC   = options.selPBC
selOFile = options.selOFile
typeFix  = options.type 
method    = options.method
fit      = options.fittraj
usage    = options.usage 
steps    = options.steps

if usage:
		print("Example: python remove-PBC-effects.py --top name.psf --itraj namd.dcd" +
			" --selPBC \"protein\" --selOFile \"all\" --otraj md_0_300.dcd" +
			" --fittraj true --type solution --method vmd --step 60")

else: 
	name_noPBC = otraj.split(".")[0] + "_noPBC." + otraj.split(".")[1]
	name_fit   = otraj.split(".")[0] + "_noPBC_fit." + otraj.split(".")[1]


	molt = MOLTransformation(itop, itraj, otraj, selPBC, selOFile, typeFix, method, steps)
	molt.pbcFix(method = method)
	if method == "mda":
		molt.writeTRAJ(otraj=name_noPBC, u = molt.u)

		if fit == "true":
			trajnoPBC = mda.Universe(itop, name_noPBC, in_memory=True, in_memory_step=1)  
			u1 = molt.fit(noPBC_traj=trajnoPBC) 
			molt.writeTRAJ(otraj=name_fit, u=u1)

	else: 
		if fit == "true":
			trajnoPBC = mda.Universe(itop, name_noPBC)#, in_memory=True, in_memory_step=1)  
			trajnoPBC.transfer_to_memory(verbose=True)
			u1 = molt.fit(noPBC_traj=trajnoPBC) 
			molt.writeTRAJ(otraj=name_fit, u=u1)









