##; Calculate contact with more of 70% along trajectory.
##; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##; @autor: Ropón-Palacios G. 
##; e-mail: groponp@gmail.com or ggdrpalacios@uesc.br 
##; date: September 14, 2021. 

import MDAnalysis as mda  
import numpy as np 
from MDAnalysis.analysis import distances
import pandas as pd

##; load data in this case for split traj 
##; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
top = "wt_150mM.hmr.psf"
traj_list = ["traj1.dcd", "traj2.dcd", "traj3.dcd"]                ##; list of traj to contact or read into universe
u = mda.Universe(top, [traj_list[0], traj_list[1], traj_list[2]])  ##; this approach for concatenate multiple traj. 
nframes = len(u.trajectory)

##; Select group atoms para measure distance and contact
##; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
proa = u.select_atoms("segid PROA and name CA")
prob = u.select_atoms("segid PROB and name CA")
proc = u.select_atoms("segid PROC and name CA")  

##; Core of script, this calculate contact 
##; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_proa = len(proa)  ##; for RBD
n_prob = len(prob)  ##; Chain B Antibody 
n_proc = len(proc)  ##; Chain C Antibody 

#++ for A --> B contact 
contact_sum_ab = np.zeros((n_proa, n_prob))

#++ for A ---> C contact 
contact_sum_ac = np.zeros((n_proa, n_proc))

#++ cut-off to determine contacts 
cutoff = 7

#++ Looping througth all frames and canculte matrix mean 
for ts in u.trajectory[0:50]:  ##; for A --> B
    p1 = proa.positions 
    p2 = prob.positions
    ts_dist_ab = distances.distance_array(p1,p2, box=u.dimensions)
    ts_dist_ab[ts_dist_ab < cutoff] = 1  ###; for contact
    ts_dist_ab[ts_dist_ab > cutoff] = 0
    contact_sum_ab +=  ts_dist_ab

for ts in u.trajectory[0:50]: ##; for A ---> C 
    p1 = proa.positions 
    p3 = proc.positions
    ts_dist_ac = distances.distance_array(p1,p2, box=u.dimensions)
    ts_dist_ac[ts_dist_ac < cutoff] = 1  ###; for contact
    ts_dist_ac[ts_dist_ac > cutoff] = 0
    contact_sum_ac +=  ts_dist_ac

    
contact_ratio_ab = contact_sum_ab/nframes 
contact_ratio_ac = contact_sum_ac/nframes 


##; Filter data and write excel file. 
##; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df0 = pd.DataFrame(contact_ratio_ab)
df1 = pd.DataFrame(contact_ratio_ac)

#cols = prob.resids 
#rows = proa.resids
#df.columns = cols 
#df.index = rows

df3 = df0.where(df > 0.7)
df4 = df1.where(df > 0.7)


df3.to_csv("contacts_filter_07.csv")
df4.to_csv("contacts_filter_07.csv")


