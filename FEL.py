#! Calculate Free Energy Lanscape - FEL. 
#! _author : Ropón-Palacios G. 
#! _date   : March 9, 2022. 
#! _e-mail : groponp@gamil.com 

import optparse 
import subprocess 

import warnings
warnings.filterwarnings("ignore")

disclaimer="""<Free Energy Landscape Tool>"""
parser = optparse.OptionParser(description=disclaimer) 
#! INPUTS options
parser.add_option("--file1", help="collective variable 1 [rgyr.dat]", type=str)
parser.add_option("--file2", help="collective variable 2 [rmsd.dat]", type=str)
parser.add_option("--temperature", help="temperature in Kelvin [i.e 310.0 ]", type=float)
parser.add_option("--bin", help="bin  [default 25 ]", type=int)

#! OUTPUT options 
#parser.add_option("--x_label", help=" r\"RMSD [$\AA$]\" ", type=str, action='store') 
#parser.add_option("--y_label", help=" r\"Rgyr [$\AA$]\" ", type=str, action='store')
parser.add_option("--ofile", help="type output name [FEL_LIG.svg]", type=str)

options, args = parser.parse_args() 


def free_energy_surface(file1, file2, temperature, i1):
    import math
    import numpy as np

    i2 = i1
    T = float(temperature)

    V = np.zeros((i1,i2))
    DG = np.zeros((i1,i2))

    kB = 3.2976268E-24 #cal/K
    An = 6.02214179E23

    minv1 = np.min(file1)
    maxv1 = np.max(file1)
    minv2 = np.min(file2)
    maxv2 = np.max(file2)

################### Data span ####################
    I1 = maxv1 - minv1
    I2 = maxv2 - minv2

####################### Binning #####################
    for i in range(len(v1)):
         for x in range(i1):
            if v1[i] <= minv1+(x+1)*I1/i1 and v1[i] > minv1+x*I1/i1:
                 for y in range(i2):
                     if v2[i] <= minv2+(y+1)*I2/i2 and v2[i] > minv2+y*I2/i2:
                         V[x][y] = V[x][y] +1
                         break
                 break

##### Finding the maximum ##############
    P = list()
    for x in range(i1):
            for y in range(i2):
                    P.append(V[x][y])

    Pmax = max(P)

##### Calculating Delta G values ##############
    LnPmax = math.log(Pmax)

    for x in range(i1):
        for y in range(i2):
            if V[x][y] == 0:
                DG[x][y] = 10
                continue
            else:
                DG[x][y] = -0.001*An*kB*T*(math.log(V[x][y])-LnPmax) #kcal/mol

    return DG

#### Funcion para generar plots de FEL
def plot_fel(file1,file2,dataframe,labels,titulo):
    z_l = r'$\Delta G$'+' [kcal/mol]' #using latex in matplotlib
    rangos = [file2.min(),file2.max(),file1.min(),file1.max()]

    ##; axs 1

    fig, ax1 = plt.subplots()

    im1 = ax1.matshow(dataframe, cmap='jet',
                      extent=rangos,
                      origin='lower',
                      interpolation='bilinear',
                      aspect='auto')

    ax1.tick_params(axis='both', labelsize=9)

    ax1.set_xlabel(labels[0], fontsize=15, labelpad=1)
    ax1.set_ylabel(labels[1], fontsize=15,labelpad=-1)
    ax1.xaxis.set_ticks_position('bottom')

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.09)
    cbar = plt.colorbar(im1, cax=cax)
    cbar.set_label(z_l,size=14, labelpad = -3)
    fig.set_size_inches(3.5, 2.5, forward=True)
    #plt.tight_layout()
    
    plt.savefig(options.ofile, dpi=300, format="png", bbox_inches='tight') 
    #plt.savefig(titulo, dpi=600, format="svg")

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
sns.set()

#plt.style.use("seaborn") 
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['figure.titlesize'] = 15

#types = { 'roto': r"$RMSD/Å$",
 #         'angl': r"$Rgyr/Å$"}

#! Define all varibales I/O 
file1 = options.file1
file2 = options.file2
ofile = options.ofile 
#xlabel = options.x_label 
#ylabel = options.y_label 
temp = options.temperature 
bin_h = options.bin 

#! Run routines
v1 = np.genfromtxt(file1, usecols=1, delimiter='\t', skip_header=1)
v2 = np.genfromtxt(file2, usecols=1, delimiter='\t', skip_header=1)

fel = free_energy_surface(v1,
                          v2, temp, bin_h)

plot_fel(v1,v2,fel,[r"Rgyr [$\AA$]", r"RMSD [$\AA$]"], ofile)
