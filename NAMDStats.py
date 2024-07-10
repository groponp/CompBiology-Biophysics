#! * Extract Thermodynamics and Kinetics observables from Equilibration.
#! * @authors : Vega-Chozo K. & Ropón-Palacios G. 
#! * date : April 2, 2022. 
#! * e-mail : groponp@gmail.com 

#! * changeLogs:  `all changes was add by Ropón-Palacios G.`
#! =========================================================
#! a) Add print, to view into screen data average. Tue 27 Dec 20:00, 2022. 
#! b) Add conditional to create plot. Tue 27 Dec 20:11, 2022. 
#! c) Add routines to save data and plotting. Tue 27 Dec 20:40, 2022.
#! d) Change initial name MDEquiStats.oy to NAMD2EqStats.py. Tue 27 Dec 20:40, 2022.
#! e) Add commanline scheme. Tue 27 Dec 21:03, 2022.
#! f) Fix calculate values, using only average date [1:] data. Thue 5 Jan 00:59 pm, 2023.

import sys
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt

def usage():
    print("[INFO   ] usage: python NAMDEqStats.py -flog <FLOG> -plot <PLOT>")
    print("[INFO   ]        -f  <FLOG> : Is log file get from MD equilibration.")
    print("[INFO   ]        -p  <PLOT> : Is an binary 0 or 1 value. Where 0 is off and 1 is on, to generate plots.")
    print("[INFO   ]        -h  <H>    : Print it message.")


def optparser(argv=sys.argv[1:]):
    args = dict() 

    if '-f' in argv:
        ifile = argv.index('-f') + 1 
        args.update({'ifile' : argv[ifile]})

    else: 
        print("[INFO   ] Pass file to read.")
        usage()
        sys.exit(2)

    if '-p' in argv:
        plot = argv.index('-p') + 1
        args.update({'plot' : int(argv[plot])})

    else: 
        print("[INFO   ] Pass plot to write or not.")
        usage()
        sys.exit(2)

    if len(argv) < 1:
        usage()

    return args

#! * return parameters

params = optparser()

#! * set empty array.

pressavg    = []
tempavg     = []
volume      = []
temperature = []
pressure    = [] 
kinetics    = []
frame       = []
etotal      = []
gpressure   = []
gpressavg   = []


#! * read file. 
with open(params["ifile"]) as file: 
    for line in file: 
            if line.startswith("ENERGY") and float(line.split()[12]) != 0 : #! it only avoid use data from energy minimization. 
                frame.append(int(line.split()[1]))                                           #! ok 
                if float(line.split()[16]) != 0: pressure.append(float(line.split()[16]))    #! ok 
                if float(line.split()[12]) != 0: temperature.append(float(line.split()[12])) #! ok 
                if float(line.split()[18]) != 0: volume.append(float(line.split()[18]))      #! ok 
                if float(line.split()[15]) != 0: tempavg.append(float(line.split()[15]))     #! ok 
                if float(line.split()[19]) != 0: pressavg.append(float(line.split()[19]))    #! ok
                if float(line.split()[11]) != 0: etotal.append(float(line.split()[11]))      #! ok 
                if float(line.split()[17]) != 0:  gpressure.append(float(line.split()[17]))  #! ok
                if float(line.split()[20]) != 0: gpressavg.append(float(line.split()[20]))   #! ok
                if float(line.split()[10]) != 0: kinetics.append(float(line.split()[10]))    #! ok

            
ns = []
for i in frame:
    ns.append((float(i) * 0.002)/1000)

def mean(x):
    return np.average(x[1:])

def error(y):
    return np.std(y[1:])/np.sqrt(len(y[1:]))

def std(z):
    return np.std(z[1:])

print("Print Average from Data:")
print("==========================")
print("Analysis   : from {} to {} timestep with n={} data.".format(frame[1], frame[-1], len(frame)-1))
print("ETOTAL     : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(etotal), error(etotal), std(etotal))) 
print("KINETICS   : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(kinetics), error(kinetics), std(kinetics))) 
print("TEMP       : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(temperature), error(temperature), std(temperature))) 
print("TEMPAVG    : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(tempavg), error(tempavg), std(tempavg))) 
print("PRESS      : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(pressure), error(pressure), std(pressure))) 
print("PRESSAVG   : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(pressavg), error(pressavg), std(pressavg))) 
print("GPRESSURE  : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(gpressure), error(gpressure), std(gpressure))) 
print("GPRESSAVG  : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(gpressavg), error(gpressavg), std(gpressavg))) 
print("VOLUME     : {0:.2f} error: {1:.2f} std: {2:.2f}".format(mean(volume), error(volume), std(volume))) 

if params["plot"] == 1: 

    def save_data(time, name, observable):
        
        name = name
        df = pd.DataFrame({"time": ns, name: observable})
        df.to_csv("{}.csv".format(name), index=False) 
        print("[INFO   ] Writting {} data.".format(name)) 


    def plotting(observable):
        #for k in range(0, len(observable)):
        a = pd.read_csv(observable+".csv", sep=",")
        a = pd.DataFrame(a)
        #block_size = 1000 # 1 ns 
        #a_block = a.groupby(a.index // block_size)
        #block_mean = a_block.mean()

        fig, ax = plt.subplots()
        ax.plot(a["time"][::10], a[observable][::10],  linewidth = 0.5, color = "black")

        minit=a["time"].min()
        maxit=a["time"].max()
        minip=a[observable].min()
        maxip=a[observable].max()

        ax.set( xlabel="time/ns", ylabel=observable)
        ax.set_xlim(minit,maxit)
        ax.set_ylim(minip,maxip)
        ax.ticklabel_format(style='plain')
        ax.ticklabel_format(useOffset=False, style='plain')
        #ax.grid()

        fig.set_size_inches(9, 5)
        fig.tight_layout()
        title = observable+ ".png"
        plt.savefig(title)

    #if "__main__" == "__name__":

    observables = {"pressavg":pressavg, "tempavg":tempavg, "volume":volume, "temperature":temperature, 
                  "pressure":pressure, "kinetics":kinetics, "etotal":etotal, "gpressure":gpressure, "gpressavg":gpressavg}

    for key, obs in observables.items():
        save_data(time=ns, name=key, observable=obs)
        plotting(observable=key)
    print("[INFO   ] NAMD2EqStats done.")


