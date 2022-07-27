#! * Extract Temperature from NAMD
#! * @authors : Vega-Chozo K. & Ropón-Palacios G. 
#! * date : April 2, 2022. 
#! * e-mail : groponp@gmail.com 

import sys
ifile = sys.argv[1] 


#! * set empty array.

pressavg =[] 
tempavg = []
volume =[]
temperature =[]
pressure = [] 
frame = []


#! * read file. 
with open(ifile) as file: 
	for line in file: 
		if line.startswith("ENERGY"):
			frame.append(line.split()[1])
			pressure.append(line.split()[16])
			temperature.append(line.split()[12])
			volume.append(line.split()[18])
			tempavg.append(line.split()[15])
			pressavg.append(line.split()[19])
			
ns = []
for i in frame:
    ns.append((float(i) * 0.002)/1000)

#! * write data
import pandas as pd 

df = pd.DataFrame({"Time":ns, "Pressure": pressure})
df.to_csv("Time_pressure", index=False)
print("[Info ] Data was save")

df = pd.DataFrame({"Time":ns, "Temperature": temperature})
df.to_csv("Time_temperature", index=False)
print("[Info ] Data was save")

df = pd.DataFrame({"Time":ns, "Volume": volume})
df.to_csv("Time_volume", index=False)
print("[Info ] Data was save")

df = pd.DataFrame({"Time":ns, "Pressavg":pressavg})
df.to_csv("Time_Pressavg", index=False)
print("[Info ] Data was save")

df = pd.DataFrame({"Time":ns, "Tempavg": tempavg})
df.to_csv("Time_Tempavg", index=False)
print("[Info ] Data was save")

#plot 

import matplotlib.pyplot as plt

mols = ["pressure"]

for i in range(0,len(mols)):
    a = pd.read_csv('Time_'+mols[i], sep = ",")
    a = pd.DataFrame(a)
	
    fig, ax = plt.subplots()
    ax.plot(a["Time"], a["Pressure"],  linewidth = 1, color = "green")
    
    minit=a.Time.min()
    maxit=a.Time.max()
    
    minip=a.Pressure.min()
    maxip=a.Pressure.max()
       
    ax.set( xlabel='time/ns', ylabel='Pressure')
    p=ax.set_xlim(minit,maxit)
    k=ax.set_ylim(minip,maxip)
    print(" plot pressure vs time ")
    print(" el limite minimo y maximo del tiempo es ", p)
    print(" el limite minimo y maximo de la presion es ",k)
    ax.grid()

    fig.set_size_inches(9, 5)
    fig.tight_layout()
    title = 'time'+mols[i][:-4] + ".png"
    plt.savefig(title)
    plt.show() 
    
   


mols = ["temperature"]

for i in range(0,len(mols)):
    a = pd.read_csv('Time_'+mols[i], sep = ",")
    a = pd.DataFrame(a)

    fig, ax = plt.subplots()
    ax.plot(a["Time"], a["Temperature"],  linewidth = 1, color = "orange")
    
    minit=a.Time.min()
    maxit=a.Time.max()
    
    minitem=a.Temperature.min()
    maxitem=a.Temperature.max()

    ax.set( xlabel='time/ns', ylabel='Temperature/K')
    p1=ax.set_xlim(minit,maxit)
    k1=ax.set_ylim(minitem,maxitem)
    print(" plot temperature vs time ")
    print(" el limite minimo y maximo del time es ", p1)
    print(" el limite minimo y maximo de la temperature es ",k1)
    ax.grid()

    fig.set_size_inches(9, 5)
    fig.tight_layout()
    title = 'time'+mols[i][:-4] + ".png"
    plt.savefig(title)
    plt.show() 


mols = ["volume"]

for i in range(0,len(mols)):
    a = pd.read_csv('Time_'+mols[i], sep = ",")
    a = pd.DataFrame(a)

    fig, ax = plt.subplots()
    ax.plot(a["Time"], a["Volume"],  linewidth = 1, color = "black")
    
    minit=a.Time.min()
    maxit=a.Time.max()
    
    miniv=a.Volume.min()
    maxiv=a.Volume.max()
    
    ax.set( xlabel='time/ns', ylabel='Volume')
    p3=ax.set_xlim(minit,maxit)
    k3=ax.set_ylim(miniv,maxiv)
    print(" plot volume vs time ")
    print(" el limite minimo y maximo del time es ", p3)
    print(" el limite minimo y maximo del volume es ",k3)
    ax.grid()

    fig.set_size_inches(9, 5)
    fig.tight_layout()
    title = 'time'+mols[i][:-4] + ".png"
    plt.savefig(title)
    plt.show() 
    


mols = ["Tempavg"]

for i in range(0,len(mols)):
    a = pd.read_csv('Time_'+mols[i], sep = ",")
    a = pd.DataFrame(a)

    fig, ax = plt.subplots()
    ax.plot(a["Time"], a["Tempavg"],  linewidth = 1, color = "red")
    
    minit=a.Time.min()
    maxit=a.Time.max()
    
    minitemg=a.Tempavg.min()
    maxitemg=a.Tempavg.max()
    
    ax.set( xlabel='time/ns', ylabel='Tempavg/K')
    p4=ax.set_xlim(minit,maxit)
    k4=ax.set_ylim(minitemg,maxitemg)
    print(" plot tempavg vs time ")
    print(" el limite minimo y maximo del time es ", p4)
    print(" el limite minimo y maximo del tempavg es ",k4)
    ax.grid()

    fig.set_size_inches(9, 5)
    fig.tight_layout()
    title = 'time'+mols[i][:-4] + ".png"
    plt.savefig(title)
    plt.show() 


mols = ["Pressavg"]

for i in range(0,len(mols)):
    a = pd.read_csv('Time_'+mols[i], sep = ",")
    a = pd.DataFrame(a)

    fig, ax = plt.subplots()
    ax.plot(a["Time"], a["Pressavg"],  linewidth = 1, color = "red")
    
    minit=a.Time.min()
    maxit=a.Time.max()
    
    minipg=a.Pressavg.min()
    maxipg=a.Pressavg.max()
    
    ax.set( xlabel='Time/ns', ylabel='Pressavg/bar')
    p5=ax.set_xlim(minit,maxit)
    k5=ax.set_ylim(minipg,maxipg)
    print(" plot Pressavg vs time ")
    print(" el limite minimo y maximo del time es ", p5)
    print(" el limite minimo y maximo del pressavg es ",k5)
    ax.grid()

    fig.set_size_inches(9, 5)
    fig.tight_layout()
    title = 'time'+mols[i][:-4] + ".png"
    plt.savefig(title)
    plt.show() 
