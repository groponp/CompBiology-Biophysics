#!/bin/bash 
#; It Script perform segment MD on HPC with time limited. 
#; __authors  = Ropón-Palacios G., Ramirez-Olivos G. and Gervacio-Villareal A.
#; __date__   = 20 May, 2022. 
#; __e-mail__ = groponp@gmail.com  


#; All code bewlow is sytaxis for HPC with SLRUM gestor. 
#; ========================================================
#SBATCH --account=def-nike-ab
#SBATCH --mail-user=icamps@gmail.com
#SBATCH --mail-type=ALL

# Allocate 10 CPUs and 1000M RAM per CPU for 5 days
#SBATCH --time=5-0:0
#SBATCH -c10 --mem-per-cpu=1000
#SBATCH --gres=gpu:v100:1
#SBATCH --job-name="viro_wt"

# Load namd-multicore module
#module load StdEnv/2020 intel/2020.1.217 namd-multicore
#module load StdEnv/2020 cuda/11.0 namd-multicore/2.14

SLURM_CPUS_PER_TASK=20

# EQ-GamD
cd 04_eq_gamd/
count_state=$(ls *.state | wc -l)

if [[ $count_state = "0" ]]
then
	namd2 +p${SLURM_CPUS_PER_TASK} +devices 0 +idlepoll  md_eq_gamd.namd > md_eq_gamd.out
	
	cp gamd-eq-wrap.restart.coor ../05_gamd/
	cp gamd-eq-wrap.restart.xsc ../05_gamd/
	cp gamd-eq-wrap.restart.gamd ../05_gamd/
	cp gamd-eq-wrap.colvars.state ../05_gamd/	
fi

cd ../

#Production-GamD
cd 05_gamd/
count_dcd=$(ls *.dcd | wc -l)
count_out=$(ls prod_*.out | wc -l)
echo "No se asuste por este aviso"

if [[ $count_dcd = "0" ]]
then
	start_md=1
	sed "s/md_${count_out}/omd/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
	sed "s/md_$((count_out-1))/gamd-eq-wrap/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
	sed "s/set restart_inicial 0/set restart_inicial 1/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
        sed "s/set restart_continuar 1/set restart_continuar 0/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
	
else
	sed "s/md_${count_out}/md_${count_dcd}/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
	sed "s/md_$((count_out-1))/md_$((count_dcd-1))/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
	namd2 +p${SLURM_CPUS_PER_TASK} +devices 0 +idlepoll  prod.namd > prod_${count_dcd}.out
	let start_md=$((count_dcd+1));
fi

steps_initial=500000000
md_steps=$(awk -v var=run '{if( $1 == var) {print $2}}' prod.namd)
dt=$(awk -v var=timestep '{if( $1 == var) {print $2/1000}}' prod.namd)

time_ns=$(echo "scale=3; $steps_initial*$dt/1000" | bc)
echo "Time in ns is :" $(echo "scale=3; $time_ns " | bc)

save_each=2 #change it by split nanoseconsd to run. example 1 eq 1ns 
n_iterations=$(echo "scale=0; $time_ns/$save_each" | bc)
echo "Number of iterations is: $n_iterations"

namd_steps=$(echo "scale=0; $save_each*1000/$dt " | bc) 

sed "s/$md_steps/$namd_steps/g" prod.namd > tmp.namd
rm prod.namd 
mv tmp.namd prod.namd 

while [ $start_md -le $n_iterations ]; do
	prev_i=$((start_md-1))

	if [[ $start_md -eq 1  ]]
	then
		sed "s/omd/md_${start_md}/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                namd2 +p${SLURM_CPUS_PER_TASK} +devices 0 +idlepoll  prod.namd > prod_${start_md}.out
        	
	elif [[ $start_md -eq 2 ]]
 	then
		sed "s/md_${prev_i}/md_${start_md}/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                sed "s/gamd-eq-wrap/md_${prev_i}/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                sed "s/set restart_inicial 1/set restart_inicial 0/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                sed "s/set restart_continuar 0/set restart_continuar 1/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                namd2 +p${SLURM_CPUS_PER_TASK} +devices 0 +idlepoll  prod.namd > prod_${start_md}.out
                
	else
		sed "s/md_${prev_i}/md_${start_md}/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                sed "s/md_$((start_md-2))/md_$((start_md-1))/g" prod.namd > tmp.namd && rm prod.namd && mv tmp.namd prod.namd
                namd2 +p${SLURM_CPUS_PER_TASK} +devices 0 +idlepoll  prod.namd > prod_${start_md}.out
	fi
	let start_md=$((start_md+1));
done
cd ../




