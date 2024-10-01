#!/bin/sh
#SBATCH -J domAddDMISel
#SBATCH --time=01:00:00
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o domAddDMISel-%j.out
#SBATCH -e domAddDMISel-%j.err

calculate() {
printf "%s\n" "$@" | bc -l;
}

#load modules
module load anaconda3/2022.10
conda activate incompSim

#######################
#EDIT THESE PARAMETERS#
#######################
parameters=domAddDMISel;
genes_file=genes_$parameters.txt
admixem_file=admixsimul_$parameters.txt

#make and empty directory so admixem starts making directories for each replicate that have the prefix number
mkdir $parameters
	
#Run job_file.slrm 1000 times
echo "starting at date on hostname"
	
#change i for number of replicates you want to simulate
for i in $(seq 100); do
(	
	echo "replicate $i"

	#make directory
	mkdir $parameters\_$i
		
	#enter directory
	cd $parameters\_$i

	#set variable to see if simulation passes or fails
	simPass=FALSE
	#repeat iteration while simulation doesnt successfully complete with all populations
	while [ $simPass == "FALSE" ];
	do

		#copy necessary files
		rsync --recursive --exclude=$parameters* .. .

		#perl -i.bak -lpe 'BEGIN { sub inc { my ($num) = @_; ++$num } } s/(RandomSeed=)(\d+)/$1 . (inc($2))/eg' $admixem_file
		#rm $admixem_file.bak

		./admixemp $admixem_file > /dev/null

		cd $parameters

		#check the sex ratios of the three pops
		if ! test -f Gen20_phenostats.txt; 
		then
			echo "Simulation did not finish, repeating replicate"
			simPass=FALSE
			cd ..
			rm -r $parameters
		elif [ $( awk 'FNR == 5 {print $3}' Gen20_phenostats.txt ) == 0[0] ] ||
		   [ $( awk 'FNR == 7 {print $3}' Gen20_phenostats.txt ) == 0[0] ] ||
		   [ $( awk 'FNR == 9 {print $3}' Gen20_phenostats.txt ) == 0[0] ]; 
		then
			echo "One or more populations went extinct before the simulation finished, repeating replicate"
			simPass=FALSE
			cd ..
			rm -r $parameters
		else 
			simPass=TRUE
			cd ..
		fi
	done
	
	#remove unnecesary files
	rm *.txt *.php admixemp *.bash
	
	#copy output from nested directory to current directory
	rsync --recursive ./$parameters/ .
	
	rm -r $parameters/
	
	rm ./Gen*\_phenotypes.txt
	rm ./Gen*\_markers.txt
	rm ./Gen*\_natselprobdump.txt
	rm ./Gen*\_phenostats.txt
	
	csplit -z -f parents Gen0_genes.txt '/id/' '{*}' > /dev/null

	for((j=1;j<=20;j++))
	do
		csplit -z -f hybrids$j Gen$j\_genes.txt '/id/' '{*}' > /dev/null	
	done 
	rm Gen*
	rm parents*
	
	cd ..

)

#rm slurm*
done
wait

echo "ended at date on hostname"
#exit 0

#Run phenotype_parser - plots mean/sd/cv and calculates proportion of overlap between hybrid and parent populations among all replicates.
Rscript domAddDMISel_genotype_parser3.R $parameters
rm -r $parameters\_{1..100}
rm -r $parameters
	
#Run violin_plots
#Rscript violin_plots_noss_nons.R 

#combine dataframe output for each initial allele freq into a single spreadsheet
#head -1 df_$parameters.csv > df_$parameters\_all.csv; tail -n +2 -q df_$parameters.* >> df_$parameters\_all.csv

mkdir output_$parameters
mv $parameters-* output_$parameters
