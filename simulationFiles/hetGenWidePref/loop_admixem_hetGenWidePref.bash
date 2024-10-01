#!/bin/sh
#SBATCH -J hetGenWidePref
#SBATCH --time=01:00:00
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o hetGenWidePref-%j.out
#SBATCH -e hetGenWidePref-%j.err

calculate() {
printf "%s\n" "$@" | bc -l;
}

#load modules
module load anaconda3/2022.10
conda activate incompSim

#######################
#EDIT THESE PARAMETERS#
#######################
parameters=hetGenWidePref;
genes_file=genes_$parameters.txt
admixem_file=admixsimul_$parameters.txt
sexsel_file=sexualsel_$parameters.txt
prefCoeff=0

#make and empty directory so admixem starts making directories for each replicate that have the prefix number
mkdir $parameters

#add empty output directory
mkdir ../../outputs/$parameters

#repeat across selection increments
for((k=1; k<=1; k++)){

	#increment dis coefficient
	if [ $k == 5 ];
	then
		increment_value=1
	else 
		increment_value=2
	fi
	Rscript hetGenWidePref_increment_selection.R $sexsel_file $increment_value $prefCoeff 
	let prefCoeff="prefCoeff+increment_value"
	echo $prefCoeff

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
		
		#copy necessary files
		rsync --recursive --exclude=$parameters* .. .

		#perl -i.bak -lpe 'BEGIN { sub inc { my ($num) = @_; ++$num } } s/(RandomSeed=)(\d+)/$1 . (inc($2))/eg' $admixem_file
		#rm $admixem_file.bak

		./admixemp $admixem_file > /dev/null
		
		#remove unnecesary files
		rm *.txt *.php admixemp *.bash
		
		#copy output from nested directory to current directory
		rsync --recursive ./$parameters/ .
		
		rm -r $parameters/
		
		rm ./Gen*\_markers.txt
		rm ./Gen*\_natselprobdump.txt
		rm ./Gen*\_phenostats.txt
		rm Gen0*
		
		cd ..

	)

	#rm slurm*
	done
	wait
		
	echo "ended at date on hostname"
	#exit 0
			
	#Run phenotype_parser - plots mean/sd/cv and calculates proportion of overlap between hybrid and parent populations among all replicates.
	Rscript hetGenWidePref_genotype_parser3.R $parameters $prefCoeff
	#rm -r $parameters\_{1..100}
	#rm -r $parameters
}

#Run violin_plots
#Rscript violin_plots_noss_nons.R 

#combine dataframe output for each initial allele freq into a single spreadsheet
#head -1 df_$parameters.csv > df_$parameters\_all.csv; tail -n +2 -q df_$parameters.* >> df_$parameters\_all.csv

mkdir output_$parameters
mv $parameters-* output_$parameters

#combine all
cd ../../outputs
echo "combining output files for different selection coefficients"
head -1 $parameters\_Fis_0.2.csv > $parameters\_Fis_All.csv; tail -n +2 -q $parameters\_Fis_0.* >> $parameters\_Fis_All.csv
head -1 $parameters\_Ho_0.2.csv > $parameters\_Ho_All.csv; tail -n +2 -q $parameters\_Ho_0.* >> $parameters\_Ho_All.csv
head -1 $parameters\_RawGeno_0.2.csv > $parameters\_RawGeno_All.csv; tail -n +2 -q $parameters\_RawGeno_0.* >> $parameters\_RawGeno_All.csv

