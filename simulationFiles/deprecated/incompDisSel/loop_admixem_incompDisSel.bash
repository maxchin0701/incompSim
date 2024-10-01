#!/bin/sh
#SBATCH -J incompDisSel
#SBATCH --time=01:00:00
#SBATCH -p RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o incompDisSel-%j.out
#SBATCH -e incompDisSel-%j.err

calculate() {
printf "%s\n" "$@" | bc -l;
}

#load modules
module load anaconda3/2022.10
conda activate incompSim

#######################
#EDIT THESE PARAMETERS#
#######################
parameters=incompDisSel;
genes_file=genes_$parameters.txt
admixem_file=admixsimul_$parameters.txt
sexsel_file=sexualsel_$parameters.txt
disCoeff=0

#make and empty directory so admixem starts making directories for each replicate that have the prefix number
mkdir $parameters

#repeat across selection increments
for((k=1; k<=5; k++)){

	#increment dis coefficient
	if [ $k == 5 ];
	then
		increment_value=1
	else 
		increment_value=2
	fi
	Rscript dis_increment_selection.R $sexsel_file $increment_value $disCoeff 
	let disCoeff="disCoeff+increment_value"
	echo $disCoeff

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
	Rscript incompDisSel_genotype_parser3.R $parameters $disCoeff
	rm -r $parameters\_{1..100}
	rm -r $parameters
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

