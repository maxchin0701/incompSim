# LTREB Incompatibility Simulations

Forward time simulations in Admixem investigating interactions between incompatibility selection and various models of mate choice.

####Incompatibility model

A single nuclear incompatibility locus is modeled, designed to mirror the mito-nuclear incompatibility empirically observed in the *Xiphophorus malinche* - *X. birchmanni* system. Individuals which are homozygous for *X. birchmanni* ancestry at this locus experience fitness reductions with S=0.91.

**Mate Choice Models**

hetGenWidePref: Preference for genome wide heterozygosity, calculated based on 24 loci distributed across the genome (1 per chrmoosome). Preference strength describes the P(mating) difference between completely heterozygous and homozygous individuals, and preference increases linearly with heterozygosity.

malGenWidePref: Preference for genome wide malinche, calculated based on 24 loci distributed across the genome (1 per chrmoosome). Preference strength describes the P(mating) difference between completely heterozygous and homozygous individuals, and preference increases linearly with heterozygosity.

Each mate choice model is run with preference strengths of 0.2, 0.4, 0.6, 0.8, and 0.91.

**Directory structure**

simulation_files: raw simulation files for running in Admixem, organized by mate choice model. To run 100 iterations of each preference model, run the loop file as a bash script or slurm job.

outputs: processed outputs for simulations, also organized by mate choice model into subdirectories. Under each subdirectory, you will find four .csv files. 

      RawGeno: Raw genotypes for all individuals
      Phenotypes: Phenotypes (genome wide percent heterozygosity and malinche ancestry) for all individuals
      Fis: Fis statistics for each replicate-population-generation combination, calculated based solely on the incompatibility locus
      Ho: Heterozygosity for each replicate-population-generation combination, calculated based solely on the incompatibility locus


