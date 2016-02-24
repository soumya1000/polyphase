# polyphase
This repository has two methods, 
 1. A C++ phasing program which phases genotyped markers for polyploids. Accepts unphased haplotypes of genotyped data and the physical positions of the markers as input and returns the phased haplotypes as the output. 
 2. A method that validates the out put of a phasing method against true haplotypes( for simulated data sets).
 Computes the switch error: The number switches required to match the haplotypes of phased out to those of true haplotypes. The lower this value, the better is the phasing method.

 1.Phasing method:
    1. This program requires the pre-installation of boost(http://www.boost.org/), DAKOTA(https://dakota.sandia.gov/sites/default/files/docs/6.0/html-ref/index.html), TBB (https://www.threadingbuildingblocks.org/) and 
	GNU multi precision libraries (http://www.mpfr.org/).
    2. Command : ./init Input_datafile
    
    Input files:
	Poly_phase_params.txt: This file allows user to modify parameters like number of clusters. This is already provided and does not need modifications in usual scenarios.
	Inputp_data file: This provides the details of ploidy, number of individuals, number of markers and the unphased genotypes of each individual.

    Output files:
	Phase_clusters.txt: This displays the array of most likely ancestral cluster, responsible for the current alleles.
	Phased_haplo.txt:  This displays the phased data of all the individuals.
    
 2. Phasing_ validation    
      1. Requires gtools package
      2. Command :  Rscript switch_error.R Phased_output.txt true_data.dat n_ploidy  n_ind n_mark, 
         where Phased_output.txt  is output from the phasing method and true_data.dat is a file with true values of haplotypes.
      3. The output  is a csv file with switch error for every individual, as shown in the example file.
