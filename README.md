# polyphase
Phases the genotyped markers for polyploids.

Installation: This package requires the installation of boost(http://www.boost.org/), DAKOTA(https://dakota.sandia.gov/sites/default/files/docs/6.0/html-ref/index.html), TBB (https://www.threadingbuildingblocks.org/) and 
GNU multi precision libraries (http://www.mpfr.org/).
Description:
This program takes unphased haplotypes of genotyped data and the physical positions of the markers as input and returns the phased haplotypes as the output. 


Input files:
Parameter file: This file allows user to modify parameters like number of clusters. 
Data file: This provides the details of ploidy, number of individuals, number of markers and the unphased genotypes of each individual.

output files:
Phase_clusters.txt: This displays the array of most likely ancestral cluster, responsible for the current alleles.
Phased_haplo.txt:  This displays the phased data of all the individuals.
