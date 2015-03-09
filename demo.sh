#!/bin/sh

# This script will compile cnF2freq using the Intel compiler (recommended, recent GCC is good, but Intel auto-vectorization is nice for performance)
g++ -O2 cnF2freq.cpp -ffast-math  -o cnF2freq -fopenmp -I boost_1_57_0 -g
icc -openmp -openmp-linkstatic -fast -ftz cnF2freq.cpp -o cnF2freq



# PlantImpute is using three demo files, a marker map, a pedigree and genotype data.
# Markers are simple top down cM locations. The pedigree is a list of individuals
# <ind id> <parent 1 id> <parent 2 id> [<generation>]
# Generation is 0 if missing. Parent ID 0 also means missing.
# Generation >=2 are those individuals that are actually processed.
# If generation 0 individuals are listed as parents, an implicit F1 generation
# is generated.
# A backcross can be generated by introducing a gen-1 individual in this list.
# Generation numbers larger than 2 mean implicit selfing, so the listed pedigree
# should be the ancestors for a generation 2 individual, and then further selfing
# is assumed from generation 2 down.
# The very last marker is a dummy, if not included, some results for the last marker
# will be wrong (sorry about that)
#
# Genotypes are specified as allele count 0,1,2. 9 for missing genotypes. An alternative
# format is <count-1>/<count-2> such as 5/0 for 5 reads of allele 1 and 0 reads for allele 2.
#
# During iterations (iteration count 10 in this example) most output is written to
# stdout. Pipe to file if this gets excessive.
# The final output file contains genotype probabilities for 11, 12, 21, 22 respectively,
# i.e. imputed individuals.
#

# g++ puts a lot of thread private storage on the stack
# This is large enough for at least some use cases. The defaults on some implementations will be too low.
export OMP_STACK_SIZE=256M
./cnF2freq demoplantimpute.map demoplantimpute.ped demoplantimpute.gen demooutput 10
