#!/bin/sh

# This script will compile cnF2freq using the Intel compiler (recommended, recent GCC is good, but Intel auto-vectorization is nice for performance)
g++ -O2 cnF2freq.cpp -ffast-math  -o cnF2freq -fopenmp -I boost_1_57_0 -g
icc -openmp -openmp-linkstatic -fast -ftz cnF2freq.cpp -o cnF2freq



# cnF2freq is the run on two demo files, a marker map and a pedigree
# The very last marker is a dummy, if not included, some results for the last marker will be wrong (sorry about that)
# During iterations (iteration count 35 in this example) most output is written to stdout. Pipe to file if this gets excessive.
# Per-marker phasing, sureness and inferred genotypes are only written to stdout
# Final output to demooutput file will contain per-marker probabilities for the offspring individuals, expressed as:
# Parent 2 haplotype * 2 + parent 1 haplotype

# g++ puts a lot of thread private storage on the stack
# This is large enough for at least some use cases. The defaults on some implementations will be too low.
EXPORT OMP_STACK_SIZE=128M
./cnF2freq halfsibdemo.map halfsibdemo.ped demooutput 35

# In this example, the parent haplotypes are 1111 1121 (sire), and 1111 2222 (d1)
# This is correctly inferred, despite genotypes lacking for d1
# In this toy example, no recombinations happen of course, but the code can handle
# arbitrary recombination patterns in the offspring individuals and infer likely
# genotypes and haplotypes in parents based on that data.

# Note that this specific version of the code is tailored to the half-sib case,
# where parent 1 is of sex 1, parent 2 of sex 2 and only the parent 2s are to be
# statistically inferred.
