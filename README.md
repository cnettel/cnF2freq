# cnF2freq

cnF2freq is a tool to compute genotype and haplotype probabilities, including imputing missing genotypes. The code has existed in multiple versions, but starting in February 2015 versioning is done on github. The first fork introduced here is the PlantImpute branch.

The code is more or less a monolithic beast of a C++ file, but the intent is to clean that up, eventually. Compile this file with whatever optimization settings suit you. Auto-vectorization and strict aliasing rules can do a lot to improve performance of the code.

For any academic use, in addition to following the BSD license, citing the following references is encouraged:
cnf2freq: Efficient determination of genotype and haplotype probabilities in outbred populations using markov models
C Nettelblad, S Holmgren, L Crooks, Ö Carlborg - Bioinformatics and Computational Biology, 2009

Inferring haplotypes and parental genotypes in larger full sib-ships and other pedigrees with missing or erroneous genotype data
C Nettelblad - BMC Genetics, 2012

Specific versions and uses of the code might have other appropriate publications.