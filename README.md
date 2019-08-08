CVrep
Package developed by Smaragda Tsairidou (Tsairidou et al. 2014, Tsairidou et al. 2019)

This package can be used to generate low density SNP panels and to perform multiple cross-validation iterations for the calculation of genomic prediction accuracies using genomic data. The package contains 3 functions:

1.	“CreateLDdata” 

Description:
Selects SNPs to create low density SNP panels. SNPs can be selected folowing different methods and a number of SNP panel replicates can be generated for each density (see Options).

Usage:
CreateLDdata(filename="scottish_genotypes3.tped", nreps=2, method=2, denstart=200, denend=500, step=100)
CreateLDdata(filename="scottish_genotypes3.tped", nreps=2, method=3, denstart=200, denend=500, step=100, gen_len=2240.19, chr_len_file="Chr_len.txt")
CreateLDdata(filename="scottish_genotypes3.tped", method=4, stepstart=5000000, stepend=9000000, step=2000000)

Requirements:
The input file is a plink '.tped' data file. 
For method 3, “chr_len_file” is an input '.txt' file with two columns, where column one is the chromosome, column 2 is the chromosome name, and column three is the chromosome length. Chromosome length can be given either in Mbp or centimorgans. 

Options:
If “method=1” (default), then SNPs are selected by random sampling across the entire genome, where the SNPs selected in the different replicates may overlap by chance.
If “method=2”, then LD SNP panels are created with random sampling across the genome, but replicates are non-overlapping. 
If “method=3”, then SNPs are sampled randomly within each chromosome but proportionally to the chromosome length. This method also generates a “snp_count.txt” file where the exact number of SNPs sampled from each chromosome is recorded.  
If “method=4”, then SNPs are selected based on their physical (or genetic) distance for pre-defined step sizes.  For example, if stepstart=1000000, stepend=11000000, and step=2000000, then it will generate datasets where the gap between selected SNPs will be 1,3,5 Mbp etc until 11Mbp. This method generates no replicates. This method always selects the first and the last SNP on each chromosome. This method also generates a “snp_count.txt” file where the exact number of SNPs sampled from each chromosome is recorded.

Other arguments:
“nreps” is the number of sampling replicates for each density. 
“denstart” and “denend” are the starting (lowest) and the highest densities that will sample for as chosen by the user. For example, if denstart=100, denend=1000, and step=100, then the function samples densities of 100, 200, 300, 400 etc until 1000 SNPs. 
“step” is either the step by which the densities increase for method 3 (see “denstart” and “denend” above), or the SNP distance step size for method 4.
“gen_len” is the genome length if using method 3. 

Output:
Generates plink format '.tped' files with SNP genotypes for each density (step) and replicate (depends on the method used, see Options).
Generates a SNP count file with the step used and the number of SNPs selected (see file “snp_count.txt”) (methods 3 and 4).

Examples:
CreateLDdata(filename="genotypes.tped", nreps=3, method=1, denstart=200, denend=500, step=100)
CreateLDdata(filename="genotypes.tped", nreps=3, method=3, denstart=200, denend=500, step=100, gen_len=2240.19, chr_len_file="Chr_len.txt")
CreateLDdata(filename="genotypes.tped", method=4, stepstart=5000000, stepend=9000000, step=2000000)
