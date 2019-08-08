CVrep
Package developed by Smaragda Tsairidou (Tsairidou et al. 2014, Tsairidou et al. 2019)

This package can be used to generate low density SNP panels and to perform multiple cross-validation iterations for the calculation of genomic prediction accuracies using genomic data. The package contains 3 functions:

1.	“CreateLDdata” 

Description:
Selects SNPs to create low density SNP panels. SNPs can be selected folowing different methods and a number of SNP panel replicates can be generated for each density (see Options).

Usage:
CreateLDdata(filename, nreps, method, denstart, denend, step, gen_len, chr_len_file)

CreateLDdata(filename, method, stepstart, stepend, step)

Requirements:
The input file is a plink '.tped' data file. 
For method 3, “chr_len_file” is an input '.txt' file with two columns, where column one is the chromosome, column two is the chromosome name, and column three is the chromosome length. Chromosome length can be given either in Mbp or centimorgans. 

Options:
If “method=1” then SNPs are selected by random sampling across the entire genome, where the SNPs selected in the different replicates may overlap by chance.
If “method=2”, then low density SNP panels are created with random sampling across the genome, but replicates are non-overlapping. 
If “method=3”, then SNPs are sampled randomly within each chromosome but proportionally to the chromosome length. This method also generates a “snp_count.txt” file where the exact number of SNPs sampled from each chromosome is recorded.  
If “method=4”, then SNPs are selected based on their physical (or genetic) distance for pre-defined step sizes. For example, if stepstart=1000000, stepend=11000000, and step=2000000, then it will generate datasets where the gap between selected SNPs will be 1,3,5 Mbp etc until 11Mbp. This method generates no replicates. This method always selects the first and the last SNP on each chromosome. This method also generates a “snp_count.txt” file where the exact number of SNPs sampled from each chromosome is recorded.

Other arguments:
“nreps” is the number of sampling replicates for each density. 
“denstart” and “denend” are the starting (lowest) and the highest densities that will sample for, as chosen by the user. For example, if denstart=100, denend=1000, and step=100, then the function samples densities of 100, 200, 300, 400 etc. until 1000 SNPs. 
“step” is either the step by which the densities increase for method 3 (see “denstart” and “denend” above), or the SNP distance step size for method 4.
“gen_len” is the genome length if using method 3. 

Output:
Generates plink format '.tped' files with SNP genotypes for each density (step) and replicate (depends on the method used, see Options).
Generates a SNP count file with the step used and the number of SNPs selected (see file “snp_count.txt”) (methods 3 and 4).

Examples:
CreateLDdata(filename="genotypes.tped", nreps=3, method=1, denstart=200, denend=500, step=100)
CreateLDdata(filename="genotypes.tped", nreps=3, method=3, denstart=200, denend=500, step=100, gen_len=2240.19, chr_len_file="Chr_len.txt")
CreateLDdata(filename="genotypes.tped", method=4, stepstart=5000000, stepend=9000000, step=2000000)


2. “CreateCVgroups” 

Description:
Randomises individuals and partitions into cross-validation groups. It repeats the process for a number of times, where each time different individuals are allocated into different cross-validation groups.

Usage:
createCVgroups(datafile, ngroups, nreps, multiPhenoFiles=F, PhenoHeader=T)
 
Output:
It saves phenotype output datafiles (see Options) where the last column shows the allocated group number, and saves R workspace to be used by function “calcAccur” (see below). 

Options: 
If “multiPhenoFiles=F” (default), a single phenotype file is produced for each replicate where the last column indicates the allocated cross-validation group number (phenotypes_grp_nrep#.txt). 

If “multiPhenoFiles=T”, then multiple phenotype files are created, as many as the cross-validation folds. For example, if it is a five-fold cross-validation (default), then 5 files are created for each repeat. Those files contain the training sets so that from each file one group is missing (the test set which corresponds to the group number shown in file name). The last column in each file shows the allocated cross-validation group number. 

Other arguments:
“nreps” is the number of cross-validation replicates, i.e. the number of 5-fold group assignment iterations, where each time, different individuals are allocated into different cross-validation groups.
If “PhenoHeader=T” (default) output phenotypic files have headers (set by the user ). 

Requirements: 
In the phenotypic data file the first column must be the recoded id so that it starts from 1 and increases to the total number of individuals (rows) in the datafile. 
In the phenotypic data file the second column must contain the phenotype trait value, i.e. the phenotype with which we will calculate the correlation (and accuracy) with the predicted breeding values (see function “calcAccur” below). 

Example:
createCVgroups("Phenotypes_for_CV.dat", nreps=10, multiPhenoFiles=F, PhenoHeader=T)

 
3.	“CalcAccur” 

Description:
It calculates correlations between the trait value and the predicted breeding values, using the output files from ASReml analysis. It calculates the genomic prediction accuracies (see Options). Creates ASReml status report file.

Usage:
calcAccur(ngroups, nreps, skipnum, genomichh=T, h2)

Options:
If “genomichh=F” (default), it uses the heritability from each cross-validation fold to calculate the accuracy
If “genomichh=T”, the user also needs to provide the genomic heritability (argument “h2”), so that the accuracy is calculated  by dividing with the square root of the genomic heritability. 

Other arguments:
“Groups” is the number of cross validation groups
“skipnum” is the number of lines to skip when reading the ASReml “.sln” files (depends on the number of fixed effects in the model and the number of their levels in the data). 

Requirements: 
It loads the R workspace saved from function “CreateCVgroups”.
It uses ASReml output files .asr, .sln, and .pvc. Those files must be named as: “cvgroup” fllowed by the group number, followed by "_nrep", the replicate number, and the type of file (for example, “cvgroup5_nrep1.asr”).

Output:
It generates a file where it records the exist status of ASReml, so that the user can check convergence for each cross-validation fold and each repeat (see file “asrem_status_report.txt”).
Generates file with correlations, heritabilities and accuracies, for each cross-validation fold and means for each repeat (file “CVResults.txt”).
Generates file only with mean correlation, heritability and accuracy for each repeat (file “rr.txt”). 

Example:
calcAccur(ngroups=5, nreps=10, skipnum=4, genomichh=F)


