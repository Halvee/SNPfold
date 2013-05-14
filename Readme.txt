NAME
  SNPfold - calculate the amount of structure change in an RNA given 
            specified SNPs
SYNOPSIS
  SNPfold.py [options] <sequence or seq file> <mutations>
DESCRIPTION
  SNPfold is a program that is built entirely off of the partition function calculation feature found in RNAfold, which utilizes an algorithm based on the work of J.S. McCaskill to calculate a matrix of base-pairing probabilities for a given sequence of RNA. Given a sequence of RNA, the RNAfold program will return a 'dot.ps' file that contains, among other things, a listing of all non-zero probabilities of base pairing for all nucleotides in the sequence.   SNPfold constructs matrices of base-pairing probabilities from these output files for multiple sequence variants of equal length and quantifies the structure change brought about by a given mutation or mutations, as well as the mutation's significance in terms of structure changing capability compared to all other possible mutations in the RNA.   SNPfold can return several output files, including partition function matrix text and graph files, files containing the basepairing probabilities for each given nucleotide, and graphic representations of average absolute structure change among multiple mutations in a given RNA.
  When running SNPfold, the sequence and mutations are required at bare minimum. The user can either paste the sequence directly into the command line command to run the program, or list the name of a two-line fasta file containing the sequence of interest.  Mutations that the user wishes to analyze must be entered in colon-delimited format, in the format of wildtype_nucleotide(position relative to 5' end)mutant_nucleotide (ex: A31U:G33C:C18A etc).   If the user wishes to analyze sequence variants involving combinations of SNPs, they can list the SNPs for a given variant separated by commas ((ex: A31U,G33C:C18A,G17C etc).
  Note: SNPfold is written for use on a Linux/Unix type of operating system.   Use of SNPfold in any instance requires the installation of the Numpy module for Python.   The production of various output graphs requires the installation of the Matplotlib module for Python.
OPTIONS
  -a, --accurate    Calculate accurate p-values for mutations indicated by user.   This option will tell SNPfold to calculate correlation coefficients for all possible single mutations in the RNA, and then find the rank of the mutation(s) in question.   This allows for an accurate assessment of a mutation's secondary structure-changing potential in comparison to all other single mutations in a given transcript.   Unfortunately, it can make the runtime go up substantially.   Take note that the calculation of accurate p-values go up on the order of N^3, where N represents the length of the RNA sequence in question.   It is not advisable to try to calculate accurate p-values for a sequence greater than 2000 nucleotides in length.   Also take note that by default, the input sequence, mutation correlation coefficients and accurate p-values, as well as all correlation coefficients calculated for every possible single mutation will be saved in a directory in the output folder.
  -h, --help    Display the contents of 'running.txt,' a help file found in the SNPfold source folder.
  -n, --nameDir= <dirName>    Allows user to assign a name to the output dir where output files will be placed.   Default setting is to assign a randomized directory name consisting of eight characters.
  -o, --output  Toggle different output options. '-o partition' or 'partit' will print the partition matrix to the command line and automatically saves the sequence to the output directory. '-o shannon' or 'shan' will calculate the Shannon Entropy value (bits) for the wild type sequence and any mutant sequences, and save those to a folder in the output directory.
  -s, --save	Save output files into a directory found in the 'output' folder.  This will allow user to save the raw calculation result files, including data files for the partition function matrices, base pairing probabilities, and correlation coefficients/p-values for mutations indicated.   If user wishes to make their own graphs of output data, they may use the files produced by this as their starting point.


REFERENCES
  M. Halvorsen, J. S. Martin, S. Broadaway, A. Laederach (2010) Disease-Associated Mutations that Alter the RNA Structural Ensemble.   
I.L.  Hofacker,  W.  Fontana,  P.F.  Stadler, S. Bonhoeffer, M. Tacker, P. Schuster (1994) Fast Folding and Comparison of RNA  Secondary  Structures.   Monatshefte  f. Chemie 125: 167-188
  M.  Zuker, P. Stiegler (1981) Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information, Nucl Acid Res 9: 133-148 
  J.S. McCaskill (1990) The equilibrium partition  function  and  base  pair  binding probabilities for RNA secondary structures, Biopolymers 29: 1105-1119 
  I.L. Hofacker & P.F. Stadler (2006) Memory Efficient Folding Algorithms for Circular RNA Secondary Structures, Bioinformatics (2006)
  A.F. BompfA1/4newerer, R. Backofen, S.H. Bernhart, J. Hertel, I.L.  Hofacker,  P.F. Stadler,  S.  Will  (2007) "Variations on {RNA} Folding and Alignment: Lessons from Benasque" J. Math. Biol.
VERSION
This man page documents SNPfold version 1.2.
AUTHORS
Matt Halvorsen, Sam Broadaway, Josh Martin, Chas Kissick, Alain Laederach
BUGS
None known (yet).