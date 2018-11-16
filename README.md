# Repertoire
Repertiore is a gene clustering program based on the works of Peter van Dam's FoEC program. 
It is made in python 2.7 and requires the following dependencies:
Biopython
https://github.com/biopython/biopython
TIRmite 
https://github.com/Adamtaranto/TIRmite
Signalp
http://www.cbs.dtu.dk/services/SignalP/

Example: 

mimp_finder:
Mimp finder finds TIR elements in a genome and records regions up/down stream of them

-dir directory containing all genomes to be examined
-c max length of TIR elements that are recorded, the default value for this is 400 bp
-w working directory for the program
-s seeder mimps for TIRmite. This is a file containing known, example mimps for TIRmite to base its search on
-d distance of the region up/downstream  of the mimp that is recorded that is recorded


gene_finder:
Gene finder finds genes or effectors using given sequence data. This can be used with the up/downstream region recorded from mimp finder.

-dir directory containing all raw sequence data
-o output directory
-min minimum length for a protein in amino acids. Defaults to 10
-max max length for a protein in amino acids. Defaults to 134
-d2m max distance between the start of the sequence and the first amino acid. Defaults to 2500
-s location of the signalp file
-sp minimum score that will result in a positive prediction of a signal peptide. Defaults to .45
