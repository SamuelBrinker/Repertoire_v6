# Repertoire
Repertiore is a local BLAST gene clustering program which builds on the works of Peter van Dam's FoEC program. Currently a work in progress, it is being designed out of the Gordon Lab at UC Davis's Plant Pathology Department.

It is made in python 2.7 and requires the following dependencies:

Biopython

https://github.com/biopython/biopython

TIRmite 

https://github.com/Adamtaranto/TIRmite

Signalp

http://www.cbs.dtu.dk/services/SignalP/

locla BLAST, R


Example: 

mimp_finder:

Mimp finder finds TIR elements in a genome and records regions up/down stream of them

-dir directory containing all genomes to be examined

-c max length of TIR elements that are recorded, the default value for this is 400 bp

-w working directory for the program

-s seeder mimps for TIRmite. This is a file containing known, example mimps for TIRmite to base its search on

-d distance of the region up/downstream  of the mimp that is recorded that is recorded. This is the region that will be searched for effectors / genes

-f do not allow program to rewrite over old data





gene_finder:

Gene finder finds genes or effectors using given sequence file. Genes must be entirely located within the file or they will not be included. This can be used with the up/downstream region recorded from mimp finder. 

-dir directory containing all raw sequence data. This would be the folder containing the up/downstream files.

-o output directory

-min minimum length for a protein in amino acids pre removal of the signal peptide. Defaults to 10

-max max length for a protein in amino acids pre removal of the signal peptide. Defaults to 134

-d2m max distance between the start of the sequence and the first amino acid. Defaults to 2500

-s directory of the signalp file

-sp minimum score that will result in a positive prediction of a signal peptide. Defaults to .45

-f do not allow program to rewrite over old data





cluster:

Cluster takes sequences and groups them together based on coverage and identity scores generated from a self blast. A high and low version of both scores can be specified. High quality clusters (h_i_h_c) are comprised of sequences that meet the the high end requirements for both coverage and identity. l_i_h_c and h_i_l_c are clusters formed when genes meet one, but not both of the higher end parameters. Singletons are genes that were unable to be grouped in h_i_h_c clusters. 

-i The fasta file containing data to be clustered. For example, the file all_putative_genes_concatenated.fasta file generated from the gene finder program

-bd a folder where all blasted files can be stored to

-dc Data is a clustered file that does not contain any description, but only a .id. Defaults to false

-b BLASTbindir (i.e. /usr/local/bin)

-x a fasta file containing clusters / sequences to examine further. Optional

-e the minimum e value for BLAST. Default to .001

-p the minimum percent of shared identity between two sequences needed for the two to be clustered. Defaults to 90

-t bare minimum of overlap needed for clustering. Defaults to .25

-w bare minimum of identity needed to for clustering. Defaults to 50

-l minimum amount of overlap two sequences need to have to cluster. Defaults to .9

-cn do not remove sequences the have more than x percentage of Ns

-n max percentage of Ns that are allowed in a sequence file. Defaults to 0

-f do not allow program to rewrite over old data

-hi generates high identity, low coverage clusters

-li generates low identity, high coverage clusters

-j generates expanded clusters

-t number of threads to run blast with

-nm each sequence can gnerate this number of blast hits. Default is 250 hits


extract:

Extract allows the user to submit the header of a sequence of interest and retrieve all of the sequences of all other genes it is clustered with.

-e the expanded_clusters.txt file generated from cluster

-d the file containing all of the preclustered genes. Ex the file all_putative_effectors_concatenated.fasta generated from the mimp finder program.

-x a file containing all genes you wish to examine 




Merge:

Merge can be used to sort sequences into predtermined clusters. Users submit the sequences to be merged, all sequences present in the clusters, and the expanded_clusters file and returns the clusters with the sequences sorted in, file showing how many sequences were sorted into each cluster, a list of sequences that could not be merged, and a tab deliminated summary file in which the first column contains the ID of the unknown sequence, the second column has ID of the cluster representative, the third column is % identity, the fourth column is the alignment length, the fifth column is % coverage query/subject, and the sixth column is % coverage subject/query. 

-c location of cluster sequences

-i file containing sequences to merge in with clusters

-t number of threads to run blast with

-a each sequence can gnerate this number of blast hits. Default is 3 hits

-b BLASTbindir

-bd a folder where all blasted files can be stored to

-e the minimum e value for BLAST, default to .001

-p the minimum percent of shared identity between two sequences needed for the two to be clustered, defaults to 80

-l minimum amount of overlap two sequences need to have to cluster, default =.8

-f do not allow program to rewrite over old data

-ex location of the expanded file for the cluster




pres_abs_var:

Presence absence variation takes a fasta file containing genes of interest and directory full of genomes of interest. It will search for the presence of each of genes inside everyone of the genomes and form table of the results. Black indicates the presence of a given gene.

-q a fasta file containing all genes of interest. For example, the high_quality_clustered_genes.fasta file from cluster

-g a folder containing all genomes to be examined

-bd a folder where all blasted files can be stored to

-b BLASTbindir (i.e. /usr/local/bin)

-o output directory

-p the minimum percent of shared identity between two sequences needed for the two to be examined. Defaults to 80

-d Build a blast db of the genome, yes/no. Defaults to yes

-r path to R script

-f do not allow program to rewrite over old data

-s gene is present/absent in all genomes will be shown




pat_match:

Pat_match takes a modified presence / absences table and will search it a specified pattern. The table must be formated in the following fashion

                pattern gene_1   gene_2     gene_3    gene_4
        
        genome_1  N       1         0         1         1

        genome_2  1       1         1         0         1

        genome_3  1       0         1         1         1

        genome_4  0       0         0         0         0

        genome_5  N       1         0         1         0

        genome_6  1       1         0         1         1

A '1' or '0' indicates that a gene must be present/absent in a given genome while an 'N' will ignore the geome. This pattern will return gene_4 

-i input file name, pres / abs table generated from pres_abs_var

-o output file name. If the file name already exists or no name is entered it will default to pat_match_output.txt
