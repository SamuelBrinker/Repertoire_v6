import sys
import os
import argparse
#############
import pres_abs_var #03
import gene_finder #01
import cluster #02
import mimp_finder #0
import extract
import pat_match

def mimp_finderParser(subparsers):
  mimp_finder_parser = subparsers.add_parser('mimp_finder', help='mimp_finder finds TIR elements in a genome and records regions up/down stream of them')
  mimp_finder_parser.add_argument('-dir', '--directory', help='directory containing all genomes to be examined', dest='directory_folder', type=str)
  mimp_finder_parser.add_argument('-c', '--cut_off', help='max length of TIR elements, default=400', dest='cut_off', type=int, default=400)
  mimp_finder_parser.add_argument('-w', '--work', help="working directory", dest='working', type=str)
  mimp_finder_parser.add_argument('-s', '--seed', help='Seeder mimps for TIRmite, use MIMP_TIR.fasta file', dest='seed_mimp', type=str, default='MIMP_TIR.fasta')
  mimp_finder_parser.add_argument('-d', '--distance', help='distance up/downstream that is recorded', dest='distance', type=int, default=2500)
  mimp_finder_parser.add_argument('-me', '--maxeval', help='max e_value for TIRmite, default=1', dest='maxeval', type=int, default=1)
  mimp_finder_parser.add_argument('-md', '--maxdist', help='max distance between TIRs in TIRmite, default=10000', dest='maxdist', type=int, default=10000)
  mimp_finder_parser.add_argument('-mc', '--mincov', help='minimum coverage of found TIRs to reference TIRs, default=.9', dest='mincov', type=int, default=.9)
  mimp_finder_parser.add_argument('-f', '--force', help='allow program to rewrite over old data, default=True', dest='force', type=bool, default=True)


  return mimp_finder_parser

class mimp_finderCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = mimp_finder.mimp_finderApp()
    return app.start(args.directory_folder,args.working, args.seed_mimp, args.cut_off, args.distance, args.maxeval, args.maxdist, args.mincov,args.force)
    # maxeval=1, maxdist=10000, mincov=.9
def gene_finderParser(subparsers):
  gene_finder_parser = subparsers.add_parser('gene_finder', help='gene_finder finds effectors / genes using given sequence data')
  gene_finder_parser.add_argument('-dir', '--directory', help='directory containing all raw sequence data', dest='directory_folder', type=str)
  gene_finder_parser.add_argument('-o', '--out', help='output directory', dest='output_dir', type=str)
  gene_finder_parser.add_argument('-min', '--min_prot_len', help='minimum length for a protein in amino acids', dest='min_prot_len', type=int, default=15)
  gene_finder_parser.add_argument('-max', '--max_prot_len', help='max length for a protein in amino acids', type=int, dest='max_prot_len', default=134)
  gene_finder_parser.add_argument('-d2m', '--max_d2m', help='max distance between mimp TIR and the first amino acid', dest='max_d2m', type=str, default=2500)
  gene_finder_parser.add_argument('-s', '--signalp', help='directory of the signalp file', dest='SignalPpath', type=str, default='signalP')
  gene_finder_parser.add_argument('-sp', '--SignalP_threshold', help='minimum score that will result in a positive prediction of a signal peptide, default .45', dest='SignalP_threshold', type=float, default= .45)
  gene_finder_parser.add_argument('-f', '--force', help='allow program to rewrite over old data, default=True', dest='force', type=bool, default=True)

  return gene_finder_parser

class gene_finderCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = gene_finder.gene_finderApp()
    return app.start(args.directory_folder, args.output_dir, args.SignalPpath, args.min_prot_len, args.max_d2m, args.max_prot_len, args.SignalP_threshold, args.force)

def pres_abs_varParser(subparsers):
  pres_abs_var_parser = subparsers.add_parser('pres_abs_var', help='pres_abs_var makes a presence absence sheet using target sequences and genomes of interest')
  pres_abs_var_parser.add_argument('-q', '--queryfile', help='file containing genes of interest', dest='queryfile', type=str)
  pres_abs_var_parser.add_argument('-g', '--genome_folder', help='folder containing all genomes to be examined', dest='genome_folder', type=str)
  pres_abs_var_parser.add_argument('-bd', '--blastdb', help='a folder where all blasted files can be stored to', dest='blastdatabasedir', type=str)
  pres_abs_var_parser.add_argument('-b', '--blastbin', help='BLASTbindir (i.e. /usr/local/bin)', dest='BLASTbindir', type=str)
  pres_abs_var_parser.add_argument('-o', '--out', help='Output directory', dest='outputdir', type=str)
  pres_abs_var_parser.add_argument('-p', '--perc', help='The minimum percent of shared identity between two sequences needed for the two to be examined, defaults to 80', type=int, dest='PERC_IDENTITY_THRESH', default=80)
  pres_abs_var_parser.add_argument('-d', '--database', help='Build blastdb, yes/no', dest='buildblastdb', type=str, default='yes')
  pres_abs_var_parser.add_argument('-r', '--r_location', help='path to R script', dest='r_location', type=str)
  pres_abs_var_parser.add_argument('-f', '--force', help='allow program to rewrite over old data, default=True', dest='force', type=bool, default=True)
  pres_abs_var_parser.add_argument('-s', '--show', help="True/False, if False a gene is present/absent in all genomes won't be shown, default=False", dest='show_all', type=bool, default=False)

  return pres_abs_var_parser

class pres_abs_varCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = pres_abs_var.pres_abs_varApp()
    return app.start(args.queryfile, args.genome_folder, args.blastdatabasedir, args.BLASTbindir, args.outputdir, args.buildblastdb, args.r_location, args.PERC_IDENTITY_THRESH, args.force, args.show_all)


def clusterParser(subparsers):
  cluster_parser = subparsers.add_parser('cluster', help='takes sequences and clusters them')
  cluster_parser.add_argument('-i', '--infile', help='the file all_putative_genes_concatenated.fasta file generated from the gene finder program', dest='infile', type=str)
  cluster_parser.add_argument('-bd', '--blastdir', help='a folder where all blasted files can be stored to', dest='blastdatabasedir', type=str)
  cluster_parser.add_argument('-dc', '--leave', help='input TRUE if you are clustering a clustered file that does not contain any description, but only a .id, default=FALSE', dest='leave_put_eff_identifiers_during_clustering', type=str, default="TRUE")
  cluster_parser.add_argument('-b', '--BLASTbindir', help='BLASTbindir (i.e. /usr/local/bin)', dest='BLASTbindir', type=str)
  cluster_parser.add_argument('-x', '--examine', help='file containing clusters / sequences to examine further', dest='examin', type=str, default="")
  cluster_parser.add_argument('-e', '--e_value', help='the minimum e value for BLAST, default to .001', dest='E_VALUE_THRESH', type=float, default=.001)
  cluster_parser.add_argument('-p', '--percent', help='the minimum percent of shared identity between two sequences needed for the two to be clustered, defaults to 90', type=int, dest='PERC_IDENTITY_THRESH', default=90.0)
  cluster_parser.add_argument('-t', '--ll_thresh', help='bare minimum of overlap needed for clustering, default=.25', dest='low_length_thresh', type=float, default=.50)
  cluster_parser.add_argument('-w', '--low_identity', help='bare minimum of identity needed to for clustering, default=50', dest='low_identity', type=int, default=70)
  cluster_parser.add_argument('-l', '--length_thresh', help='minimum amount of overlap two sequences need to have to cluster, default =.9', dest='LENGTH_THRESH', type=float, default=.9)
  cluster_parser.add_argument('-cn', '--check_n', help='remove sequences the have more than x percentage of Ns, defaults to True', dest='check_n', type=bool, default=True)
  cluster_parser.add_argument('-n', '--allowed', help='max percentage of Ns that are allowed in a sequence file', dest='n_allowed', type=float, default=0)
  cluster_parser.add_argument('-f', '--force', help='allow program to rewrite over old data, default=True', dest='force', type=bool, default=True)

################
################
  return cluster_parser

class clusterCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = cluster.clusterApp()
    return app.start(args.infile, args.blastdatabasedir, args.BLASTbindir, args.PERC_IDENTITY_THRESH,args.leave_put_eff_identifiers_during_clustering, args.E_VALUE_THRESH,args.examin, args.low_length_thresh, args.low_identity, args.LENGTH_THRESH, args.n_allowed, args.check_n, args.force) #args.all_data)

def extractParser(subparsers):
  cluster_parser = subparsers.add_parser('extract', help='input representative clusters and receive the elements of the cluster')
  cluster_parser.add_argument('-e', '--expand', help='A expanded_clusters.txt file', dest='expanded', type=str)
  cluster_parser.add_argument('-d', '--database', help='the file all_putative_effectors_concatenated.fasta generated from the mimp finder program', dest='database', type=str)
  cluster_parser.add_argument('-x', '--to_examine', help='a file containing all genes you wish to examine', dest='to_examine', type=str)
################
################
  return cluster_parser

class extractCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = extract.extractApp()
    return app.start(args.database, args.expanded, args.to_examine )

def pat_matchParser(subparsers):
  pat_match_parser = subparsers.add_parser('pat_match', help='searches tab-delimited presence absence table for specified pattern')
  pat_match_parser.add_argument('-o', '--outputfile', help='output file name', dest='outputfile', type=str, default="pat_match_output.txt")
  pat_match_parser.add_argument('-i', '--inputfile', help='input file name, pres / abs table generated from pres_abs_var', dest='inputfile', type=str)

  return pat_match_parser

class pat_matchCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = pat_match.pat_matchApp()
    return app.start(args.inputfile, args.outputfile)
   

def parseArgs():
    parser = argparse.ArgumentParser(
        description='dbcAmplicons, a python package for preprocessing of massively multiplexed, dual barcoded Illumina Amplicons', add_help=True)
    subparsers = parser.add_subparsers(help='commands', dest='command')
    mimp_finderParser(subparsers)
    gene_finderParser(subparsers)
    pres_abs_varParser(subparsers)
    clusterParser(subparsers)
    extractParser(subparsers)
    pat_matchParser(subparsers)
    args = parser.parse_args()
    return args

def main():
    lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
    if lib_path not in sys.path:
        sys.path.insert(0, lib_path)
    pres_abs_var = pres_abs_varCMD()
    cluster = clusterCMD()
    gene_finder = gene_finderCMD()
    mimp_finder=mimp_finderCMD()
    extract=extractCMD()
    pat_match=pat_matchCMD()
    commands = {'cluster': cluster, 'pres_abs_var': pres_abs_var, 'gene_finder': gene_finder, 'mimp_finder': mimp_finder, 'extract':extract, 'pat_match':pat_match}
    args = parseArgs()
    commands[args.command].execute(args)

if __name__ == '__main__':
    main()


