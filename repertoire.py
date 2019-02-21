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
import merge
import reformat

def mimp_finderParser(subparsers):
  mimp_finder_parser = subparsers.add_parser('mimp_finder', help='Mimp_finder finds TIR elements in a genome and records regions up/down stream of them')
  mimp_finder_parser.add_argument('-dir', '--directory', help='Directory containing all genomes to be examined', dest='directory_folder', type=str)
  mimp_finder_parser.add_argument('-c', '--cut_off', help='Max length of TIR elements, default=400', dest='cut_off', type=int, default=400)
  mimp_finder_parser.add_argument('-w', '--work', help="Working directory", dest='working', type=str)
  mimp_finder_parser.add_argument('-s', '--seed', help='A file containing example TIRs. TIRmite uses this to find similar TIRs throughout genomes.', dest='seed_mimp', type=str)
  mimp_finder_parser.add_argument('-d', '--distance', help='Distance up/downstream that is recorded, default=2500', dest='distance', type=int, default=2500)
  mimp_finder_parser.add_argument('-me', '--maxeval', help='Max e_value for TIRmite, default=1', dest='maxeval', type=float, default=1)
  mimp_finder_parser.add_argument('-md', '--maxdist', help='Max distance between TIRs in TIRmite, default=10000', dest='maxdist', type=int, default=10000)
  mimp_finder_parser.add_argument('-mc', '--mincov', help='Minimum coverage of found TIRs relative to reference TIRs needed to be considered, default=.9', dest='mincov', type=float, default=.9)
  mimp_finder_parser.add_argument('-f', '--force', help='Do not allow program to rewrite over old data', dest='force',action='store_false')


  return mimp_finder_parser

class mimp_finderCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = mimp_finder.mimp_finderApp()
    return app.start(args.force, args.directory_folder,args.working, args.seed_mimp, args.cut_off, args.distance, args.maxeval, args.maxdist, args.mincov)
    # maxeval=1, maxdist=10000, mincov=.9
def gene_finderParser(subparsers):
  gene_finder_parser = subparsers.add_parser('gene_finder', help='Gene_finder finds effectors / genes using given sequence data')
  gene_finder_parser.add_argument('-dir', '--directory', help='Directory containing all raw sequence data', dest='directory_folder', type=str)
  gene_finder_parser.add_argument('-o', '--out', help='Output directory', dest='output_dir', type=str)
  gene_finder_parser.add_argument('-min', '--min_prot_len', help='Minimum length for a protein in amino acids, default=15', dest='min_prot_len', type=int, default=15)
  gene_finder_parser.add_argument('-max', '--max_prot_len', help='Max length for a protein in amino acids, default=134', type=int, dest='max_prot_len', default=134)
  gene_finder_parser.add_argument('-d2m', '--max_d2m', help='Max distance between mimp TIR and the first amino acid, default=2500', dest='max_d2m', type=str, default=2500)
  gene_finder_parser.add_argument('-s', '--signalp', help='Directory of signalP', dest='SignalPpath', type=str, default='signalP')
  gene_finder_parser.add_argument('-sp', '--signalP_threshold', help='Minimum score that will result in a positive prediction of a signal peptide, default=.45', dest='SignalP_threshold', type=float, default= .45)
  gene_finder_parser.add_argument('-f', '--force', help='Do not allow program to rewrite over old data', dest='force',action='store_false')

  return gene_finder_parser

class gene_finderCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = gene_finder.gene_finderApp()
    return app.start(args.force, args.directory_folder, args.output_dir, args.SignalPpath, args.min_prot_len, args.max_d2m, args.max_prot_len, args.SignalP_threshold)

def pres_abs_varParser(subparsers):
  pres_abs_var_parser = subparsers.add_parser('pres_abs_var', help='Pres_abs_var makes a presence absence sheet using target sequences and genomes of interest')
  pres_abs_var_parser.add_argument('-q', '--queryfile', help='File containing genes of interest', dest='queryfile', type=str)
  pres_abs_var_parser.add_argument('-g', '--genome_folder', help='Folder containing all genomes to be examined', dest='genome_folder', type=str)
  pres_abs_var_parser.add_argument('-bd', '--blastdb', help='A folder where all blasted files can be stored to', dest='blastdatabasedir', type=str)
  pres_abs_var_parser.add_argument('-b', '--blastbin', help='BLAST bin directory ', dest='BLASTbindir', type=str, default='bin')
  pres_abs_var_parser.add_argument('-o', '--out', help='Output directory', dest='outputdir', type=str)
  pres_abs_var_parser.add_argument('-p', '--perc', help='The minimum percent of shared identity between two sequences needed for the two to be examined, default=80', type=float, dest='PERC_IDENTITY_THRESH', default=80)
  pres_abs_var_parser.add_argument('-d', '--database', help='Build blastdb, yes/no, default = yes', dest='buildblastdb', type=str, default='yes')
  pres_abs_var_parser.add_argument('-r', '--r_location', help='Path to R script', dest='r_location', type=str)
  pres_abs_var_parser.add_argument('-f', '--force', help='Do not allow program to rewrite over old data', dest='force',action='store_false')
  pres_abs_var_parser.add_argument('-s', '--show', help="Genes that are present/absent in all genomes will be shown", dest='show_all',action='store_true')

  return pres_abs_var_parser

class pres_abs_varCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = pres_abs_var.pres_abs_varApp()
    return app.start(args.force, args.show_all,args.queryfile, args.genome_folder, args.blastdatabasedir, args.BLASTbindir, args.outputdir, args.buildblastdb, args.r_location, args.PERC_IDENTITY_THRESH)


def clusterParser(subparsers):
  cluster_parser = subparsers.add_parser('cluster', help='Takes sequences given in an input file and uses BLAST to cluster them')
  cluster_parser.add_argument('-i', '--infile', help='File containing sequences to be clustered, such as the all_putative_genes_concatenated.fasta file generated from the gene finder program', dest='infile', type=str)
  cluster_parser.add_argument('-bd', '--blastdir', help='A folder where all blasted files can be stored to', dest='blastdatabasedir', type=str)
  cluster_parser.add_argument('-dc', '--leave', help='Sequences contain both a description and an ID', dest='leave_put_eff_identifiers_during_clustering',  action='store_true')
  cluster_parser.add_argument('-b', '--BLASTbindir', help='BLAST bin directory (i.e. /usr/local/bin)', dest='BLASTbindir', type=str, default='bin')
  cluster_parser.add_argument('-x', '--examine', help='Text file containing sequence/cluster headers. Every sequence that was clustered with a provided header will be outputed', dest='examin', type=str, default="")
  cluster_parser.add_argument('-e', '--e_value', help='The minimum e value for BLAST, default=.001', dest='E_VALUE_THRESH', type=float, default=.001)
  cluster_parser.add_argument('-p', '--percent', help='The minimum percent of shared identity between two sequences needed for the two to be clustered, default=90', type=int, dest='PERC_IDENTITY_THRESH', default=90.0)
  cluster_parser.add_argument('-lt', '--ll_thresh', help='Bare minimum of overlap needed for low length clustering, default=.25', dest='low_length_thresh', type=float, default=.50)
  cluster_parser.add_argument('-w', '--low_identity', help='Bare minimum of identity needed to for low identity clustering, default=50', dest='low_identity', type=int, default=70)
  cluster_parser.add_argument('-l', '--length_thresh', help='Minimum amount of overlap two sequences need to have to cluster, default=.9', dest='LENGTH_THRESH', type=float, default=.9)
  cluster_parser.add_argument('-cn', '--check_n', help='Do not remove sequences the have more than (-n) percentage of Ns', dest='check_n',action='store_true')
  cluster_parser.add_argument('-n', '--allowed', help='Max percentage of Ns that are allowed in a sequence file, default=0', dest='n_allowed', type=float, default=0)
  cluster_parser.add_argument('-f', '--force', help='Do not allow program to rewrite over old data', dest='force',action='store_false')
  cluster_parser.add_argument('-hi', '--hilc', help='Generate high identity, low coverage clusters', dest='hilc', action='store_true')
  cluster_parser.add_argument('-li', '--lihc', help='Generate low identity, high coverage clusters', dest='lihc', action='store_true')
  cluster_parser.add_argument('-j', '--expanded', help='Do not generate a file listing which sequences are in what cluster as well as basic cluster statistics', dest='expanded', action='store_false')
  cluster_parser.add_argument('-t', '--threads', help='Number of threads to run blast with, default=1', dest='threads', type=int, default=1)
  cluster_parser.add_argument('-nm', '--align', help='Each sequence can gnerate this number of blast results, default=250', dest='alignments', type=int, default=250)

################
################
  return cluster_parser

class clusterCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = cluster.clusterApp()
    return app.start(args.infile, args.blastdatabasedir, args.BLASTbindir,args.check_n, args.force, args.lihc, args.hilc, args.expanded, args.PERC_IDENTITY_THRESH,args.leave_put_eff_identifiers_during_clustering, args.E_VALUE_THRESH,args.examin, args.low_length_thresh, args.low_identity, args.LENGTH_THRESH, args.n_allowed, args.threads, args.alignments) #args.all_data)

def extractParser(subparsers):
  cluster_parser = subparsers.add_parser('extract', help='Input representative clusters and receive the elements of the cluster')
  cluster_parser.add_argument('-e', '--expand', help='An expanded_clusters.txt file generated from clustering', dest='expanded', type=str)
  cluster_parser.add_argument('-d', '--database', help='The file all_putative_effectors_concatenated.fasta generated from the mimp finder program', dest='database', type=str)
  cluster_parser.add_argument('-x', '--to_examine', help='A file containing all genes you wish to examine', dest='to_examine', type=str)
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
  pat_match_parser = subparsers.add_parser('pat_match', help='Searches tab-delimited presence absence table for specified pattern')
  pat_match_parser.add_argument('-o', '--outputfile', help='Output file name', dest='outputfile', type=str, default="pat_match_output.txt")
  pat_match_parser.add_argument('-i', '--inputfile', help='Input file name, pres / abs table generated from pres_abs_var', dest='inputfile', type=str)

  return pat_match_parser

class pat_matchCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = pat_match.pat_matchApp()
    return app.start(args.inputfile, args.outputfile)
   
def mergeParser(subparsers):
  merge_parser = subparsers.add_parser('merge', help='Will take a list of sequences and a list of clusters and merge the two together')
  merge_parser.add_argument('-c', '--cluster', help='Location of clusters', dest='clusters', type=str)
  merge_parser.add_argument('-i', '--infile', help='File containing sequences to merge with clusters', dest='infile', type=str)
  merge_parser.add_argument('-t', '--threads', help='Number of threads to run blast with, default=1', dest='threads', type=int, default=1)
  merge_parser.add_argument('-a', '--align', help='Each sequence can gnerate this number of blast hits, default=3', dest='alignments', type=int, default=3)
  merge_parser.add_argument('-b', '--BLASTbindir', help='BLAST bin directory (i.e. /usr/local/bin)', dest='BLASTbindir', type=str, default='bin')
  merge_parser.add_argument('-bd', '--blastdir', help='A folder where all blasted files can be stored to', dest='blastdatabasedir', type=str)
  merge_parser.add_argument('-e', '--e_value', help='The minimum e value for BLAST, default to .001', dest='E_VALUE_THRESH', type=float, default=.001)
  merge_parser.add_argument('-p', '--percent', help='The minimum percent of shared identity between two sequences needed for the two to be clustered, default=80.0', type=int, dest='PERC_IDENTITY_THRESH', default=80.0)
  merge_parser.add_argument('-l', '--length_thresh', help='Minimum amount of overlap two sequences need to have to cluster, default=.8', dest='LENGTH_THRESH', type=float, default=.8)
  merge_parser.add_argument('-f', '--force', help='Do not allow program to rewrite over old data', dest='force',action='store_false')
  merge_parser.add_argument('-ex', '--expanded', help='An expanded_clusters.txt file generated from clustering', dest='expanded', type=str)

  return merge

class mergeCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = merge.mergeApp()
    #infile, clusters, blastdatabasedir,BLASTbindir, identity=80, coverage=80,E_VALUE_THRESH=.001,threads=1,alignments=1
    return app.start(args.infile, args.clusters, args.expanded, args.blastdatabasedir, args.BLASTbindir, args.PERC_IDENTITY_THRESH, args.LENGTH_THRESH, args.E_VALUE_THRESH, args.threads, args.alignments, args.force)

def reformatParser(subparsers):
  reformat_parser = subparsers.add_parser('reformat', help="Reformats fasta file(s)'s headers so that they can be processed by the mimp_finder program")
  reformat_parser.add_argument('-i', '--infile', help='Reformat the singular file ___', default='', dest='infile', type=str)
  reformat_parser.add_argument('-id', '--identification', default='', help='Add this to the start of each sequence in each fasta file. If left blank, will put file name in header', dest='identification', type=str)
  reformat_parser.add_argument('-fo', '--folder', help='Will reformat all fasta files in specified folder',default='',  dest='folder', type=str)

  return reformat

class reformatCMD:
  
  def __init__(self):
    pass

  def execute(self, args):
    app = reformat.reformatApp()
    #infile, clusters, blastdatabasedir,BLASTbindir, identity=80, coverage=80,E_VALUE_THRESH=.001,threads=1,alignments=1
    return app.start(args.infile, args.folder, args.identification)

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Repertoire, a python program for clustering and analysing genetic data', add_help=True)
    subparsers = parser.add_subparsers(help='commands', dest='command')
    mimp_finderParser(subparsers)
    gene_finderParser(subparsers)
    pres_abs_varParser(subparsers)
    clusterParser(subparsers)
    extractParser(subparsers)
    pat_matchParser(subparsers)
    mergeParser(subparsers)
    reformatParser(subparsers)
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
    merge=mergeCMD()
    reformat=reformatCMD()
    commands = {'reformat':reformat,'merge':merge,'cluster': cluster, 'pres_abs_var': pres_abs_var, 'gene_finder': gene_finder, 'mimp_finder': mimp_finder, 'extract':extract, 'pat_match':pat_match}
    args = parseArgs()
    commands[args.command].execute(args)

if __name__ == '__main__':
    main()


