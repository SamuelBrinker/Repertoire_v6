
import sys, os, re, glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

def Translator(infile, datahandler2):
	datahandler_list2 = []
	for seq_record in SeqIO.parse(infile, 'fasta', IUPAC.ambiguous_dna):
		for i in range(3): 				#(i=0 - i=1 - i=2)
			s = seq_record.seq[i:]   #checks each open reading frame of each gene in the fasta file
			while len(s)%3 != 0: # add N to the end of the region if not devisable by 3
				s += 'N'
			translated_record = SeqRecord(seq=s.translate(table=1), id=seq_record.id, description=seq_record.description+' frame='+str(i))
			datahandler_list2.append(translated_record)

			s = s.reverse_complement()   #checks the reverse compliment of the open reading frame of each gene in the fasta file
			while len(s)%3 != 0: # add N to the end of the region if not devisable by 3
				s += 'N'
			translated_record = SeqRecord(seq=s.translate(table=1), id=seq_record.id+"_reverse_compliment", description=seq_record.description+' frame='+str(i))
			datahandler_list2.append(translated_record)

	SeqIO.write(datahandler_list2, datahandler2, 'fasta')
	print '-'*20
	print '// Translated sequence to amino acids'
	print '-'*20

def OrfFinder(datahandler2, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m):
	for seq_record in SeqIO.parse(datahandler2, 'fasta', IUPAC.protein):
		genomic_region_up = int(seq_record.description.split('region:')[1].split('_')[0])
		genomic_region_down = int(seq_record.description.split('region:')[1].split('_')[1].split(' ')[0])

		MetStop(seq_record.seq, seq_record.id, genomic_region_up, genomic_region_down, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m)
	SeqIO.write(orfs, datahandler3, 'fasta')
	#print int(1/0)

	print '// Wrote %s protein sequences with a mimp IR motif in their promoter and an ORF >%i and <%i aa to %s' % (len(orfs), min_prot_len, max_prot_len, datahandler3)
	print '-'*20

def MetStop(sequence, ident, genomic_region_up, genomic_region_down, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m): #frame removec
	

	met_location = [i for i, a in enumerate(sequence) if a == 'M']

	x=0
	while x<len(met_location):
		if met_location[x] >= 0:# and stop_location>=0:
			stop_location = sequence.find('*', met_location[x])
			#print stop_location, "stop"
			if stop_location>0:
				prot = sequence[met_location[x]:stop_location+1] 	#+1 = add STOP
				dist_to_gene_fromstop = (stop_location*3)
			
				dist_to_gene = (met_location[x]*3)		 	#genomic distance

				if len(prot) > min_prot_len:
					startpos = genomic_region_up+dist_to_gene
					endpos = genomic_region_up+dist_to_gene_fromstop
					orf_record = SeqRecord(seq=prot.strip('*'), id=str(ident).replace('downstream', 'ds').replace('upstream', 'us') +"|"+str(startpos)+'-'+str(endpos)+'|d2m:'+str(dist_to_gene)+'|len:'+str(len(prot)-1), description='')
					if len(prot) < int(max_prot_len):
						orfs.append(orf_record)

		x+=1
	


#####important
def OrfWriter(datahandler3, signalpfile, min_prot_len, proteinoutfile, SignalPpath, SignalP_threshold,outdirectory):
	RunSignalP(datahandler3, signalpfile, 'euk', SignalPpath, SignalP_threshold,outdirectory)
	SPorfs = []
	files = os.listdir(outdirectory)
	orfs = SeqIO.parse(datahandler3, 'fasta', IUPAC.protein)

	with open(datahandler3, 'r') as file:
		data=file.readlines()
		file.close()
	test =True
	test2 =True
	test3 =True
	x=0
	check =[]

	for file in files:
		if '.summary_out' in file:
			with open(outdirectory+file, 'r') as sp:
				temp_file=sp.readlines()
				sp.close()
			for line in temp_file:
					
					if 'Name=' in line:
					
						sp_name = line.split('Name=')[1].split('\tSP=')[0]
						x=0
						for orf in data:
							if '>' in orf:

								orf_id=orf.split(" ")[0].split('>')[1]

								if sp_name in orf_id:

									sp_present = line.split('\tSP=\'')[1].split('\'')[0] # SP: YES OR NO
									orf_id += '|SP=' + sp_present
									orf_id += ';D='+line.split(' D=0')[1].split(' ')[0]  # SP D-value
									
									if 'Y' in sp_present:
										y=x+1 #add all seqence parts together 
										orf_seq=''

										while y<len(data) and "|"  not in data[y]:
											orf_seq+=data[y].split("\n")[0]
											y+=1

										cleavage_site_pos1 = line.split('Cleavage site between pos. ')[1].split(' and ')[0]
										cleavage_site_pos2 = line.split('Cleavage site between pos. ')[1].split(' and ')[1].split(': ')[0]

										if 'X' not in (orf_seq[:int(cleavage_site_pos2)-1]): #prevent 'NNN' translated to 'X' in signalpeptides to end up in the list
											temp_id=''
											for i in orf_id:
												if i==',':
													temp_id+=''
												else:
													temp_id+=i
											orf_id=temp_id
											signalpeptideseq = str(orf_seq[:int(cleavage_site_pos2)-1])
											seq = str(signalpeptideseq.lower()+orf_seq[int(cleavage_site_pos2)-1:].upper())
											orfie=">"+orf.replace('\n','')+' signalpeptideseq='+signalpeptideseq+'\n'+seq+"\n"
											if orfie not in SPorfs:
												SPorfs.append(orfie)
							x+=1
	sp.close()

	list(set(SPorfs))
	print '   Total # of sequences matching the criteria: %i' % len(SPorfs)
	with open(proteinoutfile, 'w') as file:
		file.writelines(SPorfs)		
		file.close()

def RunSignalP(datahandler3, signalpfile, organism, SignalPpath, SignalP_threshold,outdirectory):
	print '// Parsing data into batches...'
	with open(datahandler3, 'r') as file:
		data =file.readlines()
		file.close()
	if len(data) >20000:
		x=1
		storage=[]
		temp=[]
		y=0
		for seq in data:
			if len(temp)>=9998 and ">" in seq:
				storage.append(temp)
				y+=1
				temp=[]
			temp.append(str(seq))
			if len(storage) >10000:
				end() 
		storage.append(temp)
		y+=1
		z=1
		for s in storage:

			with open(outdirectory+"split_data_"+str(z)+".fasta","w") as file:
				file.writelines(s)
				file.close()
			z+=1
		print "Number of batches made: "+str(y)
		print '// Running SignalP 4.1 on batches...'
		while y>0:
			cline = SignalPpath+' -t %s -f summary -u %s %s > %s' % (organism, SignalP_threshold, outdirectory+"/split_data_"+str(y)+".fasta", signalpfile+str(y)+'.summary_out')
			os.system(cline)
			y-=1
		cline="rm -f "+ outdirectory+"split_data_*.fasta"  #Removes previous iteration data
		subprocess.check_output(['bash','-c', cline]) 
	else:
		cline = SignalPpath+' -t %s -f summary -u %s %s > %s' % (organism, SignalP_threshold, datahandler3, signalpfile+'.summary_out')
		os.system(cline)

#####


def ExtractOrfToFasta(proteinsfasta, uberinfile, puteff_dnaseqs, genome, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2):
	#EXTRACT GENOMIC SEQUENSE:
	open(puteff_dnaseqs, 'wb').close() 						#clear genomic DNA sequence fastafile
	dnaoutfile = open(puteff_dnaseqs, 'a')
	open(puteff_logfile, 'wb').close()

	logheader="genome\tputeff_supercontig\tgenomic_start_pos\tgenomic_end_pos\tdist2mimp\torientation\tprotlength\tD_value\tmimp_IR_seq\tmimp_IR_pos\tsignalpeptideseq\tproteinseq\tgenomicsequence\n"
	logheader2="genome\tnr_of_complete_mimps\tnr_of_incomplete_mimps\tnr_of_put_eff_MetStop\n"
	puteff_logfile_writer = open(puteff_logfile, 'a')
	puteff_logfile_writer.write(logheader) #for each genome, write a log with a header.

	combined_putefffile = open(combined_puteff_fasta, 'a')
	combined_puteff_logfile_writer = open(combined_puteff_logfile, 'a')
	combined_puteff_logfile2_writer = open(combined_puteff_logfile2, 'a')
	if filecounter < 1: # During running of the first genome file, the log header should be added at the top in the case of a combined log file.
		combined_puteff_logfile_writer.write(logheader)
		combined_puteff_logfile2_writer.write(logheader2)
	#proteinsfastafile = SeqIO.parse(proteinsfasta, 'fasta')
	with open(proteinsfasta, 'r') as file:
		proteinsfastafile=file.readlines()
		file.close()
	n=0
	#all_contigs = list(SeqIO.parse(uberinfile, "fasta", IUPAC.unambiguous_dna))
	with open(uberinfile, 'r') as file:
		all_contigs=file.readlines()
		file.close()

	from_previous_program=False
	for seq_record in proteinsfastafile:
		if '>' in seq_record:

			puteff_supercontig =seq_record.split('contig_')[1].split('_')[0].replace('\n','')
			genomic_start_pos = int(seq_record.split('|')[1].split('-')[0])
			genomic_end_pos = int(seq_record.split('|')[1].split('-')[1].split('|')[0])
			dist2mimp = seq_record.split('d2m:')[1].split('|')[0]
			protlength = seq_record.split('|len:')[1].split('|')[0]
			signalpeptideseq = seq_record.split('signalpeptideseq=')[1].replace('\n','')
			proteinseq = proteinsfastafile[n+1]
			x=0

	 		for sc in all_contigs:
	 			sc_id=''

	 			r_start=0
	 			r_end=0
	 			if '>' in sc:
	 				sc_id =sc.split('contig_')[1].replace('\n','')
	 				if '_' in sc_id:
	 					sc_id=sc_id.split('_')[0]
	 				if 'region' in sc:
	 					from_previous_program=True
	 					r_start =int(sc.split('region:')[1].split('_')[0])
	 					r_end =int(sc.split('region:')[1].split('_')[1].replace('\n',''))

				if from_previous_program==False and sc_id == puteff_supercontig:

					genomicsequence = all_contigs[x+1][genomic_start_pos-1:genomic_end_pos-1]#, len(sc.seq[genomic_start_pos-1:genomic_end_pos-1])
					print '   contig_'+str(puteff_supercontig)+'\tposition '+str(genomic_start_pos)+'-'+str(genomic_end_pos)+'\t'+signalpeptideseq
					putEff_fastaentry = ">"+str(n).zfill(4)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+str(genomicsequence)+"\n\n"
					dnaoutfile.write(putEff_fastaentry)
																													#orientation						#D_value, mimp_IR_seq, mimp_IR_pos
					puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, protlength, signalpeptideseq, proteinseq, genomicsequence]
					putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
					puteff_logfile_writer.write(putEff_logentry)

					combined_putefffile.write(putEff_fastaentry)
					combined_puteff_logfile_writer.write(putEff_logentry)

				elif from_previous_program==True and sc_id == puteff_supercontig and r_start<=genomic_start_pos and r_end>=genomic_end_pos:
					
					genomicsequence = all_contigs[x+1][(genomic_start_pos-r_start):(genomic_end_pos-r_start)]#, len(sc.seq[genomic_start_pos-1:genomic_end_pos-1])
					print '   contig_'+str(puteff_supercontig)+'\tposition '+str(genomic_start_pos)+'-'+str(genomic_end_pos)+'\t'+signalpeptideseq
					putEff_fastaentry = ">"+str(n).zfill(4)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+"\n"+str(genomicsequence)+"\n\n"
					dnaoutfile.write(putEff_fastaentry)
																													#orientation						#D_value, mimp_IR_seq, mimp_IR_pos
					puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, protlength, signalpeptideseq, proteinseq, genomicsequence]
					putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
					puteff_logfile_writer.write(putEff_logentry)

					#combined_putefffile will collect all the output from the mimpsearch; this means there will be many redundant put effectors.
					combined_putefffile.write(putEff_fastaentry)
					combined_puteff_logfile_writer.write(putEff_logentry)

				x+=1
		n+=1
	putEff_logentry2 = genome+'\t'+str(n)+'\n'
	combined_puteff_logfile2_writer.write(putEff_logentry2)
	dnaoutfile.close() #collects inside genome out folder all DNA sequences of the putative effectors
	puteff_logfile_writer.close() #writes a log for all puteff found in the current genome (inside genome out folder)
	combined_putefffile.close() #collects inside the out folder all DNA sequences of the putative effectors of all genomes that are being processed by the script.
	combined_puteff_logfile_writer.close() #writes a general log file with more details of the putative effectors identified
	combined_puteff_logfile2_writer.close()
	print '-'*20
	print "// Finished with genome of %s; wrote %i genomic DNA sequences of putEff ORFs to %s" % (genome, n/2, puteff_dnaseqs)

def MainDef(genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir,min_prot_len,max_prot_len,max_d2m,SignalPpath,SignalP_threshold):

	print directory, "  ", folder
	infilename, infileextension = os.path.splitext(genomefastafile)
	infile      	= directory+folder+'/'+genomefastafile
	outdirectory	= combined_puteff_dir+infilename.split("_downstream_")[0]+'_mimpfinder_out/'
	if not os.path.exists(outdirectory):
		os.makedirs(outdirectory)

	datahandler2	= outdirectory+infilename.split("_downstream_")[0]+'_1_gene_finder_translateddownstreamregions.fasta'
	datahandler3	= outdirectory+infilename.split("_downstream_")[0]+'_2_gene_finder_putativeORFs.fasta'
	signalpfile 	= outdirectory+infilename.split("_downstream_")[0]+'_3_gene_finder_SignalP'
	proteinoutfile 	= outdirectory+infilename.split("_downstream_")[0]+'_4_gene_finder_proteinseq_out.fasta'
	puteff_dnaseqs 	= outdirectory+infilename.split("_downstream_")[0]+'_5_gene_finder_puteff_genomicseq_out.fasta'
	puteff_logfile 	= outdirectory+infilename.split("_downstream_")[0]+'_6_gene_finder_puteff_logfile.txt'
	#motiefje    	= 'TT[TA]TTGC..CCCACTG..'
	#motiefje_rc 	= '..CAGTGGG..GCAA[TA]AA'
	#motiefje    	= 'TT[TAC]TTGC[ACG][CTA]C[CT][CT]ACTG..' ## (Mara Bergemann, 2008)
	#motiefje_rc 	= '..CAGT[GA][GA]G[GAT][TGC]GCAA[TAG]AA'

	orfs 			= []


	#nrofcompletemimps, nrofincompletemimps = MimpFinder(infile, sc_prefix, forward_mimps, reverse_mimps, datahandler, distance, mimpsequencesfile, infilename) #ir_dict[i] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[m.end():m.end()+distance]]

	Translator(infile, datahandler2)
	OrfFinder(datahandler2, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m)
	OrfWriter(datahandler3, signalpfile, min_prot_len, proteinoutfile, SignalPpath, SignalP_threshold,outdirectory)
	ExtractOrfToFasta(proteinoutfile, infile, puteff_dnaseqs, infilename, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2)

	#return mimpsequencesfile


class gene_finderApp():

	def __init__(self):
		self.verbose = False
	def start(self, directory_folder, output_dir, SignalPpath ='signalP', min_prot_len=10, max_d2m=2500, max_prot_len=134, SignalP_threshold=.45):

	#####################
		if directory_folder[-1]=="/":
			directory_folder=directory_folder[:-1]
		if output_dir[-1]=="/":
			output_dir=output_dir[:-1]

		folder = directory_folder.split('/')[-1]
		directory = directory_folder.split(folder)[0]

		file_extensions = (".fa", ".fasta", ".fas", "fna") # Specify the suffix of the genome files (.fa, .fasta, etc)
		combined_puteff_dir = output_dir+'/01.gene_finder/'+folder+'_MetStopOut/'

		if not os.path.exists(combined_puteff_dir):
			os.makedirs(combined_puteff_dir)

		combined_puteff_fasta	= combined_puteff_dir+'all_putative_genes.fasta'
		open(combined_puteff_fasta, 'w').close()
		combined_puteff_logfile = combined_puteff_dir+'all_putative_genes_log.txt'
		open(combined_puteff_logfile, 'w').close()
		combined_puteff_logfile2 = combined_puteff_dir+'all_putative_genes_log2.txt'
		open(combined_puteff_logfile2, 'w').close()
		######################
		filecounter=0
	 	for genomefastafile in os.listdir(directory_folder):
	 		if genomefastafile.endswith(file_extensions):
				print "\n// executing mimpfinder_MetStop script for "+genomefastafile

				MainDef(genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir, min_prot_len, max_prot_len, max_d2m, SignalPpath, SignalP_threshold)
				filecounter+=1
		
		print '-'*20
		print '// THE END'
		print "No more files with extension '%s' were found in directory '%s'" % (file_extensions, (directory+folder))
		print "Executed the script for %i files." % (filecounter)
		print '-'*20

