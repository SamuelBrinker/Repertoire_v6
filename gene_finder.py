
import sys, os, re, glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

def Translator(infile, datahandler2):
	datahandler_list2 = []
	datahandler_nucl_list2 = []
	for seq_record in SeqIO.parse(infile, 'fasta'): ###PMH: removed iupac designation
		for i in range(3): 				#(i=0 - i=1 - i=2)
			s = seq_record.seq[i:]   #starts the sequence at the open reading frame "i" (0, 1, or 2)
			while len(s)%3 != 0: # add N to the end of the region if not devisable by 3
				s += 'N'
			translated_record = SeqRecord(seq=s.translate(table=1), id=seq_record.id, description=seq_record.description+' frame='+str(i))
			datahandler_list2.append(translated_record)

			s = seq_record.seq[:-i].reverse_complement()   #reverse complements the sequence
			while len(s)%3 != 0: # add N to the end of the region if not devisable by 3
				s += 'N'
			translated_record = SeqRecord(seq=s.translate(table=1), id=seq_record.id+"_reverse_complement", description="frame="+str(i))
			datahandler_list2.append(translated_record)

	SeqIO.write(datahandler_list2, datahandler2, 'fasta')

	print('-'*20)
	print('// Translated sequence to amino acids')
	print('-'*20)

def OrfFinder(datahandler2, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m): #orfs is a blank, it is not duplicating up/down regions
	u=0
	unique_orfs=[]
	for seq_record in SeqIO.parse(datahandler2, 'fasta', IUPAC.protein):
		genomic_region_start = int(seq_record.description.split('region:')[1].split('_')[0])
		genomic_region_stop = int(seq_record.description.split('region:')[1].split('_')[1].split(' ')[0])
		frame = int(seq_record.description.split('frame=')[1].split('\n')[0])
		if str(seq_record.description.split(' ')[0].split('_')[-1]) == str("complement"):
			rc = True
		else:
			rc = False

		
		MetStop(seq_record.seq, seq_record.id, genomic_region_start, genomic_region_stop, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m, frame, rc,u,seq_record)
		

		u+=1

	SeqIO.write(orfs, datahandler3, 'fasta')
	'''
	for o in orfs:
		if '_us' in o.id:
			check_name=o.id.split('_us')[0]
			check_location=o.id.split('|')[1]+'|'
		elif '_upstream' in o.id:
			check_name=o.id.split('_upstream')[0]
			check_location=o.id.split('|')[1]+'|'
		else:
			check_name=o.id.split('_ds')[0]
			check_location=o.id.split('|')[1]+'|'

		if 'reverse_complement' in o.id:
			unique_orf_record=o
			unique_orf_record.id=check_name+'_r_'+check_location
			unique_orfs.append(unique_orf_record)
		else:
			unique_orf_record=o
			unique_orf_record.id=check_name+check_location
			unique_orfs.append(unique_orf_record)
	SeqIO.write(unique_orfs, datahandler3.split('.fasta')[0]+'_unique_names.fasta', 'fasta')
	'''
	#print [all_orfs]
	#with open (datahandler3+'_test.fasta','w') as file:
	#	file.writelines(str(all_orfs))
	#	file.close()
	#print int(1/0)
	#print u
	print('// Wrote %s protein sequences with a mimp IR motif in their promoter and an ORF >%i and <%i aa to %s' % (len(orfs), min_prot_len, max_prot_len, datahandler3))
	print('-'*20)

def MetStop(sequence, ident, genomic_region_start, genomic_region_stop, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m, frame, rc,u, seq_record): #frame removec
	
	met_location = [i for i, a in enumerate(sequence) if a == 'M'] 
	#print [met_location]
	region_length = int(genomic_region_stop)-int(genomic_region_start) #for calculating start/stop locations for orfs in reverse complement of up/down region
	#print len(orfs)
	add_to_orfs=True
	x=0
	while x<len(met_location):
		
		if met_location[x] >= 0:# and stop_location>=0:
			stop_location = sequence.find('*', met_location[x])
			#print stop_location, "stop"
			if stop_location>0:
				prot = sequence[met_location[x]:stop_location+1] 	#+1 = add STOP
				if rc == False:	
					if frame == 0:
						dist_to_start = (met_location[x]*3)		 	
						dist_to_stop = (stop_location*3)
					if frame == 1:
						dist_to_start = ((met_location[x]*3)+1)	
						dist_to_stop = ((stop_location*3)+1)
					if frame == 2:
						dist_to_start = ((met_location[x]*3)+2)		 	
						dist_to_stop = ((stop_location*3)+2)
					if len(prot) > min_prot_len:
						startpos = genomic_region_start+dist_to_start
						endpos = genomic_region_start+dist_to_stop
						orf_record = SeqRecord(seq=prot.strip('*'), id=str(ident).replace('downstream', 'ds').replace('upstream', 'us') +"|"+str(startpos)+'-'+str(endpos)+'|d2m:'+str(dist_to_start)+'|len:'+str(len(prot)-1)+"|rc="+str(rc), description='')
						if len(prot) < int(max_prot_len):
							if '_us' in orf_record.id:
								check_name=orf_record.id.split('_us')[0]
								check_location='|'+str(startpos)+'-'+str(endpos)+'|'
							elif '_upstream' in orf_record.id:
								check_name=orf_record.id.split('_upstream')[0]
								check_location='|'+str(startpos)+'-'+str(endpos)+'|'
							else:
								check_name=orf_record.id.split('_ds')[0]
								check_location='|'+str(startpos)+'-'+str(endpos)+'|'
						
							t=0
							
							#while t<len(orfs):
								#if '1191186-1191414' in str(orfs[t].id) and '1191186-1191414' in check_location:
								#	print [orfs[t].id]
								#	print [check_name],[check_location]
								#	print check_name
								#	print ident.split('_downstream')[0], '\n\n'

								#if check_name in orfs[t].id and check_location in orfs[t].id:
								#	if '1191186-1191414' in str(orfs[t].id):
								#		print "huh"
								#	t=len(orfs)+3
									
								#t+=1
							for o in orfs:
								if check_name in o.id and check_location in o.id:
									#if '1191186-1191414' in str(o.id):
									#	print "huh"
									add_to_orfs=False

							if add_to_orfs:
							#if t!=len(orfs)+3:
								orf_record.id='id'+str(u)+'_'+str(x)+'_'+orf_record.id
								orfs.append(orf_record)
								

								
								#if '1191186-1191414' in str(orf_record):
								#	print 'test  ', u, seq_record.id, seq_record.description, '\n', orf_record.id
								#	print orfs[-1].id
									#print 'fail'
								#	print check_name
								#	print check_location
							
				else:
					if frame == 0:
						dist_to_start = (region_length-met_location[x]*3)
						dist_to_stop = (region_length-stop_location*3)		 	
					if frame == 1:
						dist_to_start = ((region_length-met_location[x]*3)-1)
						dist_to_stop = ((region_length-stop_location*3)-1)		 
					if frame == 2:
						dist_to_start = ((region_length-met_location[x]*3)-2)
						dist_to_stop = ((region_length-stop_location*3)-2)
					if len(prot) > min_prot_len:
						startpos = genomic_region_start+dist_to_start
						endpos = genomic_region_start+dist_to_stop
						orf_record = SeqRecord(seq=prot.strip('*'), id=str(ident).replace('downstream', 'ds').replace('upstream', 'us') +"|"+str(startpos)+'-'+str(endpos)+'|d2m:'+str(dist_to_start)+'|len:'+str(len(prot)-1)+"|rc="+str(rc), description='')
						if len(prot) < int(max_prot_len):
							if '_us' in orf_record.id:
								check_name=orf_record.id.split('_us')[0]
								check_location='|'+str(startpos)+'-'+str(endpos)+'|'
							elif '_upstream' in orf_record.id:
								check_name=orf_record.id.split('_upstream')[0]
								check_location='|'+str(startpos)+'-'+str(endpos)+'|'
							else:
								check_name=orf_record.id.split('_ds')[0]
								check_location='|'+str(startpos)+'-'+str(endpos)+'|'
					
							
							for o in orfs:
								if check_name in o.id and check_location in o.id:
									#if '1191186-1191414' in str(o.id):
									#	print "huh"
									add_to_orfs=False
							
							if add_to_orfs:
							#if t!=len(orfs)+3:
								orf_record.id='temp_id'+str(u)+'_'+orf_record.id
								orfs.append(orf_record)
								#if 'reverse_complement' in orf_record.id:
								#	unique_orf_record=orf_record
								#	unique_orf_record.id=check_name+'_r_'+check_location
								#	unique_orfs.append(unique_orf_record)
								#else:
								#	unique_orf_record=orf_record
								#	unique_orf_record.id=check_name+check_location
								#	unique_orfs.append(unique_orf_record)

								
								#if '1191186-1191414' in str(orf_record):
								#	print 'test  ', u, seq_record.id, seq_record.description, '\n', orf_record.id
									#print 'fail'
								#	print check_name
								#	print check_location
							

		x+=1
	#return orfs


#####important
def OrfWriter(datahandler3, signalpfile, min_prot_len, proteinoutfile, SignalPpath, SignalP_threshold,outdirectory):
	RunSignalP(datahandler3, signalpfile, 'euk', SignalPpath, SignalP_threshold,outdirectory)
	SPorfs = []
	files = os.listdir(outdirectory)
	orfs = SeqIO.parse(datahandler3, 'fasta', IUPAC.protein)

	with open( datahandler3, 'r') as file:
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

								orf_id=orf.split('>')[1].split(' ')[0]
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
											orfie=orf.replace('\n','')+' signalpeptideseq='+signalpeptideseq+'\n'+seq+"\n\n"
											if orfie not in SPorfs:
												SPorfs.append(orfie)
							x+=1
	sp.close()

	list(set(SPorfs))
	print('   Total # of sequences matching the criteria: %i' % len(SPorfs))
	#print [SPorfs[0:3]]
	with open(proteinoutfile, 'w') as file:
		file.writelines(SPorfs)		
		file.close()

def RunSignalP(datahandler3, signalpfile, organism, SignalPpath, SignalP_threshold,outdirectory):
	print('// Parsing data into batches...')
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
		print("Number of batches made: "+str(y))
		print('// Running SignalP 4.1 on batches...')
		while y>0:
			cline = SignalPpath+' -t %s -f summary -u %s %s > %s' % (organism, SignalP_threshold, outdirectory+"/split_data_"+str(y)+".fasta", signalpfile+str(y)+'.summary_out')
			os.system(cline)
			y-=1
#		cline="rm -f "+ outdirectory+"split_data_*.fasta"  #Removes previous iteration data
		subprocess.check_output(['bash','-c', cline]) 
	else:
		print('// Running SignalP 4.1 on batches...')
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
	#test=True
	from_previous_program=False
	for seq_record in proteinsfastafile:
		if '>' in seq_record:
			#print(seq_record)
			puteff_supercontig=''
			if '_ds_' in seq_record:
				temp_record_id=seq_record.split('_ds_')[0]
				temp_record_id=temp_record_id.split('_')[2:]
				for x in temp_record_id:
					puteff_supercontig+=x+'_'
			else:
				temp_record_id=seq_record.split('_us_')[0]
				temp_record_id=temp_record_id.split('_')[2:]
				for x in temp_record_id:
					puteff_supercontig+=x+'_'
			puteff_supercontig=puteff_supercontig[:-1]

			genomic_start_pos = int(seq_record.split('|')[1].split('-')[0])
			genomic_end_pos = int(seq_record.split('|')[1].split('-')[1].split('|')[0])
			dist2mimp = seq_record.split('d2m:')[1].split('|')[0]
			protlength = seq_record.split('|len:')[1].split('|')[0]
			signalpeptideseq = seq_record.split('signalpeptideseq=')[1].replace('\n','')
			rc = seq_record.split('=')[1].split(' ')[0] #ORF from reverse complement of upstream/downstream or not
			proteinseq = proteinsfastafile[n+1]
			x=0

			for sc in all_contigs:
				#print('*screams internally*')
				sc_id=''

				r_start=0
				r_end=0
				if '>' in sc:
					if 'downstream' in sc:
						sc_id=sc.split('_downstream_')[0].split('>')[1]
						#if test==True:
						#	print(sc_id, puteff_supercontig)
						#	test=False
					else:
						sc_id=sc.split('_upstream')[0].split('>')[1]
					if 'region' in sc:
						from_previous_program=True
						r_start =int(sc.split('region:')[1].split('_')[0])
						r_end =int(sc.split('region:')[1].split('_')[1].replace('\n',''))

				if from_previous_program==False and sc_id == puteff_supercontig:
					
					genomicsequence = all_contigs[x+1][genomic_start_pos-1:genomic_end_pos-1]#, len(sc.seq[genomic_start_pos-1:genomic_end_pos-1])
					print('   '+str(puteff_supercontig)+'\tposition '+str(genomic_start_pos)+'-'+str(genomic_end_pos)+'\t'+signalpeptideseq)
					putEff_fastaentry = ">"+str(n).zfill(4)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+"\n"+str(genomicsequence)+"\n"
					dnaoutfile.write(putEff_fastaentry)
																													#orientation						#D_value, mimp_IR_seq, mimp_IR_pos
					puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, protlength, signalpeptideseq, proteinseq, genomicsequence]
					putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
					puteff_logfile_writer.write(putEff_logentry)

					combined_putefffile.write(putEff_fastaentry)
					combined_puteff_logfile_writer.write(putEff_logentry)

				elif from_previous_program==True and sc_id == puteff_supercontig:
					
					if rc == "True": #if ORF is from reverse complement (rc) of up/downstream region
						
						genomicsequence = Seq(all_contigs[x+1][(genomic_end_pos-r_start):(genomic_start_pos-r_start)], generic_dna)# for sequences from rc: end positions < start positions
						if genomicsequence!='':
							genomicsequence = genomicsequence.reverse_complement()
							print('   '+str(puteff_supercontig)+'\tposition '+str(genomic_start_pos)+'-'+str(genomic_end_pos)+'\t'+signalpeptideseq)
							putEff_fastaentry = ">"+str(n).zfill(4)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+"\n"+str(genomicsequence)+"\n"
							dnaoutfile.write(putEff_fastaentry)
																															#orientation						#D_value, mimp_IR_seq, mimp_IR_pos
							puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, protlength, signalpeptideseq, proteinseq, genomicsequence]
							putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
							puteff_logfile_writer.write(putEff_logentry)

							#combined_putefffile will collect all the output from the mimpsearch; this means there will be many redundant put effectors.
							combined_putefffile.write(putEff_fastaentry)
							combined_puteff_logfile_writer.write(putEff_logentry)
							break

					elif rc == "False" and r_start<=genomic_start_pos and r_end>=genomic_end_pos:
						genomicsequence = all_contigs[x+1][(genomic_start_pos-r_start):(genomic_end_pos-r_start)]#, len(sc.seq[genomic_start_pos-1:genomic_end_pos-1])
						if genomicsequence!='':
							print('   '+str(puteff_supercontig)+'\tposition '+str(genomic_start_pos)+'-'+str(genomic_end_pos)+'\t'+signalpeptideseq)
							putEff_fastaentry = ">"+str(n).zfill(4)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+"\n"+str(genomicsequence)+"\n"
							dnaoutfile.write(putEff_fastaentry)
																															#orientation						#D_value, mimp_IR_seq, mimp_IR_pos
							puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, protlength, signalpeptideseq, proteinseq, genomicsequence]
							putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
							puteff_logfile_writer.write(putEff_logentry)

							#combined_putefffile will collect all the output from the mimpsearch; this means there will be many redundant put effectors.
							combined_putefffile.write(putEff_fastaentry)
							combined_puteff_logfile_writer.write(putEff_logentry)
							break
					#else:
					#	print('error')
					#	print(seq_record)
					#	print(puteff_supercontig)
					#	print(sc, sc_id)
					#	print(r_start,r_end,genomic_start_pos,genomic_end_pos)
				x+=1
		n+=1
	putEff_logentry2 = genome+'\t'+str(n)+'\n'
	combined_puteff_logfile2_writer.write(putEff_logentry2)
	dnaoutfile.close() #collects inside genome out folder all DNA sequences of the putative effectors
	puteff_logfile_writer.close() #writes a log for all puteff found in the current genome (inside genome out folder)
	combined_putefffile.close() #collects inside the out folder all DNA sequences of the putative effectors of all genomes that are being processed by the script.
	combined_puteff_logfile_writer.close() #writes a general log file with more details of the putative effectors identified
	combined_puteff_logfile2_writer.close()
	print('-'*20)
	print("// Finished with genome of %s; wrote %i genomic DNA sequences of putEff ORFs to %s" % (genome, n/2, puteff_dnaseqs))

def MainDef(genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir,min_prot_len,max_prot_len,max_d2m,SignalPpath,SignalP_threshold,force):

	print(directory, "  ", folder)
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

	if force!=True and (os.path.isfile(datahandler2)==True or os.path.isfile(datahandler3)==True or os.path.isfile(signalpfile)==True or os.path.isfile(proteinoutfile)==True or os.path.isfile(puteff_dnaseqs)==True or os.path.isfile(puteff_logfile)==True):
		pass
	else:
		print(("translator infile is: "+infile))
		Translator(infile, datahandler2)
		OrfFinder(datahandler2, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m)
		OrfWriter(datahandler3, signalpfile, min_prot_len, proteinoutfile, SignalPpath, SignalP_threshold,outdirectory)
		ExtractOrfToFasta(proteinoutfile, infile, puteff_dnaseqs, infilename, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2)

	#return mimpsequencesfile


class gene_finderApp():

	def __init__(self):
		self.verbose = False
	def start(self, force, directory_folder, output_dir, SignalPpath ='signalP', min_prot_len=10, max_d2m=2500, max_prot_len=134, SignalP_threshold=.45, working=''):

	#####################
		if working!='':
			if working[-1]!='/':
				working+='/'
			if '../' in working:
				dir = os.path.dirname(__file__)
				working = os.path.join(dir, working)
			os.chdir(working)
		else:
			working=os.path.dirname(__file__)

		if '../' in directory_folder or directory_folder[0]!='/':
			dir = os.path.dirname(working)
			directory_folder = os.path.join(dir, directory_folder)
			directory_folder = os.path.abspath(os.path.realpath(directory_folder))

		if '../' in output_dir or output_dir[0]!='/':
			dir = os.path.dirname(working)
			output_dir = os.path.join(dir, output_dir)
			output_dir = os.path.abspath(os.path.realpath(output_dir))

		if '../' in SignalPpath or SignalPpath[0]!='/':
			dir = os.path.dirname(working)
			SignalPpath = os.path.join(dir, SignalPpath)
			SignalPpath = os.path.abspath(os.path.realpath(SignalPpath))

		if directory_folder[-1]=="/":
			directory_folder=directory_folder[:-1]
		if output_dir[-1]=="/":
			output_dir=output_dir[:-1]

		folder = directory_folder.split('/')[-1]
		directory = directory_folder.split(folder)[0]

		file_extensions = (".fa", ".fasta", ".fas", "fna") # Specify the suffix of the genome files (.fa, .fasta, etc)
		combined_puteff_dir = output_dir+'/'
#		combined_puteff_dir = output_dir+'/01.gene_finder/'+folder+'_MetStopOut/' ##old version

		if not os.path.exists(combined_puteff_dir):
			os.makedirs(combined_puteff_dir)

		combined_puteff_fasta	= combined_puteff_dir+'all_putative_genes.fasta'
		#open(combined_puteff_fasta, 'w').close()
		combined_puteff_logfile = combined_puteff_dir+'all_putative_genes_log.txt'
		#open(combined_puteff_logfile, 'w').close()
		combined_puteff_logfile2 = combined_puteff_dir+'all_putative_genes_log2.txt'
		#open(combined_puteff_logfile2, 'w').close()

		if force!=True:
			x=2
			while os.path.isfile(combined_puteff_fasta)==True and os.path.isfile(combined_puteff_logfile)==True and os.path.isfile(combined_puteff_logfile2)==True:
				combined_puteff_fasta	= combined_puteff_dir+'all_putative_genes_run_'+str(x)+ '.fasta'
				combined_puteff_logfile = combined_puteff_dir+'all_putative_genes_log_run_'+str(x)+ '.txt'
				combined_puteff_logfile2 = combined_puteff_dir+'all_putative_genes_log2_run_'+str(x)+ '.txt'
				x+=1
		else:
			open(combined_puteff_fasta, 'w').close()
			open(combined_puteff_logfile, 'w').close()
			open(combined_puteff_logfile2, 'w').close()
		######################
		filecounter=0
		for genomefastafile in os.listdir(directory_folder):
			#print("Examining "+ str(genomefastafile))
			if genomefastafile.endswith(file_extensions):
				print("\n// executing mimpfinder_MetStop script for "+genomefastafile)

				MainDef(genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir, min_prot_len, max_prot_len, max_d2m, SignalPpath, SignalP_threshold,force)
				filecounter+=1
		
		print('-'*20)
		print('// THE END')
		print("No more files with extension '%s' were found in directory '%s'" % (file_extensions, (directory+folder)))
		print("Executed the script for %i files." % (filecounter))
		print('-'*20)

