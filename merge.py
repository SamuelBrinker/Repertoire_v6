import os, sys
from Bio import SeqIO
class mergeApp():

	def __init__(self):
		self.verbose = False


	def start(self,infile, clusters,expanded ,blastdatabasedir,BLASTbindir, identity=80, coverage=80,E_VALUE_THRESH=.001,threads=1,alignments=3, force=True):

		wd=os.getcwd()
		print_not_merged=False

		#check relative path
		if '../' in infile or infile[0]!='/':
			dir = os.path.dirname(__file__)
			infile = os.path.join(dir, infile)
		if '../' in clusters or clusters[0]!='/':
			dir = os.path.dirname(__file__)
			clusters = os.path.join(dir, clusters)
		if '../' in BLASTbindir or BLASTbindir[0]!='/':
			dir = os.path.dirname(__file__)
			BLASTbindir = os.path.join(dir, BLASTbindir)
		if '../' in blastdatabasedir or blastdatabasedir[0]!='/':
			dir = os.path.dirname(__file__)
			blastdatabasedir = os.path.join(dir, blastdatabasedir)
		if '../' in expanded or expanded[0]!='/':
			dir = os.path.dirname(__file__)
			expanded = os.path.join(dir, expanded)

		if '/' in clusters and clusters[-1] !='/':
			database_store = blastdatabasedir+'/'+clusters.split('/')[-1].split('.fa')[0]
		else:
			database_store = blastdatabasedir+'/'+clusters.split('.fa')[0]	


		with open(expanded,'r') as file:
			expanded_file=file.readlines()
			file.close()

		#Sperates the clusters into spaces in an array, sorts by length of cluster. Each cluster given a name based on size with largest =1
		expanded_clusters=[]
		temp_line=''
		x=0
		y=0
		for line in expanded_file:
			if '-----' in line:
				if temp_line!='':
					expanded_clusters.append(temp_line)
					temp_line=''
					
				x+=1
			else:
				temp_line+=line.replace('\n','')+'\t'
		expanded_clusters.append(temp_line)
		expanded_clusters.sort(key = len, reverse=True)

		#attaches the sequnce to the sequence header
		x=0
		seq_cluster=list(SeqIO.parse(clusters, "fasta"))
		updated_expanded_clusters=[]
		temp_string=''
		while x <len(expanded_clusters):
			cluster=expanded_clusters[x].split('\t')
			for c in cluster:
				for seq in seq_cluster:
					if c.replace('_Reversed','')==str(seq.id):
						temp_string+='>'+str(seq.id)+'_cluster'+str(x+1)+'\n'+str(seq.seq)+'\n'
			updated_expanded_clusters.append(temp_string)
			temp_string=''
			x+=1


		cmnd = BLASTbindir+'/makeblastdb -dbtype nucl -in '+clusters+' -out '+database_store
		print "---BLASTDB---\n", cmnd, os.system(cmnd), "---\n" #uncomment if you want to rebuild a blastdb #untag # if you want to build the BLAST database
		blastoutfilename = wd+'/'+infile.split('/')[-1].split('.fa')[0]+'.vs.'+clusters.split('/')[-1].split('.fa')[0]+'.blastout'
		cmnd = BLASTbindir+"/blastn -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' -query "+infile+' -db '+database_store+' -out '+blastoutfilename+' -evalue '+str(E_VALUE_THRESH)+" -num_threads "+str(threads)+" -num_alignments "+str(alignments)
												#0		1		2		3	4		5	6		7	8	9
		os.system(cmnd)
	
		with open(blastoutfilename, 'r') as file:
			lines=file.readlines()
			file.close()
		with open(infile, 'r') as file:
			raw_infile=file.readlines()
			file.close()

		not_merged=[]
		for line in raw_infile:
			if '>' in line:
				not_merged.append(line[1:].replace('\n',''))
		
		test=True
		test2=True
		merged=[]
		output_file=[]
		for line in lines:
			tabs=line.split('\t')
			q_score =abs((float(tabs[5]) - float(tabs[4]))/float(tabs[8])) #covereage calculated by q_score/s_score ratio
			s_score =abs((float(tabs[7]) - float(tabs[6]))/float(tabs[9])) # the score value = the length of the matached section divided by the total length of the effector
			
			if float(tabs[2]) >= identity and float(q_score)>=coverage and float(s_score)>=coverage:
				try:
					for line in not_merged:
						if " " in line:
							edited_line = line.split(" ")[0]
							#if test==True:
							#	print [edited_line], [">"+tabs[0]]
							#	if edited_line==">"+tabs[0]:
							#		print "!!!"
							#	test=False
						else:
							edited_line=line
						if edited_line == tabs[0]:
							not_merged.remove(line)

							#if test2==True:
							#	print "!"					
							#	test2=False

					#not_merged.remove(tabs[0])
				except:
					print tabs[0]
					pass

				if tabs[1] in merged:
					merged[merged.index(tabs[1])+1]+=tabs[0]+'\t'
					if tabs[0].replace('\n','')=='':
						print tabs[0]
				else:
					merged.append(tabs[1])
					merged.append(tabs[0].replace('\n','')+'\t')
					if tabs[0].replace('\n','')=='':
						print tabs[0]
				output_file.append(tabs[0]+'\t'+tabs[1]+'\t'+tabs[2]+'\t'+tabs[3]+'\t'+str(q_score)+'\t'+str(s_score)+'\n')

		#print len(not_merged), len(merged)
		#print [merged[2:4]]
		

		edited_merge=[]
		x=0

		seq_infile = list(SeqIO.parse(infile, "fasta"))
		seq_cluster=list(SeqIO.parse(clusters, "fasta"))

		############## new stuff

		for raw_line in output_file:
			line = raw_line.split('\t')
			for seq in seq_cluster:
				if str(seq.id) == line[1]:
					cluster_seq=seq
			for seq in seq_infile:
				if str(seq.id) == line[0]:
					to_merge_seq=seq
			x=0
			while x< len(updated_expanded_clusters):
				cluster=updated_expanded_clusters[x].split('\n')
				for c in cluster:
					if '_cluster' in c:
						if str(cluster_seq.id) == c.split('_cluster')[0].split('>')[1]:
							output_file[output_file.index(raw_line)]=line[0]+'\t'+str(c).replace('>','')+'\t'+line[2]+'\t'+str(line[3])+'\t'+str(line[4])+str(line[5])
							if str(to_merge_seq.id) not in updated_expanded_clusters[x]:
								updated_expanded_clusters[x]+='>'+str(to_merge_seq.id)+'_cluster'+c.split('_cluster')[1].split('\n')[0]+'\n'+str(to_merge_seq.seq)+'\n'
								#print 'test'
							break

				x+=1 

		##############



		#new_header_list=[]
		z=0
		f=0
		x=0
		while x<len(merged):
				
			y=0
			temp=merged[x]+'\t'
			for seq in seq_cluster:
				if str(seq.id) == merged[x]:
					old_id=seq
					cluster_id=seq
			to_merge=merged[x+1].split('\t')
			for line in to_merge:
				z+=1
				debug=True	
				for seq in seq_infile:
					if str(seq.id) == line:
						y+=1
						f+=1
						if len(seq.seq)>len(cluster_id.seq):
							cluster_id=seq

			###### This code would print new header names
			'''
			if str(old_id.id)!=str(cluster_id.id):
				for new_name in expanded_clusters:
					if new_name.split('\n')[0]==str(cluster_id.id):
						new_header_list.append('>'+str(new_name.split('\n')[0])+' (was '+str(old_id.id)+')\n'+str(cluster_id.seq)+'\n')
			else:
				for new_name in expanded_clusters:
					if new_name.split('\n')[0]==str(cluster_id.id):
						new_header_list.append('>'+str(new_name.split('\n')[0])+'\n'+str(cluster_id.seq)+'\n')
			'''


			edited_merge.append(temp+'\t'+str(y)+'\n')
			x+=2

		'''
		for seq in seq_cluster:
			if seq.id not in merged:
				new_header_list.append('>'+str(seq.id)+'\n'+str(seq.seq)+'\n')
		'''

		if len(not_merged)>0:
			edited_merge.append("Number of sequences not merged\t"+str(len(not_merged)))
			print_not_merged=True
			to_print_not_merged=[]
			x=0
			check=False
			while x<len(raw_infile):
				if check==True and '>' in raw_infile[x]:
					check=False
				elif check==True:
					to_print_not_merged.append(raw_infile[x])
				for line in not_merged:
					if line in raw_infile[x]:
						not_merged.remove(line)
						to_print_not_merged.append(raw_infile[x])
						check=True
				x+=1

		blast_results="blast_results_run1.txt"
		summery="summery_run1.txt"
		not_merged_output="not_merged_run1.txt"
		new_headers='new_cluster_names_run1.txt'
		sorted_sequences='sorted_sequences_run1.fasta'
		if force != True:
			x=2
			while os.path.isfile(summery) ==True and os.path.isfile(sorted_sequences) ==True and  os.path.isfile(blast_results) ==True and  os.path.isfile(not_merged_output) ==True and  os.path.isfile(new_headers) ==True:
				summery="summery_run"+str(x)+'.txt'
				blast_results="blast_results_run"+str(x)+'.txt'
				not_merged_output="not_merged_run"+str(x)+".txt"
				new_headers='new_cluster_names_run'+str(x)+'.txt'
				sorted_sequences=='sorted_sequences_run'+str(x)+'.fasta'
				x+=1

		with open(blast_results,'w') as file:
			file.writelines(edited_merge)
			file.close()
		#with open(new_headers,'w') as file:
		#	file.writelines(new_header_list)
		#	file.close()
		with open(summery,'w') as file:
			file.writelines(output_file)
			file.close()
		with open(sorted_sequences, 'w') as file:
			for line in updated_expanded_clusters:
				file.write(line+'\n--------\n')
			file.close()
		if print_not_merged==True:
			with open(not_merged_output,'w') as file:
				file.writelines(to_print_not_merged)
				file.close()







