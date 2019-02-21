#this script takes in the output of mimp-identification script "all_putative_effectors.fasta"
#It will take each put effector in this fastafile and blast it against itself. Then it takes these groups 
#together and clusters them with single_linkage. Each of the clusters is assessed for the longest
#sequence, and these longest sequences are written to the outfile (_clustered.fasta)


#python 02.cluster_putefflists.py all_putative_effectors_MetStop, blastdatabasedir, leave_put_eff_identifiers_during_clustering, BLASTbindir
				#sysarg1 = ../input_files_for_testing/all_putative_effectors_concatenated.fasta
				#sysarg2 = /Users/pmh/Documents/Sam_Brinker/effectors/input_files_for_testing/blastdbdir
				#sysarg3 = TRUE
				#sysarg4 = /usr/local/bin

from Bio import SeqIO
import os, sys
def remove_n(infile, n_allowed):
	with open(infile, 'r') as file:
		raw=file.readlines()
		file.close()

	edited_data=[]
	sequence=False
	header=''

	fixed_fasta=[]
	temp_line=''
	for line in raw:

		if '>' in line:
			if temp_line!='':
				fixed_fasta.append(temp_line+'\n')
			fixed_fasta.append(line)
			temp_line=''
		else:
			temp_line+=line.replace('\n','')

	if n_allowed >1:
		n_allowed=float(n_allowed/100)

	for line in fixed_fasta:
		if '\n'!=line and ''!=line:
			if sequence ==True:
				if float(line[:-2].upper().count('N'))/len(line)<=n_allowed:
					edited_data.append(header)
					edited_data.append(line)
						#if "08_fo_con_54008_MIMP_TIR_Element_31" in header:
						#	print header, line
						#if '08_fo_con_54008_MIMP_TIR_Element_66' in header:
						#	print 'tes', header, [line]
					sequence=False
				#elif '\n'!=line or ''!=line:
				#	edited_data.append(header)
				#	edited_data.append(line)
					#if "08_fo_con_54008_MIMP_TIR_Element_31" in header:
					#	print header, line
					#if '08_fo_con_54008_MIMP_TIR_Element_66' in header:
					#		print 'test1'
					#sequence=False
			if '>' in line:
				header=line 
				#if '08_fo_con_54008_MIMP_TIR_Element_66' in header:
				#	print 'tes1'
				sequence=True
	infile_filename = infile.split('/')[-1]
	infiledir = infile.split(infile_filename)[0]
	with open(infiledir+'edited_'+infile_filename,'w') as file:
		file.writelines(edited_data)
		file.close()



def cluster_homologous_effectors(threads, infile, E_VALUE_THRESH, PERC_IDENTITY_THRESH, LENGTH_THRESH, blastdatabasedir, BLASTbindir, infiledir, low_identity, low_length_thresh, lihc, hilc,all_clusters, alignments):
	
	if '/' in infile and infile[-1] !='/':
		database_store = blastdatabasedir+'/'+infile.split('/')[-1].split('.fa')[0]
	else:
		database_store = blastdatabasedir+'/'+infile.split('.fa')[0]


	with open(infile,'r') as file:
		t=file.readlines()
		file.close()
	temp=[]
	for x in t:
		if x!='\n' or x!="":
			temp+=x.replace(' sign','_signal')
	with open(infile,'w') as g:
		g.writelines(temp)
		g.close()

	cmnd = BLASTbindir+'/makeblastdb -dbtype nucl -in '+infile+' -out '+database_store
	print "---BLASTDB---\n", cmnd, os.system(cmnd), "---\n" #uncomment if you want to rebuild a blastdb #untag # if you want to build the BLAST database
	blastoutfilename = infiledir+infile.split('/')[-1].split('.fa')[0]+'.vs.'+infile.split('/')[-1].split('.fa')[0]+'.blastout'
	cmnd = BLASTbindir+"/blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -dust no -query "+infile+' -db '+database_store+' -out '+blastoutfilename+' -evalue '+str(E_VALUE_THRESH)+" -num_threads "+str(threads)+" -num_alignments "+str(alignments)
	#										0	  1	    2	   3	  4 		5	    6	   7		8	  9 	10	 11	   12	 13

	#print hilc, lihc

	effector2homologs = {}
	high_i_high_c=[] #high identity, high coverage
	high_i_low_c=[]
	low_i_high_c=[]
	singletons=[] #things that do not fit into above
	precluster_catagories = [high_i_high_c, high_i_low_c, low_i_high_c, singletons]
	cluster_catagories=[]
	i=0

	if os.system(cmnd) == 0:
		lines = open(blastoutfilename).readlines()
			#'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen    slen'
		while i< 4:
			if i==1 and hilc!=True:
				cluster_catagories.append([[]])
				i+=1
			if i==2 and lihc!=True:
				cluster_catagories.append([[]])
				i+=1

			if i==0:
				print ("Finding high identity high coverage clusters")
			elif i==1:
				print ("Finding high identity low coverage clusters")
			elif i==2:
				print ("Finding low identity high coverage clusters")
			else:
				print ("Finding singletons")
			y = 0.0

			check_25=True
			check_50=True
			check_75=True
			check_100=True

			for line in lines:
				if (y/len(lines))>=.25 and check_25==True:
					print "25%..."
					check_25=False
				elif (y/len(lines))>=.50 and check_50==True:
					print "50%..."
					check_50=False
				elif (y/len(lines))>=.75 and check_75==True:
					print "75%..."
					check_75=False
				elif (y/len(lines))>=1 and check_100==True:
					print "100%... clustering"
					check_100=False

				q_new=True
				s_new=True

				tabs = line.strip().split('\t')				
				q_score =abs((float(tabs[7]) - float(tabs[6]))/float(tabs[12])) #covereage calculated by q_score/s_score ratio
				s_score =abs((float(tabs[9]) - float(tabs[8]))/float(tabs[13])) # the score value = the length of the matached section divided by the total length of the effector
				
				c_score=1

				if i==0: #high_i_high_c, effectors have to be above the length identity threshold and above the coverage threshold
					if float(tabs[2]) >= PERC_IDENTITY_THRESH and (c_score) >= (LENGTH_THRESH) and tabs[0]!=tabs[1] and q_score>LENGTH_THRESH and s_score>LENGTH_THRESH: #ignoes self matches
						x=0
						q_new=True
						s_new=True
						while x<len(high_i_high_c):
							if q_new!=False and tabs[0] == high_i_high_c[x][0] : #checks to see if effector is already in the set to avoid duplicates 																 
								high_i_high_c[x][1].append(tabs[1]) #add all BLAST-associated hits to this entry
								q_new=False
							if s_new != False and tabs[1] == high_i_high_c[x][0]: 
								high_i_high_c[x][1].append(tabs[0]) 
								s_new=False
							if q_new== False and s_new==False:
								x=len(high_i_high_c)
							x+=1
						if q_new==True:
							high_i_high_c.append([tabs[0],[tabs[1]]])
							q_new=False
						if s_new==True:	
							high_i_high_c.append([tabs[1],[tabs[0]]])
							s_new=False

				elif i==1: #high_i_low_c
				 	t=0
				 	while t<len(high_i_high_c) and tabs[1]!=tabs[0]: 
				 		if tabs[0] == high_i_high_c[t]: #this checks to see if the effector pair were in h_i_h_c. If it was, it will be ignored
				 			if tabs[1] == high_i_high_c[t][1]: 
				 				break
				 			else:
				 				t=len(high_i_high_c)
				 				break
				 		t+=1																				#been added to list, prevents repeats
				 	if t==len(high_i_high_c) and float(tabs[2]) >= PERC_IDENTITY_THRESH and (c_score) >= (low_length_thresh) and q_score>low_length_thresh and s_score>low_length_thresh:
				 		x=0
						q_new=True
						s_new=True
				 		while x<len(high_i_low_c):
				 			if q_new!=False and tabs[0] == high_i_low_c[x][0] : 																		 
				 				high_i_low_c[x][1].append(tabs[1]) #add all BLAST-associated hits to this entry
				 				q_new=False
				 			if s_new != False and tabs[1] == high_i_low_c[x][0]: 
				 				high_i_low_c[x][1].append(tabs[0]) 
				 				s_new=False
				 			if q_new== False and s_new==False:
				 				break
				 			x+=1
				 		if q_new==True:
				 			high_i_low_c.append([tabs[0],[tabs[1]]])
				 			q_new =False
				 		if s_new==True:
				 			high_i_low_c.append([tabs[1],[tabs[0]]]) 
				 			s_new =False

				
					
				elif i==2: #low_i_high_c
					t=0
					while t<len(high_i_high_c) and tabs[1]!=tabs[0]:
						if tabs[0] == high_i_high_c[t]:
							if tabs[1] == high_i_high_c[t][1]: #this checks to see if the effector pair were in h_i_h_c. If it was, it will be ignored
								break
							else:
								t=len(high_i_high_c)
								break
						t+=1
					if t==len(high_i_high_c) and (c_score) >= (LENGTH_THRESH) and float(tabs[2]) >= low_identity and q_score>LENGTH_THRESH and s_score>LENGTH_THRESH:
						x=0
						q_new=True
						s_new=True

						while x<len(low_i_high_c):
							if q_new==True and tabs[0] == low_i_high_c[x][0] : 																		 
								low_i_high_c[x][1].append(tabs[1]) #add all BLAST-associated hits to this entry
								q_new=False
							if s_new == True and tabs[1] == low_i_high_c[x][0]: 
								low_i_high_c[x][1].append(tabs[0]) 
								s_new=False
							if q_new== False and s_new==False:
								break
							x+=1
						if q_new==True:
							low_i_high_c.append([tabs[0],[tabs[1]]])
							q_new =False
						if s_new==True:
							low_i_high_c.append([tabs[1],[tabs[0]]])
							s_new=False
					

				elif i==3: #singletons
					x=0
					
					if lihc==True:
						while x<len(low_i_high_c): #checks through all of the clusters to make sure the effector is not present in them
							if tabs[0]  in low_i_high_c[x] :
								q_new=False
							if tabs[1] in low_i_high_c[x]:
								s_new=False
								
							x+=1
							if q_new== False and s_new==False:
								x==len(low_i_high_c)
						x=0
					else:
						while x<len(high_i_high_c):
							if tabs[0]  in high_i_high_c[x]:
								q_new=False

							if tabs[1]  in high_i_high_c[x]:
								s_new=False
							x+=1
							if q_new== False and s_new==False:
								x==len(high_i_high_c)
						x=0
					'''
					while x<len(high_i_high_c):
						if tabs[0]  in high_i_high_c[x]:
							q_new=False

						if tabs[1]  in high_i_high_c[x]:
							s_new=False
						x+=1
						if q_new== False and s_new==False:
							x==len(high_i_high_c)
					x=0
					
					while x<len(high_i_low_c):
						if tabs[0] in high_i_low_c[x]:
							q_new=False
						if tabs[1] in high_i_low_c[x]:
							s_new=False
						x+=1
						if q_new== False and s_new==False:
							x==len(high_i_low_c)
					'''
					if q_new==True:
						if tabs[0] not in singletons:
							singletons.append(tabs[0])
							#if '01_fo_lyc_4287_MIMP_TIR_Element_41' in tabs[0]:
							#	r=[tabs[0]]
							#	print 'test', r
					if s_new==True:
						if tabs[1] not in singletons:
							singletons.append(tabs[1])
							#if '01_fo_lyc_4287_MIMP_TIR_Element_41' in tabs[1]:
							#	r=[tabs[1]]
							#	print 'test', r
				y+=1		
					
				
			if i<3:
				if i==1:
					r=0
				cluster_catagories.append(single_linkage(precluster_catagories[i],all_clusters)) #cluseters together the groups
			else:
				print "Number of singletons", len(singletons), '\n'
				cluster_catagories+=[singletons]
				with open("all_clusters_test.txt", "a") as file:
					file.write("singletons\n")
					#print singletons[0]
					for s in singletons:
						file.writelines(s+'\n')
			i+=1	

	return cluster_catagories
def single_linkage(node_partners,all_clusters):
		
	clusters=[]
	i=0
	nodes=[]
	repeats=[]
	with open(all_clusters, "a") as file:
		z=1
		while i<len(node_partners):
			nodes.append(node_partners[i][0])
			i+=1
		i=0 #sets nodes to all of the keys, which should be every sample with a match
		while i< len(nodes): #cycles each of the nodes and tries to cluster them
			no_repeated_clusters=True
			x=0
			while x<len(clusters):
				if nodes[i] in clusters[x]:
					no_repeated_clusters=False
					break
				x+=1

			if no_repeated_clusters==True:
				cluster =[]
				cluster=update_cluster(node_partners, nodes[i]) #sends all of the master dictionary, the key, and an empty set to update_cluster
				file.write("cluster_"+str(z)+'\n')
				for c in cluster:
					file.writelines(c+'\n')
				file.write('\n')
				z+=1
				clusters.append(cluster)
				
			i+=1
		file.close()
	print "number of clusters: ", len(clusters)
		
	return clusters

def update_cluster(node_partners, temp):
	todo=[]
	todo.append(temp)
	cluster =[]
	while True:
		new_todo=[]
		for t in todo: #itterates through all of the keys, clusters together all of the ones that can be linked together
			partners=[]
			i=0
			while i<len(node_partners):
				if node_partners[i][0]==t:
					partners+=node_partners[i][1] #sets variable partners to all of the deffinitions for said key
					break
				i+=1
			if isinstance(partners, list) and len(partners)>1:
				for sample in partners:
					if sample not in new_todo:
						new_todo.append(sample) #forms one big array of all unique deffinitions for all keys
			elif isinstance(partners, list):
				new_todo+=partners
			cluster.append(t)
				
		x=0
		while x<len(cluster):
			if cluster[x] in new_todo:
				new_todo.remove(cluster[x]) #removes the key from the new_todo array
			x+=1
		if len(new_todo)>0: 
			todo=new_todo
		else:
			break

	return cluster

class clusterApp():

	def __init__(self):
		self.verbose = False


	def start(self, infile, blastdatabasedir, BLASTbindir,check_n, force, lihc, hilc, expanded, PERC_IDENTITY_THRESH=90.0,leave_put_eff_identifiers_during_clustering="TRUE", E_VALUE_THRESH=.001,examin="", low_length_thresh=.5, low_identity=70, LENGTH_THRESH=.9, n_allowed=0, threads=1,alignments=250): #, all_data="False"):
		
		####### ADDD PROPER NAME
		all_clusters='all_clusters_run_1.txt'
		if force!=True:
			x=2
			while os.path.isfile(all_clusters) ==True:
				all_clusters=all_clusters.split('run_')[0]+'run_'+str(x)+'.txt'
				x+=1
		with open(all_clusters, 'w') as file:
			file.close()

		if '../' in infile or infile[0]!='/':
			dir = os.path.dirname(__file__)
			infile = os.path.join(dir, infile)

		if '../' in BLASTbindir or BLASTbindir[0]!='/':
			dir = os.path.dirname(__file__)
			BLASTbindir = os.path.join(dir, BLASTbindir)

		if '../' in blastdatabasedir or blastdatabasedir[0]!='/':
			dir = os.path.dirname(__file__)
			blastdatabasedir = os.path.join(dir, blastdatabasedir)

		infile_filename = infile.split('/')[-1]

		infiledir = infile.split(infile_filename)[0]
		outfile = infile.split('.fa')[0]
		to_write=[]
		cluster_present=[]
		species_list=[]
		
		clusters_to_examine=[]
		clusters_to_write=""

		#test=True

		if examin!= "": #retrieve effectors to look at if they exist
			with open(examin, 'r') as file:
				temp=file.readlines()
				for sample in temp:
					if 'down' in sample:
						line=(sample.split("down")[0])
						if '.' in line:
							clusters_to_examine.append(line.split(".")[-1])
							#print sample.split(".")[-1]
						else:
							clusters_to_examine.append(line)
					elif 'up' in sample:
						line = (sample.split("up")[0])
						if '.' in line:
							clusters_to_examine.append(line.split(".")[-1])
						else:
							clusters_to_examine.append(line)
					else:
						clusters_to_examine.append(sample.split("\n")[0].replace('.',''))

				file.close()

		if check_n==True:
			remove_n(infile, n_allowed)
			infile=infiledir+'edited_'+infile_filename #to reflect the edited fasta file where n's are removed

		print '\n'
		print "// Running clustering script on file %s" % infile

		#print hilc, lihc

		catagory = cluster_homologous_effectors(threads, infile, E_VALUE_THRESH, PERC_IDENTITY_THRESH,LENGTH_THRESH, blastdatabasedir, BLASTbindir, infiledir, low_identity, low_length_thresh, lihc, hilc,all_clusters,alignments)

		print '\n', "Generating cluster information \n"


		blastoutfilename = infiledir+infile.split('/')[-1].split('.fa')[0]+'.vs.'+infile.split('/')[-1].split('.fa')[0]+'.blastout'
		blastout = open(blastoutfilename).readlines()
		
		##### check if files exist
		high_quality_clustered_genes_file="high_quality_clustered_genes_run_1.fasta"
		expanded_clusters_file="expanded_clusters_run_1.txt"
		all_genes_file='all_clustered_genes_run_1.fasta'
		pres_abs_file='cluster_pres_abs_run_1.txt'

		if force!=True:
			x=2
			#print "test"
			#if os.path.isfile(high_quality_clustered_genes_file)==True:
				#print "test"
			while os.path.isfile(high_quality_clustered_genes_file)==True and os.path.isfile(expanded_clusters_file)==True and os.path.isfile(all_genes_file)==True and os.path.isfile(pres_abs_file)==True:
				high_quality_clustered_genes_file=high_quality_clustered_genes_file.split('run')[0]+'_run_'+str(x)+".fasta"
				expanded_clusters_file=expanded_clusters_file.split('run')[0]+'_run_'+str(x)+'.txt'
				all_genes_file=all_genes_file.split('run')[0]+'_run_'+str(x)+".fasta"
				pres_abs_file=pres_abs_file.split('run')[0]+'_run_'+str(x)+'.txt'
				x+=1
			


		#####

		all_expanded=[]
		i=0
		d=0
		if expanded:
			expanded_clusters=''
		#debug=True
		file = list(SeqIO.parse(infile, "fasta"))
		for clusters in catagory: #catagories -> clusters -> cluster -> effector
			longest_elements = []
			cluster=0 
			while cluster<len(clusters): #cycles through each cluster and then each effector in the cluster. The longest effector will
				longest_elem = [0, '']	 #represent the cluster as a whole
				temp_array=[""]


				if i==3: #retrieves the sequence for singeltons
					for sample in file:
						if clusters[cluster] == sample.id:
							#if debug==True:
							#	print sample
							#	debug=False
							longest_elements.append(sample)
							all_expanded.append('\n--------singleton_'+sample.id+'\n'+str(sample.id)+'\n')
							if examin!='':
								for effector in clusters_to_examine:
									if d<5:
										print effector, '\n', clusters[cluster]
										d+=1
									if effector in clusters[cluster]:
										clusters_to_write+=">"+str(sample.id)+"\n"+str(sample.seq)+"\n\n"
										#print test

										break

							break 
					try:
						temp_split= clusters[cluster].split('.', 1)[1].split('_len', 1)[0].split('_')
					except:
						temp_split= clusters[cluster]
					temp_species=" "
					z=0
					while z < len(temp_split):	
						if z!=0 and z !=1 and z != len(temp_split)-1:
							temp_species+=temp_split[z]+"_"
						z+=1
					if temp_species not in species_list:
						species_list.append(temp_species)

					temp_array=[str(clusters[cluster])]
					for species in species_list:
						if species == temp_species:
							temp_array.append(1)
						else:
							temp_array.append(0)
					cluster_present.append(temp_array)
					
					
				else:  #Things that are not singletons
					if expanded:
						expanded_clusters=''
					elm_in_cluster=[]

					for elem in clusters[cluster]: #clusters are an array of samples
						if expanded:
							expanded_clusters+=str(elem)+'\n'

						if examin!=[]:  #if there are clusters to examine, this will will add all elements of the cluster to a list
							elm_in_cluster.append(elem)

						
						elem_len =0
				

						for sample in file:
							if elem == sample.id:
								elem_len = len(sample.seq)
								break 

						if elem_len > longest_elem[0]: #the longest sequence of the cluster will be what represents the cluster
							longest_elem[0] = elem_len
							longest_elem[1] = elem

						try:
							temp_split= elem.split('.', 1)[1].split('_len', 1)[0].split('_')
						except:
							try:
								temp_split=elem.split('\n')[0]
							except:
								temp_split=elem

						temp_species=""
						z=0
						while z < len(temp_split):
							if z!=0 and z != len(temp_split)-1:
								temp_species+=temp_split[z]+"_"
							z+=1
						if temp_species not in species_list:
							species_list.append(temp_species)

						
						z=1
						for species in species_list: #this will write a pres/abs sheet for the clusters that will tell how many effectors comprise a 
							if species in temp_species: #cluster and where those effectors are from
								if len(temp_array)<=z: 
									temp_array.append(1)
								else:
									temp_array[z]+=1
							else:
								if len(temp_array)<z: 
									temp_array.append(0)
							z+=1

					#####################################################################################################
					if expanded:
						#cmnd = BLASTbindir+"/blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -query "+infile+' -db '+database_store+' -out '+blastoutfilename+' -evalue '+str(E_VALUE_THRESH)
						#											0	  1	    2	   3	  4 		5	    6	   7	8	  9 	10	 11	   12	 13
						difference=[100,'',0,'',' gene(s) do not match: '] #difference identity, coverage
						new_cluster=[]
						genes_not_present=0
						for header in clusters[cluster]: ########### Refine the cluster
							matches_longest=False
							longest_check =True

							for line in blastout:
								if longest_elem[1] in line and header in line and longest_elem[1] != header and header!='\n':
									matches_longest=True
									tabs = line.strip().split('\t')	
									q_score =abs((int(tabs[7]) - int(tabs[6]))/float(tabs[12])) #covereage calculated by q_score/s_score ratio
									s_score =abs((int(tabs[9]) - int(tabs[8]))/float(tabs[13])) 
									if q_score>=s_score:
										c_score=(s_score/q_score)
									else:
										c_score=(q_score/s_score)

		
									if longest_elem[1] in tabs[0] and header in tabs[1]:
										if float(tabs[2]) >= PERC_IDENTITY_THRESH and (c_score >= (LENGTH_THRESH) and c_score <= (2-LENGTH_THRESH)):
											if (int(tabs[8])>int(tabs[9]) and longest_check==True) or (int(tabs[8])<int(tabs[9]) and longest_check==True and tabs[0]+'_Reversed' in expanded_clusters):
												temp_species=header+'_Reversed'

												if temp_species+'\n' not in expanded_clusters and tabs[0]!=tabs[1]:
													try:
														expanded_clusters=expanded_clusters.split(header+'\n')[0]+temp_species+'\n'+expanded_clusters.split(header+'\n')[1]
													except:
														print header, '\n', '\n', expanded_clusters, '\n', i, '\n'
									
											else:
												temp_species=header

											longest_check=False #In blast results, the longest match will be first. This should tell orientation 

										#else:
										#	if tabs[6]>tabs[7] and header+'_Reversed' not in expanded_clusters:
										#		temp_species=header+'_Reversed'
										#		expanded_clusters=expanded_clusters.split(header)[0]+temp_species+expanded_clusters.split(header)[1]

										#	else:
										#		temp_species=header

											if abs(c_score-1)>= difference[2]:
												difference[2] = abs(c_score-1)
												difference[3] = temp_species


											if float(tabs[2]) <= difference[0]:
												difference[0] = float(tabs[2])
												difference[1] = temp_species
							
							if matches_longest==False and header.replace('\n','') not in longest_elem[1]:
								genes_not_present+=1
								new_cluster.append(header)
								difference[4]+=header+' '

					#####################################################################################################
					if expanded and longest_element[1] !='':
						string=''
						if i==0:				############## Write expanded to all_expanded
							string='\n--------h_i_h_c_'+longest_elem[1]+" smallest shared identity: "+str(difference[0])+' '+difference[1]+" largest difference in coverage: "+str(difference[2]*100)+'% '+difference[3] 
						elif i==1:
							string='\n--------h_i_l_c_'+longest_elem[1]+" smallest shared identity: "+str(difference[0])+' '+difference[1]+" largest difference in coverage: "+str(difference[2]*100)+'% '+difference[3]	
						elif i==2:
							string='\n--------l_i_h_c_'+longest_elem[1]+" smallest shared identity: "+str(difference[0])+' '+difference[1]+" largest difference in coverage: "+str(difference[2]*100)+'% '+difference[3]
						if difference[4] == ' gene(s) do not match: ':
							expanded_clusters=string+'\n'+expanded_clusters
						else:
							expanded_clusters=string+' '+str(genes_not_present)+difference[4]+'\n'+'\n'+expanded_clusters
						all_expanded.append(expanded_clusters)
					cluster_present.append(temp_array)
		
					for sample in file: #add longest element to an output list
						if sample.id ==longest_elem[1]:
							longest_elements.append(sample)
							break
						
					if examin!='': # this will check to see if any of the effectors in the clusters are something that should be looked 
						loop=True  				# at further. If it should, will output the sequences of the effectors in fasta format
						while loop ==True:
							for effector in clusters_to_examine:	
								for elm in elm_in_cluster:
									if effector in elm:
										clusters_to_write+="--- "+effector+" ---\n"
										for sample in file:
											for r in elm_in_cluster:
												if sample.id==r:
													clusters_to_write+="\n>"+str(sample.id)+"\n"+str(sample.seq)
										clusters_to_write+="\n\n"
										loop =False
										print 'test'
							loop=False

				cluster+=1

			for z, longest_element in enumerate(longest_elements): #Add ID names to clusters
				if leave_put_eff_identifiers_during_clustering == 'TRUE':
					if i==0 and 'h_i_h_c_' not in longest_elements[z].id:
						longest_elements[z].id = "h_i_h_c_"+longest_elements[z].id #in case you're clustering a clustered file (that does not contain any description, but only a .id)
					elif i==1 and 'h_i_l_c_' not in longest_elements[z].id:
						longest_elements[z].id = "h_i_l_c_"+longest_elements[z].id
					elif i==2 and 'l_i_h_c_' not in longest_elements[z].id:
						longest_elements[z].id = "l_i_h_c_"+longest_elements[z].id
					elif i==3 and 'singleton_' not in longest_elements[z].id:
						longest_elements[z].id = "singleton_"+longest_elements[z].id 
				else:
					try:
						if i==0 and 'h_i_h_c_' not in longest_elements[z].id: #formats the longest elements
							longest_elements[z].id = "h_i_h_c_"+longest_element.description.split('.', 1)[1].split('_len', 1)[0] +"\n"
						
						elif i==1 and all_data==True and 'h_i_l_c_' not in longest_elements[z].id:
							longest_elements[z].id = "h_i_l_c_"+longest_element.description.split('.', 1)[1].split('_len', 1)[0]+"\n"
						elif i==2 and all_data ==True and 'l_i_h_c_' not in longest_elements[z].id:
							longest_elements[z].id = "l_i_h_c_"+longest_element.description.split('.', 1)[1].split('_len', 1)[0]+"\n"						
						elif i==3 and 'singleton_' not in longest_elements[z].id:
							longest_elements[z].id = "singleton_"+longest_element.description.split('.', 1)[1].split('_len', 1)[0]+"\n"

						else:
							print "Error"
					except:
						if i==0 and 'h_i_h_c_' not in longest_elements[z].id: #formats the longest elements
							longest_elements[z].id = "h_i_h_c_"+longest_element.description +"\n"
						
						elif i==1 and all_data==True and 'h_i_l_c_' not in longest_elements[z].id:
							longest_elements[z].id = "h_i_l_c_"+longest_element.description+"\n"
						elif i==2 and all_data ==True and 'l_i_h_c_' not in longest_elements[z].id:
							longest_elements[z].id = "l_i_h_c_"+longest_element.description+"\n"
						
						elif i==3 and 'singleton_' not in longest_elements[z].id:
							longest_elements[z].id = "singleton_"+longest_element.description+"\n"
							if longest_elements[z].id =='01_fo_lyc_4287_MIMP_TIR_Element_41':
								print 'writing error'
						else:
							print "Error"
				longest_elements[z].description = ""
			to_write+=longest_elements
			
			i+=1

		a=0
		with open(high_quality_clustered_genes_file, 'w') as outfilewriter:		
			while a<len(to_write): #outputs the clusters
				try:
					if "h_i_h_c_" in str(to_write[a].id) or "singleton_" in str(to_write[a].id):
						outfilewriter.write(">"+str(to_write[a].id)+"\n")
						outfilewriter.write(str(to_write[a].description))
						outfilewriter.write(str(to_write[a].seq)+"\n")
				except:
					if "h_i_h_c_" in str(to_write[a].id) or "singleton_" in str(to_write[a].id):
						outfilewriter.write(">"+str(to_write[a].id))
						outfilewriter.write(str(to_write[a].seq)+"\n")	
				a+=1
			outfilewriter.close()
		
		a=0
		with open(all_genes_file , 'w') as file:	
			while a<len(to_write): #outputs the clusters
				try:
					file.write(">"+str(to_write[a].id)+"\n")
					file.write(str(to_write[a].description))
					file.write(str(to_write[a].seq)+"\n")
				except:
					file.write(">"+str(to_write[a].id))
					file.write(str(to_write[a].seq)+"\n")	
				a+=1
			file.close()

		#outfilewriter= open(pres_abs_file, 'w') #writes the pres/ abs table
		x= len(species_list)+1
		temp_string ="Cluster vs species\t"
		for species in species_list:
			temp_string+=str(species)+"\t"
		outfilewriter.write(temp_string+"Total number of effectors in the cluster \n")
		for sample in cluster_present:
			effector_total=0
			while len(sample)<x:
				sample.append("0")
			temp_string =""
			for effector_count in sample:
				if isinstance(effector_count, int):
					effector_total+=effector_count
				temp_string+=str(effector_count)+"\t"
			temp_string+=str(effector_total)+"\n"
			outfilewriter.write(temp_string)
		outfilewriter.close()

		if examin!="": #outputs the examined clusters if there were any
			with open('clusters_examined.txt', 'w') as outfilewriter:
				outfilewriter.writelines(clusters_to_write)
				outfilewriter.close()
		
		if expanded:	
			with open(expanded_clusters_file,'w') as file:
				file.writelines(all_expanded)
				file.close()

		print '-'*30
		print '-'*30
		
	
			###76_fo_lag_Lag11_MIMP_TIR_20_elements
			###76_fo_lag_Lag11_MIMP_TIR_20_elements.fasta
			