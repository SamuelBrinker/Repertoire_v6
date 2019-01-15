import os
from Bio import Seq

class extractApp():

	def __init__(self):
		self.verbose = False


	def start(self, database, expanded, to_examine):

		if '../' in database or database[0]!='/':
			dir = os.path.dirname(__file__)
			database = os.path.join(dir, database)
		if '../' in expanded or expanded[0]!='/':
			dir = os.path.dirname(__file__)
			expanded = os.path.join(dir, expanded)
		if '../' in to_examine or to_examine[0]!='/':
			dir = os.path.dirname(__file__)
			to_examine = os.path.join(dir, to_examine)

		infile_filename = expanded.split('/')[-1]
		infiledir = expanded.split(infile_filename)[0]
		extracted_path=infiledir+'extracted_genes/'
		if not os.path.exists(extracted_path):
			os.makedirs(extracted_path)	
		print extracted_path
		with open(to_examine, 'r') as file:
			all_examine=file.readlines()
			file.close()

		#print all_examine[:4]

		examin_headers=[]
		for line in all_examine:
			if '>' in line:
				examin_headers.append(line.split('>')[1])
			elif '\n'!=line:
				examin_headers.append(line)
		#print len(examin_headers)

		with open(expanded, 'r') as file:
			all_exapnded =file.readlines()
			file.close()
		#print all_exapnded[:4]
		x=0
		while x <len(all_examine):
			#print [all_examine[x]]
			all_examine[x]=str(all_examine[x]).replace('\t','').replace('\n','')

			x+=1

		#print all_exapnded[0:3]
		clustered_headers=[]
		clusters=[]
		record=False
		y=0
		t=True
		print examin_headers[:4]
		for header in examin_headers:
			
			for line in all_exapnded:
				#if t==True:
				#	print 'x'
				if record==True:
					if '--------' in line:
						
						record=False
						clustered_headers.append(clusters)

						clusters=[]
					else:
						clusters.append(line)

				if ' 'in header:
					#if t==True:
					#	t=False
					#	print '-'+header.split(' ')[0].replace('\n','')

					if '-'+header.split(' ')[0].replace('\n','') in line:
						clusters.append(header.split(' ')[0])
						record=True
				else:
					#if t==True:
					#	t=False
					#	print '-'+header.split(' ')[0].replace('\n','')

					if '-'+header.replace('\n','')  in line:
						clusters.append(header.replace("signalal","signal"))
						#print "test"
						record=True
		clustered_headers.append(clusters)	

		#print len(clustered_headers)


		with open(database,'r') as file:
			all_genes=file.readlines()
			file.close()

		record=False
		#print clustered_headers[:2]
		#print len(clustered_headers)
		
		
		to_edit = all_genes
		all_genes=[]
		for gene in to_edit:
			if ' 'in gene:
				all_genes.append(gene.split(' ')[0]) #turn this into a \n replace as well as a " " replace with _. problem is signal peptide spelled signalal
			else:
				all_genes.append(gene)

		#print [clustered_headers[-1][0]]
		print [all_genes[0]]
		print len(all_genes)

		reverse=False
		check=True
		for cluster in clustered_headers:
			z=0
			to_write=''

			while z<len(cluster):
				x=0
				
				if '_Reversed' in cluster[z]:
					reverse=True
					to_reverse=''
				while x<len(all_genes):

					if record==True and reverse==False:	
						if '>' in all_genes[x]:
							record=False
							to_write+='\n'							
						else:
							to_write+=all_genes[x]

					elif record==True and reverse ==True:	
						if '>' in all_genes[x]:
							record=False
							reverse=False
							my_seq=Seq.Seq(str(to_reverse))
							rrr=my_seq.reverse_complement()

							to_write+=str(rrr).replace('\n','')+'\n'						
						else:
							to_reverse+=all_genes[x].replace('\n','')
							
					if reverse==True:

						if ">"+cluster[z].split('_Reversed')[0] == all_genes[x].split(" ")[0] and cluster[z]!='\n':
							to_write+='>'+cluster[z]
							record=True
					else:
						if ">"+cluster[z].replace('\n',"") == all_genes[x].split(" ")[0] and cluster[z]!='\n':
							to_write+='>'+cluster[z]
							record=True
					x+=1
				z+=1
			if cluster != []:
				#print extracted_path
				with open(extracted_path+cluster[0].replace('\n','')+'_extracted.fasta','w') as file:
					file.write(to_write)
					file.close
				
	
		

