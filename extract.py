import os
from Bio import Seq

class extractApp():

	def __init__(self):
		self.verbose = False


	def start(self, database, expanded, to_examine):
		
		infile_filename = expanded.split('/')[-1]
		infiledir = expanded.split(infile_filename)[0]
		extracted_path=infiledir+'extracted_genes/'
		if not os.path.exists(extracted_path):
			os.makedirs(extracted_path)		
		with open(to_examine, 'r') as file:
			all_examine=file.readlines()
			file.close()

		examin_headers=[]
		for line in all_examine:
			if '>' in line:
				examin_headers.append(line.split('>')[1])
		#print len(examin_headers)
		with open(expanded, 'r') as file:
			all_exapnded =file.readlines()
			file.close()
		clustered_headers=[]
		cluster=[]
		record=False
		y=0

		for header in examin_headers:
			for line in all_exapnded:
				if record==True:
					if '--------' in line:
						
						record=False
						clustered_headers.append(cluster)

						cluster=[]
					else:
						cluster.append(line)

				if ' 'in header:
					if '-'+header.split(' ')[0] in line:
						cluster.append(header.split(' ')[0])
						record=True
				else:
					if '-'+header.replace('\n','') in line:
						cluster.append(header)
						record=True
		clustered_headers.append(cluster)	

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
				all_genes.append(gene.split(' ')[0])
			else:
				all_genes.append(gene)

		#print [clustered_headers[-1][0]]
		#print [all_genes[0]]

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

						if ">"+cluster[z].split('_Reversed')[0] == all_genes[x] and cluster[z]!='\n':
							to_write+='>'+cluster[z]
							record=True
					else:
						if ">"+cluster[z].replace('\n',"") == all_genes[x] and cluster[z]!='\n':
							to_write+='>'+cluster[z]
							record=True
					x+=1
				z+=1
			with open(extracted_path+cluster[0].replace('\n','')+'_extracted.fasta','w') as file:
				file.write(to_write)
				file.close
			
	
		

