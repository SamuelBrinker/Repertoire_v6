import os
from Bio import SeqIO
from string import printable

def reformat_file(infile, identification):
		#with open(infile, 'r') as file:
		#	raw_data = file.readlines()
		#	file.close()
		all_sequences=list(SeqIO.parse(infile, "fasta"))
		x=0
		y=1
		file_name=''

		if infile.count('.')==1:
			try:
				file_name=infile.split('/')[-1].split('.')[0]
			except:
				file_name=infile.split('.')[0]
		else:
			try:
				for z in infile.split('/')[-1].split('.')[:-1]:
					file_name+=z+'.'
				file_name=file_name[:-1]
			except:
				for z in infile.split('.')[:-1]:
					file_name+=z+'.'
				file_name=file_name[:-1]

		while x<len(all_sequences):
			header=str(all_sequences[x].id)
			stripped_header = ''.join(char for char in header if char in printable)
			if identification !='':
				new_header=identification+'_contig_'+str(y)
			else:
				new_header=file_name+'_contig_'+str(y)
			all_sequences[x].id=new_header
			all_sequences[x].description=''
			y+=1	
			x+=1

		extracted_path=infile.split(file_name)[0]+'/extracted_path/'
		if not os.path.exists(extracted_path):
				os.makedirs(extracted_path)	

		with open(extracted_path+file_name+'_reformated.fasta','w') as file:
			SeqIO.write(all_sequences, file, "fasta")
		return

class reformatApp():

	def __init__(self):
		self.verbose = False


	def start(self, infile='', folder='', identification='',working=''): 

		if working!='':
			if working[-1]!='/':
				working+='/'
			if '../' in working:
				dir = os.path.dirname(__file__)
				working = os.path.join(dir, working)
			os.chdir(working)
		else:
			working=os.path.dirname(__file__)

		if infile!='':
			if '../' in infile or infile[0]!='/':
				dir = os.path.dirname(working)
				infile = os.path.join(dir, infile)
			reformat_file(infile,identification)

		if folder!='':
			if '../' in folder or folder[0]!='/':
				dir = os.path.dirname(working)
				folder = os.path.join(dir, folder)
				folder = os.path.abspath(os.path.realpath(folder))

			files = os.listdir(folder)
			for file in files:
				try:
					reformat_file(file,identification)
				except:
					pass

		