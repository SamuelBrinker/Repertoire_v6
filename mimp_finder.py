import os
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import gene_finder
#import pybedtools


class mimp_finderApp():

	def __init__(self):
		self.verbose = False

	def start(self, gene, force,run_gff, directory_folder, seed_mimps,working='', bedtools='bedtools',min_length=200,
	 max_length=400, distance=2500, maxeval=1, maxdist=10000, mincov=.9, signalP='signalP'):
		
		if working!='':
			if working[-1]!='/':
				working+='/'
			if '../' in working:
				dir = os.path.dirname(__file__)
				working = os.path.join(dir, working)
			os.chdir(working)
		else:
			working=os.path.dirname(__file__)
		#print os.getcwd()

		if '../' in directory_folder or directory_folder[0] !='/':
			#dir = os.path.dirname(__file__)
			#directory_folder = os.path.join(dir, directory_folder)
			dir = os.path.dirname(working)
			directory_folder = os.path.join(dir, directory_folder)
			directory_folder = os.path.abspath(os.path.realpath(directory_folder))

		if ('../' in bedtools or bedtools[0] !='/') and bedtools!='bedtools':
			dir = os.path.dirname(working)
			bedtools = os.path.join(dir, bedtools)

		if '../' in seed_mimps or seed_mimps[0] !='/':
			dir = os.path.dirname(working)
			seed_mimps = os.path.join(dir, seed_mimps)
			seed_mimps = os.path.abspath(os.path.realpath(seed_mimps))

		if '../' in signalP or signalP[0]!='/':
			dir = os.path.dirname(working)
			signalP = os.path.join(dir, signalP)
			signalP = os.path.abspath(os.path.realpath(signalP))


		if working[-1]=='/' :
			working=working[:-1]
		
		if directory_folder[-1]!='/':
			directory_folder=directory_folder+'/'

		models = os.listdir(seed_mimps)
		if len(models) <1:
			print("Need to a mimp model to continue")
			stop()


		files = os.listdir(directory_folder) #/Users/pmh/Documents/Sam_Brinker/mimp/ec
		#print len(files), files[2]
		number_of_models=0
		gff_files=[]
		os.chdir(working) #/Users/pmh/Documents/Sam_Brinker/mimp/
			
		second_run_through=False
		cheat=False
		
		print("\n----Finding mimps----\n")
		
		#bashCommand="mkdir -p tirmite_mimps" #creates a temp directory to work in
		#subprocess.check_output(['bash','-c', bashCommand])

		for model in models:
			if '.fasta' in model:
				number_of_models+=1
				
				output_path=working+'/'+model.split('.fasta')[0]+'/downstream_upstream/'
				if not os.path.exists(output_path):
					os.makedirs(output_path)
				os.chdir(working+'/'+model.split('.fasta')[0])

				for file in files:

					if 'fasta' in file:
						stream_path=output_path+ str(file.split(".fasta")[0])+"_downstream_upstream.fasta"
						#print stream_path
						if os.path.isfile(stream_path)==False or (os.path.isfile(stream_path)==True and force==True):
							with open(stream_path, 'w') as stream:
								stream.close()

							print(("\\ Processing "+str(file)+' with model '+model.split('.fasta')[0]))

							bashCommand = "tirmite --alnFile "+seed_mimps+model+" --alnFormat fasta --stableReps 2 --outdir "+output_path.split('downstream_upstream/')[0] + "tirmite_mimps -v --maxeval "+str(maxeval)+" --maxdist "+ str(maxdist)+" --genome "+str(directory_folder)+str(file)+" --prefix "+str(file.split(".fasta")[0] +" --mincov "+str(mincov))
							#print bashCommand
							subprocess.check_output(['bash','-c', bashCommand])

						
							print("\\ Finding down stream / up stream regions")



							try:
								with open("tirmite_mimps/"+str(file.split(".fasta")[0])+"_"+model.split(".")[0].split('/')[-1]+"_elements.fasta", 'r') as f:
									tir_elements=f.readlines()
									f.close()
							except:
								#print file
								try: 
									with open("tirmite_mimps/"+str(file.split(".fasta")[0]).replace('-','')+"_"+model.split(".")[0].split('/')[-1]+"_elements.fasta", 'r') as f:
										tir_elements=f.readlines()
										f.close()
								except:

									with open("tirmite_mimps/"+str(file.split(".fasta")[0]).replace('-','').replace(':','').replace('.','')+"_"+model.split(".")[0].split('/')[-1]+"_elements.fasta", 'r') as f:
										tir_elements=f.readlines()
										f.close()
							

							elements_to_examin=[]
							for elem in tir_elements:
								if '>' in elem:
									if int(elem.split('len=')[1].split("]")[0]) <=max_length and int(elem.split('len=')[1].split("]")[0]) >= min_length:
										elements_to_examin.append(elem)

							up_down_stream=[]
							#genome=[]

							#with open(directory_folder+file, 'r') as f:
							#	genome=f.readlines()
							#	f.close()
							genome = list(SeqIO.parse(str(directory_folder)+str(file), "fasta"))
							i=0


							for contig in genome:
								for elem in elements_to_examin:
									contig_to_examin=''.join(elem.split('[')[1].split(':')[:-1])
									
									if contig_to_examin in contig.id and (len(str(contig.id).split(contig_to_examin))==1 or str(contig.id).split(contig_to_examin)[1] ==',' or str(contig.id).split(contig_to_examin)[1] ==''):
										#print 'test'
										if elem.count(':') >=2:
											start = int(elem.split(":")[-1].split("_")[0])
											end = int(elem.split(":")[-1].split("_")[1].split(" ")[0])
										else:
											start = int(elem.split(":")[1].split("_")[0])
											end = int(elem.split(":")[1].split("_")[1].split(" ")[0])

										if start-distance>=1:
											to_start=start-distance
										else:
											to_start=1		

										if end+distance<=len(str(contig.seq)):
											to_end=end+distance
										else:
											to_end=len(str(contig.seq))-1					
										x='a'
										if start!=to_start:
											up_down_stream.append(elem.split(model.split(".")[0].split('/')[-1])[0]+"_"+contig_to_examin+"_upstream_region:"+str(to_start+1)+'_'+str(start+1)+"\n")
											x=str(contig.seq)[to_start:start]
											up_down_stream.append(x+'\n')
										if to_end!=end:
											up_down_stream.append(elem.split(model.split(".")[0].split('/')[-1])[0]+"_"+contig_to_examin+"_downstream_region:"+str(end+1)+'_'+str(to_end+1)+"\n")
											up_down_stream.append(str(contig.seq)[end:to_end]+'\n')
										sequence=''
										if x=="":
											print("Error")
											print(to_start, start, "1")
											print(end, to_end, "2")
											print(contig.id, contig_to_examin, "3")
											print(str(contig.seq)[1:3], "4")
											print(len(contig.seq), '5')
											#print genome[i+1][start]
											#print genome[i+1][end]
											#print genome[i+1][to_end]
											exit()
										

										with open(stream_path, 'a') as stream:
											stream.writelines(up_down_stream)
											up_down_stream=[]
											stream.close()
									
								i+=1
					# End of looping through the genome files
				
				##################
				print("\\ Combining found elements for model "+model.split('.fasta')[0])
				#bashCommand = "cat "+os.getcwd()+"/tirmite_mimps/*elements*.fasta > "+os.getcwd() +"/tirmite_mimps/all_elements.fasta"
				#subprocess.check_output(['bash','-c', bashCommand])
				elements= os.listdir(os.getcwd() +"/tirmite_mimps/")
				combined_elements=[]
				for element in elements:
					if 'elements' in element and 'all_elements.fasta' not in element and '.gff' not in element:
						with open(os.getcwd() +"/tirmite_mimps/"+element) as e:
							combined_elements +=e.readlines()
							e.close()
						#print element
						#combined_elements+=t
				#print [elements]

				with open(os.getcwd() +"/tirmite_mimps/all_elements.fasta",'w') as a:
					a.writelines(combined_elements)
					a.close()

				gff=[]
				print("\\ Producing GFF file")
				for line in combined_elements:
					if ">" in line:
						try:
							data=''.join(line.split('[')[1].split(':')[:-1])+'\t'+model.split('.fasta')[0]+'\tmimp\t'+str(int(line.split(":")[-1].split('_')[0]))+'\t'+ str(int(line.split(":")[-1].split('_')[1].split(' ')[0]))+'\t.\t0\t;'+'\n'
						except:
							print(line)
							data=''.join(line.split('[')[1].split(':')[:-1])+'\t'+model+'\tmimp\t'+str(int(line.split(":")[-1].split('_')[0]))+'\t'+ str(int(line.split(":")[-1].split('_')[1].split(' ')[0]))+'\t.\t0\t;'+'\n'

						gff.append(data)
				with open(os.getcwd()+"/tirmite_mimps/all_elements_"+model.split('.fasta')[0]+".gff","w") as file:
					file.writelines(gff)
					gff_files.append(os.getcwd()+"/tirmite_mimps/all_elements_"+model.split('.fasta')[0]+".gff")
					file.close()

		#End of looping through models
		#print gff



		######################
		if run_gff != False:
			print("\\ Finding unique TIR elements")
			if not os.path.exists(working+'/unique_TIR_elements/'):
				os.makedirs(working+'/unique_TIR_elements/')
			os.chdir(working+'/unique_TIR_elements/')

			index=0
			for gff in gff_files:
				#index_counter=index+1
				#while index_counter <len(gff_files):
				#print gff, gff_files[0]
				if index==0:
					pass

				elif index==1:
					unique_id=str(gff.split('/')[-1].split('.gff')[0])+'_verses_'+str(gff_files[index-1]).split('/')[-1].split('.gff')[0]+'.gff'
					bashCommand = bedtools+" intersect -v -a "+gff+' -b '+ str(gff_files[index-1])+' > '+unique_id
					subprocess.check_output(['bash','-c', bashCommand])
					with open(unique_id, 'r') as file:
						unique_seq=file.readlines()
						file.close()
					with open(gff_files[index-1], 'r') as file:
						not_unique_seq=file.readlines()
						file.close()
					with open('combined_unique_elements.gff','w') as file:
						file.writelines(not_unique_seq)
						file.writelines(unique_seq)
						file.close()
				else:
					unique_id=gff.split('/')[-1].split('.gff')[0]+'_verses_'+gff_files[index-1].split('/')[-1].split('.gff')[0]+'.gff'
					bashCommand = bedtools+" intersect -v -a "+gff+' -b combined_unique_elements.gff > '+unique_id
					subprocess.check_output(['bash','-c', bashCommand])
					with open(unique_id, 'r') as file:
						unique_seq=file.readlines()
						file.close()
					with open('combined_unique_elements.gff','a') as file:
						file.writelines(unique_seq)
						file.close()
					

				index+=1


			print('\\ Turning GFF file into a .fasta')
			#print force

			with open('combined_unique_elements.gff','r') as file:
				unique_gff_sequences=file.readlines()
				file.close()
			files = os.listdir(directory_folder)
			up_down_stream=[]
			unique_elements=[]
			for file in files:
				if '.fasta' in file:
					print('\\ Extracting genes from genome: '+ file.replace('.fasta',''))
					genome = list(SeqIO.parse(str(directory_folder)+str(file), "fasta"))
					i=0
					

					for contig in genome:
						for elem in unique_gff_sequences:
							contig_to_examin=elem.split('\t')[0]
							if contig_to_examin in contig.id and (len(str(contig.id).split(contig_to_examin))==1 or str(contig.id).split(contig_to_examin)[1] ==',' or str(contig.id).split(contig_to_examin)[1] ==''):
								start = int(elem.split("\t")[3])-1
								end = int(elem.split("\t")[4])-1					
								x='a'
								unique_elements.append('>'+elem.split('\t')[0]+"_"+elem.split('\t')[1]+"_region:"+str(start+1)+'_'+str(end+1)+"\n")
								x=str(contig.seq)[start:end]
								unique_elements.append(x+'\n')
								sequence=''
								unique_gff_sequences.remove(elem)
								if x=="":
									print("Error")
									print(to_start, start, "1")
									print(end, to_end, "2")
									print(contig.id, contig_to_examin, "3")
									print(str(contig.seq)[1:3], "4")
									print(len(contig.seq), '5')
									#print genome[i+1][start]
									#print genome[i+1][end]
									#print genome[i+1][to_end]
									exit()

								##################	record up_down

								if start-distance>=1:
									to_start=start-distance
								else:
									to_start=0		

								if end+distance<=len(str(contig.seq)):
									to_end=end+distance
								else:
									to_end=len(str(contig.seq))-1					
								x='a'
								if start!=to_start:
									up_down_stream.append('>'+elem.split('\t')[0]+"_upstream_region:"+str(to_start+1)+'_'+str(start+1)+"\n")
									x=str(contig.seq)[to_start:start]
									up_down_stream.append(x+'\n')
								if to_end!=end:
									up_down_stream.append('>'+elem.split('\t')[0]+"_downstream_region:"+str(end+1)+'_'+str(to_end+1)+"\n")
									up_down_stream.append(str(contig.seq)[end:to_end]+'\n')
								sequence=''
								if x=="":
									print("Error")
									print(to_start, start, "1")
									print(end, to_end, "2")
									print(contig.id, contig_to_examin, "3")
									print(str(contig.seq)[1:3], "4")
									print(len(contig.seq), '5')
									#print genome[i+1][start]
									#print genome[i+1][end]
									#print genome[i+1][to_end]
									exit()
								

			
			output_path=working+'/unique_TIR_elements/up_down_stream/'					
			if not os.path.exists(output_path):
				os.makedirs(output_path)		
			
			with open(output_path+'up_down_stream.fasta', 'w') as stream:
				stream.writelines(up_down_stream)
				stream.close()

			with open('all_unique_elements.fasta', 'a') as stream:
				stream.writelines(unique_elements)
				unique_elements=[]
				stream.close()

			if gene!=False:
				print("Running gene_finder")
				#	def start(self, force, directory_folder, output_dir, SignalPpath ='signalP', min_prot_len=10, max_d2m=2500, max_prot_len=134, SignalP_threshold=.45, working=''):

				run_gene=gene_finder.gene_finderApp()
				run_gene.start(True, output_path, working,signalP, 50, 2500, 134,.45, working)






