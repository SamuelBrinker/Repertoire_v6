import os
import subprocess

class mimp_finderApp():

	def __init__(self):
		self.verbose = False

	def start(self, force, directory_folder, working, seed_mimp="MIMP_TIR.fasta", cut_off=400, distance=2500, maxeval=1, maxdist=10000, mincov=.9):
		
		if '../' in directory_folder or directory_folder[0] !='/':
			dir = os.path.dirname(__file__)
			directory_folder = os.path.join(dir, directory_folder)

		if '../' in working or working[0] !='/':
			dir = os.path.dirname(__file__)
			working = os.path.join(dir, working)
		if '../' in seed_mimp or seed_mimp[0] !='/':
			dir = os.path.dirname(__file__)
			seed_mimp = os.path.join(dir, seed_mimp)


		if working[-1]=='/' :
			working=working[:-1]
		
		if directory_folder[-1]!='/':
			directory_folder=directory_folder+'/'

		output_path=working+'/downstream_upstream/'

		if not os.path.exists(output_path):
			os.makedirs(output_path)
		#elif force !=True:	
		#	x=2
		#	while os.path.exists(output_path)==True:
		#		output_path=working+'/downstream_upstream_number_'+x+'/'
		#		x+=1
		files = os.listdir(directory_folder) #/Users/pmh/Documents/Sam_Brinker/mimp/ec
		

		os.chdir(working) #/Users/pmh/Documents/Sam_Brinker/mimp/
		
		second_run_through=False
		cheat=False
		
		print("\n----Finding mimps----\n")
		bashCommand="mkdir -p tirmite_mimps" #creates a temp directory to work in
		subprocess.check_output(['bash','-c', bashCommand])
		for file in files:

			if 'fasta' in file:
				stream_path=output_path+ str(file.split(".fasta")[0])+"_downstream_upstream.fasta"
				
				if os.path.isfile(stream_path)==False or (os.path.isfile(stream_path)==True and force==True):
					with open(stream_path, 'w') as stream:
						stream.close()
					print("\\ Processing "+str(file))
					bashCommand = "tirmite --alnFile "+str(seed_mimp)+" --alnFormat fasta --stableReps 2 --outdir tirmite_mimps -v --maxeval "+str(maxeval)+" --maxdist "+ str(maxdist)+" --genome "+str(directory_folder)+str(file)+" --prefix "+str(file.split(".fasta")[0] +" --mincov "+str(mincov))
					subprocess.check_output(['bash','-c', bashCommand])

				
					print("\\ Finding down stream / up stream regions")



					try:
						with open("tirmite_mimps/"+str(file.split(".fasta")[0])+"_"+seed_mimp.split(".")[0].split('/')[-1]+"_elements.fasta", 'r') as f:
							tir_elements=f.readlines()
							f.close()
					except:
						#print file
						try: 
							with open("tirmite_mimps/"+str(file.split(".fasta")[0]).replace('-','')+"_"+seed_mimp.split(".")[0].split('/')[-1]+"_elements.fasta", 'r') as f:
								tir_elements=f.readlines()
								f.close()
						except:

							with open("tirmite_mimps/"+str(file.split(".fasta")[0]).replace('-','').replace(':','').replace('.','')+"_"+seed_mimp.split(".")[0].split('/')[-1]+"_elements.fasta", 'r') as f:
								tir_elements=f.readlines()
								f.close()
					

					elements_to_examin=[]
					for elem in tir_elements:
						if '>' in elem:
							if int(elem.split('len=')[1].split("]")[0]) <=cut_off:
								elements_to_examin.append(elem)

					up_down_stream=[]
					genome=[]

					with open(directory_folder+file, 'r') as f:
						genome=f.readlines()
						f.close()

					i=0


					for contig in genome:
						for elem in elements_to_examin:
							if ',' in elem:
								if 'unitig_' in elem:
									contig_to_examin = "unitig_"+elem.split("unitig_")[1].split(",:")[0]
								else:
									contig_to_examin = "contig_"+elem.split("contig_")[1].split(",:")[0]
							else:
								if 'unitig_' in elem:
									contig_to_examin = "unitig_"+elem.split("unitig_")[1].split(":")[0]
								else:
									
									try:
										contig_to_examin = "contig_"+elem.split("contig_")[1].split(":")[0]
									except:
										print elem

							if ('contig' in contig and contig_to_examin == str('contig_'+contig.split('contig_')[1].replace(',','').replace('\n',''))) or ('unitig' in contig and contig_to_examin == str('unitig'+contig.split('unitig_')[1].replace(',','').replace('\n',''))):
								start = int(elem.split('contig')[1].split(":")[1].split("_")[0])-1
								end = int(elem.split('contig')[1].split(":")[1].split("_")[1].split(" ")[0])-1

								if start-distance>=0:
									to_start=start-distance
								else:
									to_start=0
								if end+distance<=len(genome[i+1]):
									to_end=end+distance
								else:
									to_end=len(genome[i+1])-1								
								x='a'
								if start!=to_start:
									up_down_stream.append(elem.split(seed_mimp.split(".")[0].split('/')[-1])[0]+"_"+contig_to_examin+"_upstream_region:"+str(to_start+1)+'_'+str(start+1)+"\n")
									x=genome[i+1][to_start:start]
									up_down_stream.append(x+'\n')
								if to_end!=end:
									up_down_stream.append(elem.split(seed_mimp.split(".")[0].split('/')[-1])[0]+"_"+contig_to_examin+"_downstream_region:"+str(end+1)+'_'+str(to_end+1)+"\n")
									up_down_stream.append(genome[i+1][end:to_end]+'\n')
								if x=="":
									print "Error"
									print to_start, start
									print end, to_end
									print genome[i], contig_to_examin
									print genome[i+1][0:3]
									print genome[i+1][start]
									print genome[i+1][end]
									print genome[i+1][to_end]
									exit()
								

								with open(stream_path, 'a') as stream:
									stream.writelines(up_down_stream)
									up_down_stream=[]
									stream.close()
							
						i+=1

