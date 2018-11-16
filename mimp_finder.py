import os
import subprocess

class mimp_finderApp():

	def __init__(self):
		self.verbose = False

	def start(self, directory_folder, working, seed_mimp="MIMP_TIR.fasta", cut_off=400, distance=2500):
		output_path=working+'/downstream_upstream/'
		if not os.path.exists(output_path):
			os.makedirs(output_path)
		files = os.listdir(directory_folder) #/Users/pmh/Documents/Sam_Brinker/mimp/ec
		
		os.chdir(working) #/Users/pmh/Documents/Sam_Brinker/mimp/
		
		second_run_through=False
		cheat=False
		
		print("\n----Finding mimps----\n")
		bashCommand="mkdir -p tirmite_mimps" #creates a temp directory to work in
		subprocess.check_output(['bash','-c', bashCommand])
		for file in files:

			if 'fasta' in file:
				print("\\ Processing "+str(file))
				bashCommand = "tirmite --alnFile "+seed_mimp+" --alnFormat fasta --stableReps 2 --outdir tirmite_mimps -v --maxeval 1 --maxdist 10000 --genome "+str(directory_folder)+str(file)+" --prefix "+str(file.split(".fasta")[0] +" --mincov .9")
				subprocess.check_output(['bash','-c', bashCommand])

			
				print("\\ Finding down stream / up stream regions")



				try:
					with open("tirmite_mimps/"+str(file.split(".fasta")[0])+"_MIMP_TIR_elements.fasta", 'r') as f:
						tir_elements=f.readlines()
						f.close()
				except:
					try: 
						with open("tirmite_mimps/"+str(file.split(".fasta")[0]).replace('-','')+"_MIMP_TIR_elements.fasta", 'r') as f:
							tir_elements=f.readlines()
							f.close()
					except:

						with open("tirmite_mimps/"+str(file.split(".fasta")[0]).replace('-','').replace(':','').replace('.','')+"_MIMP_TIR_elements.fasta", 'r') as f:
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

				stream_path=output_path+ str(file.split(".fasta")[0])+"_downstream_upstream.fasta"
				with open(stream_path, 'w') as stream:
					stream.close()

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
								up_down_stream.append(elem.split("_MIMP_TIR_")[0]+"_"+contig_to_examin+"_upstream_region:"+str(to_start+1)+'_'+str(start+1)+"\n")
								x=genome[i+1][to_start:start]
								up_down_stream.append(x+'\n')
							if to_end!=end:
								up_down_stream.append(elem.split("_MIMP_TIR_")[0]+"_"+contig_to_examin+"_downstream_region:"+str(end+1)+'_'+str(to_end+1)+"\n")
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

