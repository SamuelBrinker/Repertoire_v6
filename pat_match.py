import os

def pat_match(pattern, text):  #  matching algorithm, returns either '1' (match) or '0' (no match)  

	for i in range(len(pattern)):  #  compares pattern and text by individual characters
		match = True
		if pattern[i] == text[i] or pattern[i] == 'N':  #  loop continues and match remains equal to true if pattern matches  
			continue                                    #  text at given position or pattern at the given position is 'N'
		else:  #  if pattern does not equal match at given position, match is set to false and breaks out of the loop
			match = False
			break		
	if match == True:
		m = 1
	else:
		m = 0

	return m

class pat_matchApp():

	def __init__(self):
		self.verbose = False

	def start(self, inputfile, outputfile='pat_match_output.txt', working=''):

		if working!='':
			if working[-1]!='/':
				working+='/'
			if '../' in working:
				dir = os.path.dirname(__file__)
				working = os.path.join(dir, working)
			os.chdir(working)
		else:
			working=os.path.dirname(__file__)

		if '../' in inputfile or inputfile[0]!='/':
			dir = os.path.dirname(working)
			inputfile = os.path.join(dir, inputfile)


		path = os.getcwd()

		with open(inputfile, 'r') as f:
			file = f.readlines()
			f.close()


		l = len(file[0].split('\t'))  #  length variable for loop below
		gene_lst = file[0].split('\t')[2:]  #  list of gene headers in order given by file
		text_lst = []
		pattern_line = True

		for x in range(1, l):  #  parses tab deliminated text file for pattern and generates a list of texts to search 
			temp = []
			for line in file:  #  splits ever line of file into a list by tabs, assigns each item on each line at a given position to a list
				number = line.split('\t')[x] 
				if number[-1] == '\n':  #. removes \n characters
					temp.append(number[:-1])
				else:
					temp.append(number)
			temp.pop(0)  # removes title of each column 

			if pattern_line == True:     #  pattern_line == True only for first instance of this loop so pattern  
				pattern = ''.join(temp)  #  variable is only assigned once
				pattern_line = False
			else:   
				text_lst.append(''.join(temp))  #  joins list of items into a string and adds it to a list 


		match_lst = []
		matches = []


		for text in text_lst:  #  loop for matching function that compares pattern to each text in the list of texts
			match_lst.append(pat_match(pattern, text))  #  outputs a '1' or '0' and appends to a list 

		#  the final result of the above loop is a list of '0's and '1's whose position in the list 
		#  matches the the position of its corresponding gene in the list of genes	

		for x in range(len(gene_lst)):  #  indexes list of matches for instances of '1' and extracts the corresponding gene header at the
		    if match_lst[x] == 1:       #  same position in the list of gene headers.  Assigns genes to final lisst (matches) 
		    	matches.append(gene_lst[x])


		outputname = path+'/'+outputfile

		if os.path.exists(outputname) == True:
			x = 2
			while os.path.exists(outputname) == True:
				outputname = path+'/pat_match_output_'+str(x)+'.txt'
				x += 1 


		with open(outputname, 'w') as output:  #  creates a new file with every line being a gene header + '\n'
		    for gene in matches:
		    	if gene[-1] == '\n':
		    		output.write(gene)
		    	else:
		    		output.write(gene+'\n')
		    output.close()
	

