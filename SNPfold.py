"""
SNPfold.py
Authors: Matt Halvorsen (Original code, definition packaging/annotating 12-15-2010)
		 Sam Broadaway  (Modifications)
		 J.S. Martin	(Modifications 7-12-2010)
"""

import os
import string
import sys
import numpy
import random

def main(argv):
	import getopt
	argv=argv[1:]
	###check to see if output dir exists.. if it doesn't, it is created
	haveOutputDir=os.path.isdir('output')
	if haveOutputDir==False:
		os.system('mkdir output')
	#establish acceptable nts in seq and mutations listed
	nts=['A','T','U','G','C','I']
	SNPs = []
	
	#INITIALIZATION
	accurateCalc=False
	DirName=None
	help=None
	thingsToSave=None
	
	
	
	try:
		#get options and arguements from user input
		opts, args = getopt.getopt(argv, "ao:n:h",\
		["accurate","output=","nameDir=","help"])
		#print 'options : '+str(opts)
		#print "arguements : "+str(args)
	except getopt.GetoptError:
		#halt program if invalid options are passed by user
		usage()
		sys.exit(2)
	
	#set paramaters based on options specified by user
	options=[]
	options_input=[]
	for item in opts:
		opt = item[0]
		arg = item[1]
	
		if opt == '-a' or opt=='--accurate':
			#if user indicates for it, turn on accurate p-val calculations
			accurateCalc = True
			if DirName == None:
				DirName = True
		elif opt == '-n' or opt=='--nameDir=':
			#if user indicates for it, give the output directory a user-specified name
			DirName = arg
		elif opt=='-h' or opt == '--help':
			#if user indicates for it, call up the directions for running SNPfold
			help = getUnixCommandOutput('more running.txt')
			print help
			sys.exit(2)
		elif opt == '-o' or opt == '--output=':
			#save output files that are specified by user
			thingsToSave = arg
			if DirName == None:
				DirName = True
			#if the argument passed by user are not valid, halt program
			if thingsToSave not in ['vals','graphs','all']:
				usage()
				sys.exit(2)
	
	#if after parameters call for the establishment of output files but a dirName is unestablished,
	#create a random 8 character long string for the dirName
	if DirName==True:
		DirName=generate_random_string(8)
	
	#obtain wildtype sequence and mutations.. halt program if one or the other isn't there
	if len(args)!=2:
		usage()
		sys.exit(2)		
	try:
		_wild_type=args[0]
		mutants=args[1]
	except:
		usage()
		sys.exit(2)
	#is 'wild_type' a FASTA filename?
	if os.path.isfile(_wild_type)==True:
		try:
			fileText=read_data(_wild_type)
			#assumes a FASTA setup:
			fileText=fileText.splitlines()
			seqInfo=fileText[0]
			wild_type="".join(fileText[1:])
		except:
			usage()
			sys.exit(2)
	#if 'wild_type' is not a filename, establish it as the wildtype sequence
	else:
		wild_type=_wild_type
	#if non-nt characters are found in wild_type sequence, halt program
	for char in wild_type:
		if char not in nts:
			usage()
			sys.exit(2)			
	#check to see if SNPs are actually SNPs:
	
	SNPs=mutants.split(':')
	for SNP in SNPs:
		subSNPs=SNP.split(',')		
		for subSNP in subSNPs:
			try:
				####RNA editing TEMPORARY SOLUTION : TREAT AN INOSINE LIKE A GUANINE
				subSNP.replace("I","G")
				####			
				origNT=subSNP[0]
				newNT=subSNP[-1]
				nts.index(origNT)
				nts.index(newNT)
				int(subSNP[1:-1])
			except:
				usage()
				sys.exit(2)				
	
	#establish a list of parameter settings to pass onto the 'calculate_data' definition
	options=[accurateCalc,DirName,thingsToSave]	
	#Carry out the SNPfold algorithm according to user-specified input and parameters
	calculate_data(wild_type, SNPs,options)
	exit(0)

def usage():
	#Standard output upon detection of user input error
	print 'SNPfold [-h, --help] [-a, --accurate] \
[-n, --nameDir= outputDirName] [-o, --output= vals|graphs|all] <sequence or seq file> <mutations>'

def getUnixCommandOutput(command):
	#Get direct stdout from Unix command
	child=os.popen(command)
	data=child.read()
	err=child.close()
	if err:
		data=None
	return data

def read_data(filename):
	#read in data from a specified file
	file=open(filename,'r')
	file_contents=file.read()
	file.close()
	return file_contents


def format_to_RNA(RNA_chars):
	#seperates out newlines, replaces T's with U's in an RNA seq or a string of RNA seqs
	if isinstance(RNA_chars,str):
		RNA_chars=RNA_chars.replace('\n','')
		RNA_chars=RNA_chars.replace('T','U')
	elif isinstance(RNA_chars,list):
		count=0
		while count<len(RNA_chars):		
			RNA_chars[count]=RNA_chars[count].replace('\n','')
			RNA_chars[count]=RNA_chars[count].replace('T','U')
			count+=1
	return RNA_chars


def calculate_data(wild_type, SNPs,options):

	#unpack list of user specified parameters
	[accurateCalc,dirName,thingsToSave]=options

			
	
	#format the RNA sequence string
	wild_type=format_to_RNA(wild_type)
	wt_partit=None
	#create directory to place user-specified results in, as well as 
	#a file containing the wildtype RNA sequence string for the job being carried out
	if dirName!=None:
		os.system('mkdir output/'+dirName)
		outfile_name='output/'+dirName+'/SNPFold_output.txt'
		os.system('echo ">SNPfold sequence : " > output/'+dirName+'/SNPfold_sequence.txt')
		os.system('echo "'+wild_type+'" >> output/'+dirName+'/SNPfold_sequence.txt')
	out_data=[]
	#format the RNA SNPs
	SNPs=format_to_RNA(SNPs)
	if accurateCalc==True:
		headerLine='Mutation corr_coeff rank acc_p_val'
		#remove 'accurateCalc' from options list, as we no longer need it
		options=options[1:]
		#find corr_coeffs, ranks, and accurate p-values for each user-specified SNP in the RNA
		(corr_coeffs,ranks,p_values)=get_accurate_p_value(wild_type, SNPs, \
		'all_corr_coeffs.txt',options)
		
		for SNP,corr_coeff,rank,p_value in zip(SNPs,corr_coeffs,ranks,p_values):
			out_data.append([SNP,corr_coeff,rank,p_value])
			print str(SNP)+':'+str(corr_coeff)+':'+str(rank)+':'+str(p_value)
		os.system("rm dot.ps rna.ps")
	else:			
		#default SNPfold mode (p-values are estimated)
		headerLine='Mutation corr_coeff rank est_p_val'
		#remove 'accurateCalc' from options list, as we no longer need it
		options=options[1:]
		#find corr_coeffs, ranks, and estimated p-values for each user-specified SNP in the RNA

		out_data=get_estimate_p_value(wild_type, SNPs, options)


	#if indicated to do so, write the results to an output 
	#directory of either random or user-specified name
	if dirName!=None:
	
		write_output_to_file(outfile_name,headerLine,out_data)

		#tell user where to find output files, if they decided to save any
		print
		print 'output files sent to location : '
		print 'output/'+str(dirName)+'/ '
		print 'in your SNPfold source folder'
		print
	
	exit(0)
	
def get_estimate_p_value(wild_type, SNPs, options):
	out_data=[]
	corr_coeffs=ranks=p_values=[]
	[dirName,thingsToSave]=options
	wt_matrix = get_matrix('wildtype',"dot.ps",wild_type)
	wt_partit = get_partit(wt_matrix)
	#Did user specify they wished to save output files?
	if thingsToSave=='vals':
		#save wt matrix and basepairing prob textfiles
		save_matrix(wt_matrix,"output/"+dirName+"/wildtype.txt")
		save_Basepairing_Probabilities("output/"+dirName+"/wildtype_bpProbs.txt",wild_type,wt_partit)
	elif thingsToSave=='graphs':
		#save wt matrix image
		make_matrix_image(wt_matrix,"output/"+dirName+"/wildtype.eps")
	elif thingsToSave=='all':
		#save wt matrix and base pairing prob textfiles, as well as wt matrix image
		save_matrix(wt_matrix,"output/"+dirName+"/wildtype.txt")
		save_Basepairing_Probabilities("output/"+dirName+"/wildtype_bpProbs.txt",wild_type,wt_partit)
		make_matrix_image(wt_matrix,"output/"+dirName+"/wildtype.eps")

	#get full RNA sequences for mutants of interest
	mutants = get_mutant_sequences(wild_type, SNPs)
	all_mutant_partits = []
	SNP_positions = []
	
	for mutant in mutants:
		#Get each mutant partition function matrix + matrix column sums
		mt_matrix=get_matrix(mutant[0],"dot.ps",mutant[1])
		mt_partit = get_partit(mt_matrix)
		all_mutant_partits.append(mt_partit)
		
		SNP_position=get_position(mutant[0])
		SNP_positions.append(SNP_position)		
		all_mutant_partits.append(mt_partit)

		#calculate correlation coefficient between wildtype and mutant part. func. column sums	
		raw_data = numpy.corrcoef(wt_partit,mt_partit)	
		corr_coeff = raw_data[0][1]
		
		#Estimate the rank and p-value of the mutation 
		rank_p_value = estimate_p_value(str(corr_coeff),len(mutant[1]))
		rank = rank_p_value[0]
		p_value = rank_p_value[1]
		
		out_data.append([mutant[0],corr_coeff,rank,p_value])
		#Did user specify they wished to save output files?
		if thingsToSave=='vals':
			#save mutant matrix and basepairing prob textfiles
			save_matrix(mt_matrix,"output/"+dirName+"/"+mutant[0]+".txt")
			save_Basepairing_Probabilities("output/"+dirName+"/"+mutant[0]+"_bpProbs.txt",mutant[1],mt_partit)
		elif thingsToSave=='graphs':
			#save matrix image, as well as graph of basepairing probs (wildtype vs. mutant)
			make_matrix_image(mt_matrix,"output/"+dirName+"/"+mutant[0]+".eps")
			make_Basepairing_Probs_graph("output/"+dirName+"/"+mutant[0]+'_bpProbs_vs_wt.eps',wt_partit,mt_partit)				
		elif thingsToSave=='all':
			#save output textfiles plus image graphs, as well as graphs
			save_matrix(mt_matrix,"output/"+dirName+"/"+mutant[0]+".txt")
			make_matrix_image(mt_matrix,"output/"+dirName+"/"+mutant[0]+".eps")
			save_Basepairing_Probabilities("output/"+dirName+"/"+mutant[0]+"_bpProbs.txt",mutant[1],mt_partit)
			make_Basepairing_Probs_graph("output/"+dirName+"/"+mutant[0]+'_bpProbs_vs_wt.eps',wt_partit,mt_partit)
	#clear out files at the end of job that aren't needed	
	os.system("rm dot.ps rna.ps estimated_p_value.txt")
	#print results on command line
	for data in out_data:
		print(str(data[0])+':'+str(data[1])+':'+str(data[2])+':'+str(data[3]))
	
	#if wt_matrix and wt_partit were overwritten, recalculate them
	if wt_partit==None:
		wt_matrix=get_matrix('wildtype',"dot.ps",wild_type)
		wt_partit=get_partit(wt_matrix)
	if (thingsToSave=='all' or thingsToSave=='graphs') and len(SNPs)>1:
		#create a graph of average base accessibility change
		#across the RNA for the mutations of interest
		make_average_change_graph(wt_partit, all_mutant_partits,SNP_positions,dirName)
		
	return out_data

def write_output_to_file(filename,headerLine,out_data):
	#writes creates a file, writes a header line, and then writes
	#to it things in the list 'out_data', where each item in the list is a line 
	os.system('echo "'+headerLine+'" > '+filename)
	for dataLine in out_data:
		#print dataLine
		dataLine=[str(item) for item in dataLine]
		dataLine=string.join(dataLine,' ')
		dataLine=dataLine.replace('     ','')
		os.system('echo "'+dataLine+'" >> '+filename)
	return
	
def get_mutant_sequences(wild_type, SNPs):
	#given a wildtype sequence and a list of formatted SNPs, returns a mutant sequence for each SNP
	#Note: format = wt_nuc+position+mut_nuc (ex: A36G, U22C, G5A)
	
	#intialize
	mutants = []
	errors = False
	bad_SNPs = []

	
	for SNP in SNPs:
		if SNP.find(',')>-1:
			#if commas are found in string, then there are multiple SNPs in listed variant..
			#split comma delimited SNPs in string into a list
			dSNP=SNP.replace(',','\n').splitlines()
			mutant=wild_type
			for subSNP in dSNP:
				#go through each sub-SNP listed in variant
				pos=int(get_position(subSNP))-1
				wt_nt = subSNP[0]
				#if wildtype seq at given position does not equal SNP inputed by user, error
				if (wild_type[pos] != wt_nt):
					#If user inputed wt_nucleotide doesn't equal input sequence, 
					#collect the faulty SNP 
					errors = True
					bad_SNPs.append(subSNP)
				mutant_nt = subSNP[len(subSNP)-1]
				#mutant sequence = wildtype (up to position) + new_nt + wildtype (after position)
				mutant = mutant[0:pos]+mutant[pos].replace(wt_nt,mutant_nt)+mutant[pos+1:]
			#put the SNP itself and mutant as an entry in a tuple
			mutants.append((SNP,mutant))
		else:
			#if a single SNP is in the listed variant
			#similiar to above, except no need to split into subSNPs
			SNP_position=get_position(SNP)
			wt_nt = SNP[0]
			pos = ""
			c_pos = 1
			while (c_pos < len(SNP) - 1):
				pos += SNP[c_pos]
				c_pos += 1
			pos = int(pos)-1
			if (wild_type[pos] != wt_nt):
				#If user inputed wt_nucleotide doesn't equal input sequence, 
				#collect the faulty SNP
				errors = True
				bad_SNPs.append(SNP)
			mutant_nt = SNP[len(SNP)-1]
			
			mutant = wild_type[0:pos]+wild_type[pos].replace(wt_nt,mutant_nt)+wild_type[pos+1:]
	
			mutants.append((SNP,mutant))
	
	if errors:
		#If errors are detected, then return error message and exit SNPfold
		message = "Error: Wild-type nucleotide in SNPs "
		count = 0
		while (count < len(bad_SNPs)-1):
			bad_SNP = bad_SNPs[count]
			message += bad_SNP+", "
			count += 1
		#report which input SNPs are incorrectly formatted
		message += bad_SNPs[-1]+" do not match wild-type sequence."
		os.system("echo "+message)
		exit(1)

	return mutants




def get_matrix(SNP,file_name,RNAseq):
	#calculate the matrix of base-pairing probabilities for a given sequence
	os.system("echo "+RNAseq+" | RNAfold -p > /dev/null")
	#get nucCOunt form the length of the RNA sequence
	nucCount = len(RNAseq)
	#create a matrix of zeros of RNA length X RNA length
	myMatrix = numpy.zeros([nucCount,nucCount])

	#get partition function matrix data from RNAfold output
	file = open(file_name,"r")
	results_data = file.read()
	file.close()
	
	#format the results data from the output file... find where the pertinent data is
	startLocation=string.find(results_data,'%data starts here')
	stopLocation=string.find(results_data,'showpage')
	results_data = results_data[(startLocation+18):stopLocation]
	rows = results_data.splitlines()

	for row in rows:
		#for every row, identify the x position, y position and partit value
		values=row.split()
		x=values[0]
		y=values[1]
		partit=values[2]
		# If the line had 'ubox' in it (which means that it is a part of the partit. func.:
		# Fill in at [position x, position y] of the matrix the value partit^2
		# Fill in at [position y, position x] of the matrix the value partit^2
		if values[3]=='ubox':
			myMatrix[int(y)-1,int(x)-1]=(float(partit)*float(partit))
			myMatrix[int(x)-1,int(y)-1]=(float(partit)*float(partit))
	

	# These values should be the same, since the matrix is hermetian
	# Need s for sum of values in columns
	s=numpy.sum(myMatrix,axis=1)
	s=numpy.sum(myMatrix,axis=0)

	# Store sums in the list named 'partit'
	partit = []
	count=0
	while(count < nucCount):
		partit.append(float(s[count]))
		count += 1
	return myMatrix

def get_partit(myMatrix):
	# These values should be the same, since the matrix is hermetian
	# Need s for sum of values in columns
	s=numpy.sum(myMatrix,axis=1)
	s=numpy.sum(myMatrix,axis=0)
	nucCount=len(s)
	

	# Store sums in the list named 'partit'
	partit = []
	count=0
	while(count < nucCount):
		partit.append(float(s[count]))
		count += 1
	return partit


	
def get_position(mutation):
	#get position from a mutation in the format wt_nuc+position+mut_nuc (ex: A36G, U22C, G5A)
	mutation=mutation.replace("A","")
	mutation=mutation.replace("U","")
	mutation=mutation.replace("T","")
	mutation=mutation.replace("C","")
	mutation=mutation.replace("G","")
	return mutation


def get_accurate_p_value(RNA, mutsOfInterest, outfile_name,options):
	#get accurate p-value for each listed mutation in the RNA of interest given 
	#the user specified options, and output files if indicated by user
	[dirName,thingsToSave]=options
	
	#create an output file for mutations looked at
	headerLine='Mutation corr_coeff'
	os.system('echo "'+headerLine+'" >'+outfile_name)

	#create a matrix of zeros of appropriate size for creating a partit func matrix for the RNA
	numpy.set_printoptions(threshold=sys.maxint)
	nucCount=len(RNA)
	myMatrix=numpy.zeros([nucCount,nucCount])
	#change seq from DNA->RNA

	count=0
	wild_type=RNA
	#do partit func calculations for wildtype seq
	wt_matrix=get_matrix('wild_type',"dot.ps",wild_type)
	wt_partit=get_partit(wt_matrix)

	if thingsToSave=='vals':
		#save wt matrix and basepairing prob textfiles
		save_matrix(wt_matrix,"output/"+dirName+"/wildtype.txt")
		save_Basepairing_Probabilities("output/"+dirName+"/wildtype_bpProbs.txt",wild_type,wt_partit)
	elif thingsToSave=='graphs':
		#save wt matrix image
		make_matrix_image(wt_matrix,"output/"+dirName+"/wildtype.eps")
	elif thingsToSave=='all':	
		#save wt matrix and base pairing prob textfiles, as well as wt matrix image
		save_matrix(wt_matrix,"output/"+dirName+"/wildtype.txt")
		save_Basepairing_Probabilities("output/"+dirName+"/wildtype_bpProbs.txt",wild_type,wt_partit)
		make_matrix_image(wt_matrix,"output/"+dirName+"/wildtype.eps")





	#create various empty lists for collecting all possible single mutations in the RNA of interest
	all_mutant_partits=[]
	all_mutant_partits_of_interest=[]
	SNP_positions=[]
	SNP_positions_of_interest=[]
	nt_list='AUGC'

	
	while (count<len(wild_type)):
		#do a partit func calculation for all possible single mutations in the RNA of interest
		variants=nt_list.replace(wild_type[count],'')
		origNT=wild_type[count]
		for variant in variants:
			if variant!=origNT:
				#for every possible single nt variant at every position, form mt_seq, get mt folding
				mutant=wild_type[0:count]+variant+wild_type[count+1:]
				mt_name=str(origNT)+str(count+1)+str(variant)
				mt_matrix=get_matrix(mt_name,"dot.ps",mutant)
				mt_partit=get_partit(mt_matrix)

				if thingsToSave in ['graphs','all']:
					# if user wished to create an average change graph, collect 
					#the mutation positions and corr coeffs into the aforementioned lists
					SNP_position=mt_name[1:-1]
					SNP_positions.append(SNP_position)
					all_mutant_partits.append(mt_partit)
					
					
				if mt_name in mutsOfInterest:
					#if, as you're going through all possible mutations, you find a mutation
					#listed in the input, then if the user indicated so, save the appropriate files
					if thingsToSave=='vals':
						save_matrix(mt_matrix,"output/"+dirName+"/"+mt_name+".txt")
						save_Basepairing_Probabilities("output/"+dirName+"/"+mt_name+"_bpProbs.txt",mutant,mt_partit)
					elif thingsToSave=='graphs':
						make_matrix_image(mt_matrix,"output/"+dirName+"/"+mt_name+".eps")
						make_Basepairing_Probs_graph("output/"+dirName+"/"+mt_name+'_bpProbs_vs_wt.eps',wt_partit,mt_partit)
					elif thingsToSave=='all':
						save_matrix(mt_matrix,"output/"+dirName+"/"+mt_name+".txt")
						save_Basepairing_Probabilities("output/"+dirName+"/"+mt_name+"_bpProbs.txt",mutant,mt_partit)
						make_matrix_image(mt_matrix,"output/"+dirName+"/"+mt_name+".eps")
						make_Basepairing_Probs_graph("output/"+dirName+"/"+mt_name+'_bpProbs_vs_wt.eps',wt_partit,mt_partit)
					
					if thingsToSave in ['graphs','all']:
						SNP_positions_of_interest.append(mt_name[1:-1])
						all_mutant_partits_of_interest.append(mt_partit)
	
				#calculate the correlation coefficient between the wildtype and
				#mutant sequence partition function matrix column sums
				raw_data = numpy.corrcoef(wt_partit,mt_partit)
				corr_coeff = raw_data[0][1]	
				#put the the corr coeff in an output file
				os.system("echo '"+mt_name+" "+str(corr_coeff)+"' >> "+outfile_name)
		count+=1
	

	
	
	# create empty lists to put the correlation coefficients, ranks, and p-values into
	corr_coeffs=[]
	ranks=[]
	p_values=[]
	for mutOfInterest in mutsOfInterest:
		#get correlation coefficients for mutants of interest from all corr. coeffs. calculated
		mutSeqOfInterest=get_mutant_sequences(wild_type,[mutOfInterest])
		formatted_mutOfInterest=mutOfInterest[1:-1]+'|'+mutOfInterest[0]+'>'+mutOfInterest[-1]
		#will work if the mutation is a single mutation
		mutOfInterest_corr_coeff=getUnixCommandOutput(\
		'grep -w "'+mutOfInterest+'" '+outfile_name)
		if mutOfInterest_corr_coeff==None:
			mt_matrix=get_matrix(mutOfInterest,"dot.ps",mutSeqOfInterest[0][1])
			mt_name=mutSeqOfInterest[0][0]
			mt_partit=get_partit(mt_matrix)
			raw_data=numpy.corrcoef(wt_partit,mt_partit)
			mutOfInterest_corr_coeff=raw_data[0][1]
			os.system("./accurate_p_value.sh "+str(mutOfInterest_corr_coeff)+" "+outfile_name+" > actual_p_value.txt")
			if thingsToSave=='vals':
				save_matrix(mt_matrix,"output/"+dirName+"/"+mt_name+".txt")
				save_Basepairing_Probabilities("output/"+dirName+"/"+mt_name+"_bpProbs.txt",mutant,mt_partit)
			elif thingsToSave=='graphs':
				make_matrix_image(mt_matrix,"output/"+dirName+"/"+mt_name+".eps")
				make_Basepairing_Probs_graph("output/"+dirName+"/"+mt_name+'_bpProbs_vs_wt.eps',wt_partit,mt_partit)
			elif thingsToSave=='all':
				save_matrix(mt_matrix,"output/"+dirName+"/"+mt_name+".txt")
				save_Basepairing_Probabilities("output/"+dirName+"/"+mt_name+"_bpProbs.txt",mutant,mt_partit)
				make_matrix_image(mt_matrix,"output/"+dirName+"/"+mt_name+".eps")
				make_Basepairing_Probs_graph("output/"+dirName+"/"+mt_name+'_bpProbs_vs_wt.eps',wt_partit,mt_partit)

			if	thingsToSave in ['graphs','all']:
				SNP_positions_of_interest.append(mutSeqOfInterest[0][0])
				all_mutant_partits_of_interest.append(mt_partit)
		else:
			mutOfInterest_position=mutOfInterest_corr_coeff.split(' ')[0][1:-1]
			mutOfInterest_corr_coeff=mutOfInterest_corr_coeff.split(' ')[1]
			mutOfInterest_corr_coeff=mutOfInterest_corr_coeff.replace('\n','')

			#run accurate_p_value.sh to find the ranking and p-value of the mutation of interest
			os.system("./accurate_p_value.sh "+str(mutOfInterest_corr_coeff)+" "+outfile_name+" > actual_p_value.txt")
			
		out_file_handle = open("actual_p_value.txt","r")
		#get rank and p-value from output of accurate_p_value.sh
		data = out_file_handle.read().split(':')
		rank = data[0]
		p_value = data[1]
		p_value = p_value.replace('\n','')
		#print rank
		#print corr_coeff
		corr_coeffs.append(mutOfInterest_corr_coeff)
		ranks.append(rank)
		p_values.append(p_value)
	if dirName!=None:
		#copy corr coeff + rank + p-value to the outfile for the job at hand
		os.system('cp '+outfile_name+' output/'+dirName)
	
	if	thingsToSave in ['graphs','all']:
		#save a graph of average absolute base accessibility change for
		#the SNPs of interest, as well as all possible single mutations
		make_average_change_graph(wt_partit, all_mutant_partits_of_interest,SNP_positions_of_interest,dirName)
		make_average_change_graph(wt_partit, all_mutant_partits,SNP_positions,dirName)
	#remove unneeded files at the end of the job
	os.system('rm actual_p_value.txt '+outfile_name)
	return (corr_coeffs,ranks,p_values)

def estimate_p_value(corr_coeff, length):
	#lengths of RNA for which there are random RNAs with pre-calculated corr coeffs to refer to
	standard_lengths = [10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,\
	450,500,550,600,650,700,750,800,850,900,950,1000,1500,2000]

	#identify which of the standard lengths is closest to the length of RNA of interest
	smallest_diff = 1000000
	closest = 0
	for standard in standard_lengths:
		diff = abs(length-standard)
		if (diff < smallest_diff):
			smallest_diff = diff
			closest = standard

	#run estimate_p_value.sh to find the rank of the corr_coeff in the
	#RNA of closest length to the RNA of interest
	os.system("./estimate_p_value.sh "+corr_coeff+" "+str(closest)+\
		" > estimated_p_value.txt")
	#get the rank and p-value from the output for estimate_p_value.sh
	out_file_handle = open("estimated_p_value.txt","r")
	data = out_file_handle.read().split(':')
	rank = data[0]
	p_value = data[1]
	p_value = p_value.replace('\n','')
	return (rank, p_value)


def make_matrix_image(myMatrix,filename):
	#Make a graphic representation of the number matrix inputted and save as filename
	name=filename.split('/')[2]
	import matplotlib.pyplot as plt
	#plt.matshow(myMatrix)
	plt.matshow(myMatrix, cmap=plt.cm.hot)
	plt.colorbar(orientation='vertical')
	plt.title('Base-pairing probabilities in '+name[:-4])
	plt.xlabel('Residue')
	plt.ylabel('Residue')	
	plt.savefig(filename,dpi=300)
	plt.clf()
	
def save_matrix(myMatrix,filename):
	#save the number matrix as a text file
	numpy.savetxt(filename,myMatrix)
	return

def save_Basepairing_Probabilities(filename,sequence,bp_probs):
	#save a list of base pairing probabilities for a given sequence as filename
	os.system('echo ">'+sequence+'" > '+filename)
	os.system('echo "'+str(bp_probs)+'" >> '+filename)
	return
	
def make_Basepairing_Probs_graph(filename,wt_bp_probs,mt_bp_probs):
	#make a graph of base-pairing probability as a function of
	#position in RNA for wildtype, mutant RNA
	import matplotlib.pyplot as plt
	x1=range(len(wt_bp_probs))
	x1=[val+1 for val in x1]
	x2=range(len(mt_bp_probs))
	x2=[val+1 for val in x2]
	plt.figure()
	plt.plot(x2,mt_bp_probs,'r',x1,wt_bp_probs,'k')
	plt.xlabel('Residue')
	plt.ylabel('Base-pairing Probability')
	plt.savefig(filename,format='eps')
	plt.clf()

def generate_random_string(stringLength):
	#generate a string of eight random characters, pulled from a string of possible characters
	allChars='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
	randomString=''
	count=0
	while count<stringLength:
		newChar=random.choice(allChars)
		randomString+=newChar
		count+=1
	return randomString
		


def make_average_change_graph(wild_type, SNPs,SNP_positions,dirName):
	#create a line graph of average absolute change in base-pairing probability
	#for a given RNA and a list of SNPs
	import matplotlib.pyplot as plt
	
	
	num_SNPs = len(SNPs)
	
	averages = []

	wt_pos = 0
	while (wt_pos < len(wild_type)):
		averages.append(0)
		c_SNP_index = 0
		while (c_SNP_index < len(SNPs)):
			c_SNP = SNPs[c_SNP_index]
			averages[wt_pos] += abs(wild_type[wt_pos]-c_SNP[wt_pos])
			c_SNP_index += 1
		wt_pos += 1

	count = 0
	while (count < len(averages)):
		averages[count] = averages[count]/num_SNPs
		count += 1

	out_file_handle = open('average_changes.txt','w')
	out_file_handle.write("Position\tAverage partition function change\n")

	#find average absolute change in base-pairing probability per position
	count = 0
	max = 0
	x=[]
	y=[]
	positions=[]
	average_bping_changes=[]
	while (count < len(averages)):
		if (averages[count] > max):
			max = averages[count]
		out_file_handle.write(str(count+1)+"\t"+str(averages[count])+"\n")
		x.append(count+1)
		y.append(averages[count])
		count += 1
	out_file_handle.close()
	#using matplotlib create a linegraph
	plt.figure()
	plt.plot(x,y,'r')
	plt.xlim(0,len(x)+1)
	plt.ylim(0,1)
	plt.xlabel('Residue')
	plt.ylabel('Average Structure Change')
	#On default, indicate the SNP positions in the graph via vertical green lines.
	#If all possible SNP positions are listed, then omit SNP position lines,
	#and name the end file accordingly
	#print SNP_positions
	#print 'afafafafaf'
	if len(SNP_positions) != len(wild_type)*3:
		for SNP_position in SNP_positions:
			if SNP_position.find(',')!=-1:
				#print SNP_position
				subSNP_positions=SNP_position.split(',')
				for subSNP_position in subSNP_positions:
					subSNP_position=subSNP_position[1:-1]
					plt.axvline(int(subSNP_position),color='g')
			else:
				plt.axvline(int(SNP_position),color='g')		
		filename= 'average_changes_in_basepairing.eps'
	else:
		filename='average_changes_in_basepairing_all_SNPs.eps'
	
	plt.savefig('output/'+dirName+'/'+filename)
	plt.clf()
	return

if __name__ == '__main__':
	#if SNPfold.py is run from the command line, the commence running the SNPfold program..
	#program can not be run via importing of the file in another python script.
	main(sys.argv)
