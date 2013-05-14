"""
Matt Halvorsen (Original code, definition packaging/annotating 12-15-2010)
		 Sam Broadaway  (Modifications)
		 J.S. Martin	(Modifications 2010-07-12)
		 Chas Kissick	(Modifications 2011-10-05)
"""
import getopt
import os
import string
import sys
import math
import numpy
import random
from shutil import rmtree
nts=['A','T','U','G','C','I','-']

import datetime

def main(argv):
	global wt_name
	
	argv=argv[1:]
	#check to see if output dir exists.. if it doesn't, it is created
	haveOutputDir=os.path.isdir('output')
	if haveOutputDir==False:
		os.system('mkdir output')
	#establish acceptable nts in seq and mutations listed
	nts=['A','T','U','G','C','I','-']
	SNPs = []
	
	#INITIALIZATION
	wt_name="wt"
	need_dir=False
	printPartit=False
	accurateCalc=False
	DirName=None
	help=None
	save_vals=False
	shannon=False
	windowMode=False
	
	import argparse
	opts = argparse.ArgumentParser()
	opts.add_argument("-a","--accurate",
		help="Calculate accurate p-values for mutations indicated by user.", action="store_true")
	opts.add_argument("-o","--output",choices=['partit','partition'],
		help="Specify files to be output.")
	opts.add_argument("-sh","--shannon", help="Specify files to be output.", action="store_true")
	opts.add_argument("-w","--window", nargs='?', type=int, const=150, default=0,
		help="Calculate the results based on only a window of the sequence centered around one SNP.")
	opts.add_argument("-n","--nameDir", help="Name the output directory.")
	opts.add_argument("-s","--save", help="Save all inputs and calculations.", action="store_true")
	opts.add_argument("wild_type", help="Nucleotide sequence")
	opts.add_argument("mutants", help="SNPs separated by colon, e.g. A5G:T18C")
	args = opts.parse_args()
	
	#Assign user requested options
	accurateCalc = args.accurate
	if args.output:
		if 'partit' in args.output or 'partition' in args.output:
			printPartit = True
		else:
			printPartit = False
	if args.shannon:
		shannon = True
		need_dir = True
	else:
		shannon = False
	windowMode = args.window
	if args.window:
		size = int(args.window)
		need_dir=True
	DirName = args.nameDir
	if args.save:
		save_vals=True
		need_dir=True
	_wild_type=args.wild_type
	mutants=args.mutants
					
	#if after parameters call for the establishment of output files but a dirName is unestablished,
	#create a random 8 character long string for the dirName
	if DirName==None and need_dir==True:
		DirName=generate_random_string(8)
	
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
	#make 'wild_type' uppercase
	wild_type=str.upper(wild_type)
	#if non-nt characters are found in wild_type sequence, halt program
	
	#is 'mutants' a filename?
	if os.path.isfile(mutants)==True:
		try:
			mutants=read_data(mutants)
		except:
			usage()
			sys.exit(2)
			
	#make contents of 'mutants' uppercase
	mutants=str.upper(mutants)
	#check to see if SNPs are actually SNPs:
	if wild_type.count("/")==1:
		(wt_name,mutation)=wild_type.split("/")
	if "-" in mutants and accurateCalc:
		usage(None,"Error: P values cannot be calculated with deletion or insertion mutations.")
		sys.exit(2)
	
	SNPs=mutants.split(":")
	
	#run window calculations before printing the sequence
	if windowMode:
		if len(SNPs) > 1:
			usage(None,"Error: Window mode can only be used with one SNP at a time.")
			sys.exit(2)
		#Make sure the segment is even so both sides of the window are of equal length
		if size:
			if size % 2 == 1:
				segment = size + 1
		#else:
		#	size = 150
		segment = size/2
		center = int(SNPs[0].strip('ACGTU-'))
		windowMode = SNPs[0]
		if center <= segment:
			first = 0
			last = segment*2 + 1
		elif len(wild_type)-center <= segment:
			first = len(wild_type) - segment*2 - 1
			last = len(wild_type)
			SNPs[0]=SNPs[0][0]+str(len(wild_type)-segment)+SNPs[0][-1]
		else:
			first = center - segment - 1
			last = center + segment
			SNPs[0]=SNPs[0][0]+str(segment+1)+SNPs[0][-1]
		wild_type = wild_type[first:last]
		print wild_type
	'''
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
				print "mut character not a nucleotide"
				usage()
				sys.exit(2)
	'''
	#establish a list of parameter settings to pass onto the 'calculate_data' definition
	options=[accurateCalc,DirName,save_vals,printPartit,shannon,windowMode]
	#Carry out the SNPfold algorithm according to user-specified input and parameters
	calculate_data(wild_type, SNPs,options)
	return 

def usage(faultyString=None,message=None):
	if faultyString!=None:
		print "error with : "
		print faultyString
	if message!=None:
		print message
	#Standard output upon detection of user input error
	print 'SNPfold [-h, --help] [-a, --accurate] \
[-n, --nameDir= outputDirName] [-o, --output (partition or shannon)] \
[-s, --save] [-w, --window= window size] <sequence or seq file> <mutations>'
	#exit(0)
	
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
	#separates out newlines, replaces T's with U's in an RNA seq or a string of RNA seqs
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
	[accurateCalc,dirName,save_vals,printPartit,shannon,window]=options
	#format the RNA sequence string
	wild_type=format_to_RNA(wild_type)
	#format the RNA SNPs
	SNPs=format_to_RNA(SNPs)
	check_mutant_accuracy(wild_type,SNPs)
	wt_partit=None
	
	#create directory to place user-specified results in, as well as 
	#a file containing the wildtype RNA sequence string for the job being carried out
	if dirName!=None:
		if os.path.isdir('output/'+dirName)==True:
			direxists(dirName)
		else: os.system('mkdir output/'+dirName)
		if shannon:
			os.system('mkdir output/'+dirName+'/bp_prob_derived')
		outfile_name='output/'+dirName+'/SNPFold_output.txt'
		os.system('echo ">SNPfold sequence: " > output/'+dirName+'/SNPfold_sequence.txt')
		os.system('echo "'+wild_type+'" >> output/'+dirName+'/SNPfold_sequence.txt')
	out_data=[]
	if accurateCalc==True:
		headerLine='Mutation corr_coeff rank acc_p_val'
		#remove 'accurateCalc' from options list, as we no longer need it
		options=options[1:]
		#find corr_coeffs, ranks, and accurate p-values for each user-specified SNP in the RNA
		(corr_coeffs,ranks,p_values)=get_accurate_p_value(wild_type, SNPs, \
		'output/'+dirName+'/all_corr_coeffs.txt',options)
		
		for SNP,corr_coeff,rank,p_value in zip(SNPs,corr_coeffs,ranks,p_values):
			out_data.append([SNP,corr_coeff,rank,p_value])
			print str(SNP)+':'+str(corr_coeff)+':'+str(rank)+':'+str(p_value)
		os.system("rm dot.ps rna.ps")
	else:			
		#default SNPfold mode (p-values are estimated)
		headerLine='Mutation corr_coeff'
		#remove 'accurateCalc' from options list, as we no longer need it
		options=options[1:]
		#find corr_coeffs, ranks, and estimated p-values for each user-specified SNP in the RNA
		out_data=no_p_value(wild_type, SNPs, options)

	#if indicated to do so, write the results to dir, either random or user-specified name
	if save_vals==True:
		if printPartit==False:
			write_output_to_file(outfile_name,headerLine,out_data)
	#tell user where to find output files
	if dirName!=None:
		print 'Output files are located at:\n%s'%(os.getcwd()+'/output/'+dirName)
	return
	
def indel_formatting(partit,sequence,type="wt"):
	formatted_partit=partit
	sequenceCount=0
	partitCount=0
	#KEEP TRACK OF DELETE POSITIONS
	deletes=[]
	#sequence length can't be less than partit length
	while sequenceCount<len(sequence):
		if type=="wt" and sequence[sequenceCount]=="-":
			formatted_partit.remove(partit[partitCount])
			deletes.append(sequenceCount)
			sequenceCount+=1
		elif type=="mut" and sequence[sequenceCount].islower():
			formatted_partit.remove(partit[partitCount])
			deletes.append(sequenceCount)
			sequenceCount+=1
			
		else:
			partitCount+=1
			sequenceCount+=1
	return [formatted_partit,deletes]

def window(fullsequence, coord, windowradius, windowoffset=0):
	coord=coord-1

	windowsequence=None
	if (2*windowradius>len(fullsequence)):
		usage()
	elif coord<windowradius:
		windowsequence=fullsequence[:coord+windowradius]
	elif (coord+windowradius)>len(fullsequence):
		windowsequence=fullsequence[coord-windowradius:]
	else:
		windowsequence=fullsequence[coord-windowradius:coord+windowradius+1]
	return windowsequence


def get_seq_eval(seq):
	eval_indeces=[]
	originalPositionCount=0
	positionCount=0
	mask=""
	original_seq_positionCount=0
	while positionCount < len(seq):
		if seq[positionCount]=="-":
			eval_indeces.append("D")
			originalPositionCount+=1

		elif seq[positionCount].islower():
			eval_indeces.append("I")
		else:
			eval_indeces.append(originalPositionCount)
			originalPositionCount+=1
		positionCount+=1
	return eval_indeces
	
def process_wt_seq(_wild_type):
	if (_wild_type.count("/")>1):
		#print "1"
		usage()
	elif (_wild_type.count("/")==1):
		[wild_type,wt_variant]=_wild_type.split("/")
		if wt_variant.find(":")!=-1:
			#print "2"
			usage()
		wt_variant=get_mutant_sequences(wild_type, [wt_variant])

		new_wild_type=wt_variant[0][1]
		
	else:
		wild_type=_wild_type
		new_wild_type=_wild_type
	for char in wild_type:
		if char not in nts:
			print "wt character not a nucleotide"
			#print "3"
			usage()
			sys.exit(2)

	wt_matrix = get_matrix('wildtype',"dot.ps",new_wild_type, 0)
	wt_partit = get_partit(wt_matrix)	
	#print new_wild_type
	#print wt_matrix
	#print wt_partit
	#print len(wt_partit)
	#quit()
	return [new_wild_type,wt_matrix,wt_partit]


def readjust_SNP_positions(SNPs,positionInfos):
	import re
	newSNPlist=[]
	for SNP in SNPs:
		subSNPlist=[]
		subSNPs=SNP.split(",")
		for subSNP in subSNPs:
			positionOrig=re.findall(r"\d+",subSNP)
			positionOrig=int(positionOrig[0])
			if positionOrig in positionInfos:
				positionNew=positionInfos.index(positionOrig)
				subSNP=subSNP.replace(str(positionOrig),str(positionNew))
				subSNPlist.append(subSNP)		
		subSNPsString=",".join(subSNPlist)
		newSNPlist.append(subSNPsString)

	return newSNPlist
	
	
def no_p_value(wild_type, SNPs, options):
	corr_coeffs=[]
	[dirName,save_vals,printPartit,shannon,window]=options
	
	[new_wild_type,wt_matrix,wt_partit]=process_wt_seq(wild_type)
	wild_type=wild_type.split("/")[0]
	positionInfos_wt=range(0,len(new_wild_type))	
	positionInfos=get_seq_eval(new_wild_type)
	formatted_wt_partit=list(wt_partit)	
	if save_vals==True:
		#save wt matrix and basepairing prob textfiles
		save_matrix(wt_matrix,"output/"+dirName+"/wildtype.txt")
		save_Basepairing_Probabilities("output/"+dirName+"/wildtype_bpProbs.txt",wild_type,wt_partit)
	
	if "NONE" in SNPs:
		mutants=[]
	else:
		#get full RNA sequences for mutants of interest
		mutants = get_mutant_sequences(new_wild_type, SNPs)
		all_mutant_partits = []
		SNP_positions = []

	if printPartit==False:
		out_data=[]
	else:
		formatted_wt_partit=[float("%.4f" % val) for val in formatted_wt_partit]		
		formatted_wt_partit=[str(val) for val in formatted_wt_partit]
		positionInfos_wt=[str(val) for val in positionInfos_wt]
		outputMaterial=str(string.join(formatted_wt_partit,","))+":"+str(string.join(positionInfos_wt,","))
		out_data=[['WT',outputMaterial]]

	if shannon:
		shanWt = shannon_entropy(dirName, wt_matrix, 'wt', 0)
		
	for mutant in mutants:
		(mutation,mut_seq)=mutant
		#delete wildtype partit vals where there is a deletion that occurs in mutant
		#print mut_seq
		
		positionInfosMut=get_seq_eval(mut_seq)
		mut_seq_2fold=mut_seq.replace("-","")
		while(1):
			try:
				positionInfosMut.remove("D")
			except:
				break

		masking_array=[wild_type,positionInfos,mut_seq,positionInfosMut]
		if window:
			mt_matrix=get_matrix(mutation,"dot.ps",mut_seq_2fold,1)
		else:
			mt_matrix=get_matrix(mutation,"dot.ps",mut_seq_2fold,0)
		mt_partit = get_partit(mt_matrix)
		formatted_mut_partit=list(mt_partit)
		
		if shannon:
			if window:
				shannon_entropy(dirName, mt_matrix, window, shanWt)
			else:
				shannon_entropy(dirName, mt_matrix, mutation, shanWt)
		
		partit_adjusted_wt=[]
		partit_adjusted_wt=[]
		
		masterSet=[]
		partit_adjusted_wt=[]
		partit_adjusted_mut=[]
		count=0
		while count<len(wild_type):
			if count in positionInfos and count in positionInfosMut:
				partitIndex=positionInfos.index(count)
				partit=formatted_wt_partit[partitIndex]
				partit_adjusted_wt.append(partit)
				partitIndex=positionInfosMut.index(count)
				partit=formatted_mut_partit[partitIndex]
				partit_adjusted_mut.append(partit)

				masterSet.append(count+1)			
			count+=1
		
		#print partit_adjusted_wt
		#print partit_adjusted_mut
		#print masterSet
		#print len(partit_adjusted_wt)
		#print len(partit_adjusted_mut)
		#print len(masterSet)

		#remove deletion symbols from mutant sequence before folding of seq
		###mut_seq=mut_seq.replace("-","")
		#fold mutant sequence (with insertions if any are present
		'''
		mt_matrix=get_matrix(mutation,"dot.ps",mut_seq)
		mt_partit = get_partit(mt_matrix)
		#delete partits corresponding to insertions in mutant sequence if there are any
		[formatted_mut_partit,mut_deletes]=indel_formatting(mt_partit,mut_seq,type="mut")

		if printPartit==True:
			mt_partit_2join=[str(val) for val in mt_partit]
			print str(mutant[0])+":"+string.join(mt_partit_2join,",")

		all_mutant_partits.append(mt_partit)
		'''
		SNP_position=mutant[0].strip('ACGTU-')
		SNP_positions.append(SNP_position)		
		all_mutant_partits.append(partit_adjusted_mut)

		#calculate correlation coefficient between wildtype and mutant part. func. column sums	
		raw_data = numpy.corrcoef(partit_adjusted_wt,partit_adjusted_mut)	
		corr_coeff = raw_data[0][1]
		if printPartit==False:
			out_data.append((wt_name+"/"+str(mutant[0]),corr_coeff))
		else:
			partit_adjusted_mut=[float("%.4f" % val) for val in partit_adjusted_mut]
			partit_adjusted_mut=[str(val) for val in partit_adjusted_mut]
			
			formatted_mut_partit=[float("%.4f" % val) for val in formatted_mut_partit]
			formatted_mut_partit=[str(val) for val in formatted_mut_partit]
			positionInfosMut=[str(val) for val in positionInfosMut]
			outputMaterial=str(string.join(formatted_mut_partit,","))+":"+str(string.join(positionInfosMut,","))
			out_data.append((mutant[0],outputMaterial))
		if save_vals==True:
			#save mutant matrix and basepairing prob textfiles
			if window:
				save_matrix(mt_matrix,"output/"+dirName+"/"+window+".txt")
				save_Basepairing_Probabilities("output/"+dirName+"/"+window+"_bpProbs.txt",
					mutant[1],mt_partit)
			else:
				save_matrix(mt_matrix,"output/"+dirName+"/"+mutant[0]+".txt")
				save_Basepairing_Probabilities("output/"+dirName+"/"+mutant[0]+"_bpProbs.txt",
					mutant[1],mt_partit)

	#clear out files at the end of job that aren't needed	
	os.system("rm dot.ps rna.ps")
	#print results on command line
	#print out_data
	for data in out_data:
		print ((data[0])+':'+str(data[1]))
			
	#if wt_matrix and wt_partit were overwritten, recalculate them
	if wt_partit==None:
		wt_matrix=get_matrix('wildtype',"dot.ps",wild_type,0)
		wt_partit=get_partit(wt_matrix)
		
	return out_data

def write_output_to_file(filename,headerLine,out_data):
	#creates a file with a header line and the list 'out_data', where each item is a line 
	os.system('echo "'+headerLine+'" > '+filename)
	for dataLine in out_data:
		#print dataLine
		dataLine=[str(item) for item in dataLine]
		dataLine=string.join(dataLine,' ')
		dataLine=dataLine.replace('     ','')
		os.system('echo "'+dataLine+'" >> '+filename)
	return

def check_mutant_accuracy(wild_type, SNPs):
	import re
	mutants = []
	bad_SNPs=[]

	if "NONE" in SNPs:
		return
	
	for SNP in SNPs:
		errors=False
		subSNPs=SNP.split(",")
		for subSNP in subSNPs:	
			positionset=re.findall(r"\d+",subSNP)
			if len(positionset)==1:
				position=positionset[0]
				pos=int(position)-1
				[wt_nt,mut_nt]=string.split(subSNP,position)
				if (wt_nt != wild_type[pos:pos+len(wt_nt)] and wt_nt!="-"):
					errors=True
			else:
				errors=True
			
		if errors==False:
			mutants.append(SNP)
		else:
			bad_SNPs.append(SNP)


	if len(bad_SNPs)>0:
		print "ERROR WITH THE FOLLOWING SNPs : "
		for item in bad_SNPs:
			print item
		sys.exit(2)
	else:
		pass
	return

def get_mutant_sequences(wild_type, SNPs):
	#print SNPs
	import re
	#given a wildtype sequence and a list of formatted SNPs, returns a mutant sequence for each SNP
	#Note: format = wt_nuc+position+mut_nuc (ex: A36G, U22C, G5A)
	#intialize
	mutants = []
	errors = False
	bad_SNPs = []

	for SNP in SNPs:
		#SNP=SNP.replace("-","X")
		SNPlist=string.split(SNP,",")
		mutant=wild_type
		indel_offset=0
		for subSNP in SNPlist:
						
			#if commas are found in string, then there are multiple SNPs in listed variant..
			#split comma delimited SNPs in string into a list
			#go through each sub-SNP listed in variant
			positionset=re.findall(r"\d+",subSNP)
			if len(positionset)>1:
				usage(subSNP,"improper formatting of entered SNP number")
				sys.exit(2)
			#print positionset
			
			pos_relativeToWt=int(positionset[0])-1
			pos=int(positionset[0])-1+indel_offset
			ntSet=re.findall(r"\D+",subSNP)
			
			#print ntSet
			if len(ntSet)!=2:
				usage(subSNP)
				sys.exit(2)
			[wt_nt,mut]=ntSet
			indelstatus=[0,0]
			if wt_nt=="-":
				indelstatus[0]=1
			if mut=="-":
				indelstatus[1]=1

			if sum(indelstatus)==0 and len(wt_nt)==len(mut):
				#SNP mode
				ntAtCoord=wild_type[pos_relativeToWt:pos_relativeToWt+len(wt_nt)]
				if (ntAtCoord != wt_nt and ntAtCoord != mut):
					errors=True
					bad_SNPs.append(subSNP)
					print "--------------"
					print ntAtCoord
					print wt_nt
					print mut
					print wild_type
					print "--------------"
					print "my"
					#usage(subSNP)
					#sys.exit(2)
				if mutant[pos]!="-":
					mutant=mutant[0:pos]+mut+mutant[len(mut)+pos:]
				else:
					mutant=mutant
					
			elif indelstatus[0]==1:
				mutant=mutant[0:pos]+mut.lower()+mutant[pos:]
				indel_offset+=len(mut)
				#insertion
				pass
				
			elif indelstatus[1]==1:
				mutant=mutant[0:pos]+"-"*len(wt_nt)+mutant[len(wt_nt)+pos:]
				#deletion
				pass
			else:
				usage(subSNP)
				sys.exit(2)
			
		mutants.append((SNP,mutant))
			
		#if wildtype seq at given position does not equal SNP inputed by user, error
	if errors:
		print "ERROR WITH FOLLOWING SNPS : "
		print ",".join(bad_SNPs)
		sys.exit()
	return mutants

def get_matrix(SNP,file_name,RNAseq,windowMode):
	#calculate the matrix of base-pairing probabilities for a given sequence
	os.system("echo "+RNAseq+" | RNAfold -p > /dev/null")
	#get nucCount from the length of the RNA sequence
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
	

	# These values should be the same, since the matrix is hermitian
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

def shannon_entropy(dirName, bp_matrix, label, wt):
	# Calculate the Shannon entropy values for each nucleotide while summing them
	nucCount=len(bp_matrix)
	shanValues=[0]*nucCount
	bp_probs=[0]*nucCount
	count=0
	while(count < nucCount):
		shanValue=0
		for p in bp_matrix[count]:
			if p > 0:
				shanValue = shanValue + p*math.log(p,2)
		shanValues[count]=shanValue
		bp_probs[count]=sum(bp_matrix[count])
		count += 1

	#n.b. "if wt" checks if the wt values were already calculated, i.e. selects for mutants
	if wt:
		wt_shan, wt_bp = zip(*wt)
		shan_cc = numpy.corrcoef(wt_shan, shanValues)[0,1]
		bp_cc = numpy.corrcoef(wt_bp, bp_probs)[0,1]
		cc = numpy.corrcoef(bp_probs, shanValues)[0,1]	
	
	#prepare the table to be printed
	nucleotides=[str(i)+'	' for i in range(1,len(shanValues))]
	bp_table=[str(i)+'	' for i in bp_probs]
	sh_table=[str(i)+'	' for i in shanValues]
	table=zip(nucleotides, bp_table, sh_table)
	
	if wt:
		table.append(['WT/Mutant Correlation Coefficients:\n',
		'base pairing probability: %s\n'%bp_cc,'shannon entropy: %s'%shan_cc])
		
	# output a chart containing the nucleotide number, base pairing probability and entropy	
	write_output_to_file('output/'+dirName+'/bp_prob_derived/'+label+'.txt',
		'nuc	base pairing prob	entropy (bits)',table)
		
	if label == 'wt':
		wt = zip(shanValues, bp_probs)
		return wt

def get_partit(myMatrix):
	# These values should be the same, since the matrix is hermitian
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

def interval_search_simpleAlg(allCorrCoeffValsSorted,val):
	rank=None
	if val<allCorrCoeffValsSorted[0]:
		rank=1
	elif val>allCorrCoeffValsSorted[-1]:
		rank=len(allCorrCoeffValsSorted)+1
	else:

		count=1
		while count<=len(allCorrCoeffValsSorted):
			lowerbound=allCorrCoeffValsSorted[count-1]
			upperbound=allCorrCoeffValsSorted[count]
			if val>=lowerbound and val<=upperbound:
				rank=count+1
				break
	return rank
			

def find_accurate_p_val(corrCoeffDict,muts):
	#muts is in list form
	totalNumberMuts=len(corrCoeffDict)
	
	allSingleMuts_ordered=sorted(corrCoeffDict,key=corrCoeffDict.get)	
	allSingleMutCCs_ordered=[]
	for indiv_mut in allSingleMuts_ordered:
		allSingleMutCCs_ordered.append(corrCoeffDict[indiv_mut])
	
	ranks=[]
	pvals=[]
	orphan_muts=[]
	for mut in muts:
		if mut in allSingleMuts_ordered:
			rank=allSingleMuts_ordered.index(mut)+1
			pval=float(rank)/totalNumberMuts
			pval=round(pval,5)
			pvals.append(pval)
			ranks.append(rank)
		else:
			interval_search_simpleAlg(allSingleMuts_ordered,mut)
	return (ranks,pvals)		
		

def get_accurate_p_value(RNA, mutsOfInterest, outfile_name,options):
	#get accurate p-value for each listed mutation in the RNA of interest given 
	#the user specified options, and output files if indicated by user
	[dirName,save_vals,printPartit,shannon,window]=options
	
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
	[new_wild_type,wt_matrix,wt_partit]=process_wt_seq(wild_type)
	#print len(wt_partit)
	if save_vals==True:
		#save wt matrix and basepairing prob textfiles
		save_matrix(wt_matrix,"output/"+dirName+"/wildtype.txt")
		save_Basepairing_Probabilities("output/"+dirName+"/wildtype_bpProbs.txt",wild_type,wt_partit)
	if shannon:
			shanWt = shannon_entropy(dirName, wt_matrix, 'wt', 0)
	#create various empty lists for collecting all possible single mutations in the RNA of interest
	all_mutant_partits=[]
	all_mutant_partits_of_interest=[]
	SNP_positions=[]
	SNP_positions_of_interest=[]
	nt_list='AUGC'

	allSingleMuts={}
	while (count<len(new_wild_type)):
		#do a partit func calculation for all possible single mutations in the RNA of interest
		variants=nt_list.replace(wild_type[count],'')
		origNT=wild_type[count]
		for variant in variants:
			if variant!=origNT:
				#for every possible single nt variant at every position, form mt_seq, get mt folding
				mutant=new_wild_type[0:count]+variant+new_wild_type[count+1:]
				#print mutant
				#print len(mutant)
				mt_name=str(origNT)+str(count+1)+str(variant)
				mt_matrix=get_matrix(mt_name,"dot.ps",mutant,0)
				mt_partit=get_partit(mt_matrix)
					
				if mt_name in mutsOfInterest:
					#if, as you're going through all possible mutations, you find a mutation
					#listed in the input, then if the user indicated so, save the appropriate files
					if save_vals==True:
						if window:
							save_matrix(mt_matrix,"output/"+dirName+"/"+window+".txt")
							save_Basepairing_Probabilities("output/"+dirName+"/"+window+"_bpProbs.txt",
								mutant,mt_partit)
						else:
							save_matrix(mt_matrix,"output/"+dirName+"/"+mt_name+".txt")
							save_Basepairing_Probabilities("output/"+dirName+"/"+mt_name+"_bpProbs.txt",
								mutant,mt_partit)
					if shannon:
						if window:
							shannon_entropy(dirName, mt_matrix, window, shanWt)
						else:
							shannon_entropy(dirName, mt_matrix, mt_name, shanWt)

				#calculate the correlation coefficient between the wildtype and
				#mutant sequence partition function matrix column sums
				raw_data = numpy.corrcoef(wt_partit,mt_partit)
				corr_coeff = raw_data[0][1]	
				#put the the corr coeff in an output file
				os.system("echo '"+mt_name+" "+str(corr_coeff)+"' >> "+outfile_name)
				print mt_name+" "+str(corr_coeff)
				allSingleMuts[mt_name]=corr_coeff
		count+=1


	# create empty lists to put the correlation coefficients, ranks, and p-values into
	corr_coeffs=[]
	ranks=[]
	p_values=[]
	for mutOfInterest in mutsOfInterest:
		#get correlation coefficients for mutants of interest from all corr. coeffs. calculated
		mutSeqOfInterest=get_mutant_sequences(new_wild_type,[mutOfInterest])
		print "MUT : "
		print mutSeqOfInterest
		print mutOfInterest
		print len(mutSeqOfInterest)
		print "WT : "
		print new_wild_type
		print outfile_name
		formatted_mutOfInterest=mutOfInterest[1:-1]+'|'+mutOfInterest[0]+'>'+mutOfInterest[-1]
		#will work if the mutation is a single mutation
		mutOfInterest_corr_coeff=getUnixCommandOutput(\
		'grep -w "'+mutOfInterest+'" '+outfile_name)
		print mutOfInterest_corr_coeff
		if mutOfInterest_corr_coeff==None:
			mt_matrix=get_matrix(mutOfInterest,"dot.ps",mutSeqOfInterest[0][1],0)
			mt_name=mutSeqOfInterest[0][0]
			mt_partit=get_partit(mt_matrix)
			raw_data=numpy.corrcoef(wt_partit,mt_partit)
			mutOfInterest_corr_coeff=raw_data[0][1]
			os.system("./accurate_p_value.sh "+str(mutOfInterest_corr_coeff)+" "+outfile_name+" > actual_p_value.txt")
			if save_vals==True:
				if window:
					save_matrix(mt_matrix,"output/"+dirName+"/"+window+".txt")
					save_Basepairing_Probabilities("output/"+dirName+"/"+window+"_bpProbs.txt",
						mutant,mt_partit)
				else:
					save_matrix(mt_matrix,"output/"+dirName+"/"+mt_name+".txt")
					save_Basepairing_Probabilities("output/"+dirName+"/"+mt_name+"_bpProbs.txt",
						mutant,mt_partit)
		else:
			mutOfInterest_position=mutOfInterest_corr_coeff.split(' ')[0][1:-1]
			mutOfInterest_corr_coeff=mutOfInterest_corr_coeff.split(' ')[1]
			mutOfInterest_corr_coeff=mutOfInterest_corr_coeff.replace('\n','')
			#run accurate_p_value.sh to find the ranking and p-value of the mutation of interest
			([rank],[p_value])=find_accurate_p_val(allSingleMuts,[mutOfInterest])
			#os.system("./accurate_p_value.sh "+str(mutOfInterest_corr_coeff)+" "+outfile_name+" > actual_p_value.txt")
			
			
			
		#out_file_handle = open("actual_p_value.txt","r")
		#get rank and p-value from output of accurate_p_value.sh
		#data = out_file_handle.read().split(':')
		#rank = data[0]
		#p_value = data[1]
		#p_value = p_value.replace('\n','')
		#print rank
		#print corr_coeff
		corr_coeffs.append(mutOfInterest_corr_coeff)
		ranks.append(rank)
		p_values.append(p_value)
	if dirName!=None:
		#copy corr coeff + rank + p-value to the outfile for the job at hand
		os.system('cp '+outfile_name+' output/'+dirName)
		
	#remove unneeded files at the end of the job
	#os.system('rm actual_p_value.txt '+outfile_name)
	return (corr_coeffs,ranks,p_values)

def get_bpProbs(filename):
	file=open(filename,"r")
	lines=file.readlines()
	file.close()
	bpProbs = lines[1]
	bpProbs = bpProbs.replace("\n","")
	bpProbs = bpProbs.replace("[","")
	bpProbs = bpProbs.replace("]","")
	bpProbs = bpProbs.split(", ")
	bpProbs = [float(x) for x in bpProbs]
	return bpProbs
	
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
	allChars=(string.ascii_letters+string.digits)
	randomString=''
	count=0
	while len(randomString)<stringLength:
		newChar=random.choice(allChars)
		randomString+=newChar
	return randomString
	
def direxists(dirName):
	print 'That directory already exists. Do you want to overwrite it?\n\
Everything will be deleted. Type "y" to overwrite, "n" to cancel.'
	overwrite = raw_input("> ")
	
	if overwrite in ['y','Y']:
		rmtree(os.getcwd()+'/output/'+dirName)
		print '%s was emptied.'%(os.getcwd()+'/output/'+dirName)
		os.system('mkdir output/'+dirName)
	elif overwrite in ['n','N']:
		sys.exit('Operation aborted!')
	else:
		print 'Invalid input.'
		direxists(dirName)

if __name__ == '__main__':
	#If SNPfold.py is run from the command line, commence the SNPfold program.
	#The program can not be run by importing the file in another python script.
	main(sys.argv)
	exit(0)
