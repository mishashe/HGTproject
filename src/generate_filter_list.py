# -*-coding:utf8 -*
import os
import sys
import pickle
#import numpy as np
import random
import re
from collections import defaultdict
from optparse import OptionParser
import itertools
################################



parser = OptionParser()
parser.add_option( "--species1",dest="species1",action="store",
		  help="Ref Species")
parser.add_option( "--country1",dest="country1",action="store",
		  help="Species country")
parser.add_option( "--threshold","-t",dest="threshold",action="store",default=0.002,
		  help="Distance Threshold below which files should be filtered out")




(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
#ori=defaultdict(list)
species1 = options.species1
country1 = options.country1

#mumFile = options.mumFile
threshold = float(options.threshold)
myDir = "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/dist/"

def AssembToRemove(species,country,threshold=threshold):

	filterFile="/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/processed/toFilterMash/"+species+"/"+species+"-"+country+"/toFilter.txt"

	AssembToRemove=[]
	comm="mkdir -p /cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/processed/toFilterMash/"+species+"/"+species+"-"+country
	os.system(comm)
	fileRegex = re.compile(species+"-"+country)
	testFileSp = lambda x: fileRegex.search(x)
	myFiles = filter(testFileSp,
	os.listdir(myDir+species))


	for i in myFiles:
		print(i)
		print(species+"-"+country+"-"+country)
		if i == species+"-"+country+"-"+country :
			print("do not use this file to filter")
		else:
			myPath=myDir+species+"/"+i+"/dist-close.txt"
			if os.path.exists(myPath):
				op = open(myPath, 'r')
				thisCount = 0
				for line in op.readlines():
					line = line.rstrip()
					row = line.split('\t')
					thisAssemb1 = row[0]
					thisAssemb2 = row[1]
					dist = float(row[2])

					if dist < threshold : 
						thisCount = thisCount + 1
						thisAssemb1=re.sub(r".+\/(.*).fasta",r"\1",thisAssemb1)
						thisAssemb2=re.sub(r".+\/(.*).fasta",r"\1",thisAssemb2)
						AssembToRemove.append(thisAssemb1)
	#					AssembToRemove.append(thisAssemb2)
				print("count: "+str(thisCount))
		
	AssembToRemove=set(AssembToRemove)
	op=open(filterFile,"w")
	op.write("\n".join(i for i in AssembToRemove)+"\n")
	op.close()

	return(AssembToRemove)

AssembToRemoveSP1=AssembToRemove(species1,country1,threshold=threshold)





