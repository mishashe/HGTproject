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
parser.add_option( "--species2",dest="species2",action="store",
		  help="Query Species")
parser.add_option( "--country1",dest="country1",action="store",
		  help="country Ref Species")
parser.add_option( "--country2",dest="country2",action="store",
		  help="country Query Species")
parser.add_option( "--file","-f",dest="mumFile",action="store",
		  help="Query Species")


(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
#ori=defaultdict(list)
species1 = options.species1
species2 = options.species2
country1 = options.country1
country2 = options.country2

mumFile = options.mumFile

#coordFileRegex = re.compile(".coo")
#testFileCoord = lambda x: coordFileRegex.search(x)

#mumFileRegex = re.compile(".mum")
#testFileMum = lambda x: mumFileRegex.search(x)


#myDir = "/cluster/CBIO/data1/data3/fmassip/HGT/ProjectMisha/HGTnew/data/processed/nucmer/"
myDir = "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/dist/"

def AssembToRemove(species,country):

	filterFile="/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/data/processed/toFilterMash/"+species+"/"+species+"-"+country+"/toFilter.txt"

	AssembToRemove=[]
	op=open(filterFile,"r")
	for line in op.readlines():
		line=line.rstrip()
		AssembToRemove.append(line)
	return(AssembToRemove)



AssembToRemoveSP1=AssembToRemove(species1,country1)
AssembToRemoveSP2=AssembToRemove(species2,country2)


print(species1+" "+country1+"\t".join(i for i in AssembToRemoveSP1))
print(species2+" "+country2+"\t".join(i for i in AssembToRemoveSP2))

histo = dict()
op = open(mumFile, 'r')
countRMSP1=0
countRMSP2=0

f2=open("AssembSp1.txt","w")
for line in op.readlines():
	if line[0]== ">":
		i=2
		while i < len(line) and line[i] != ' ':
			i=i+1
		AssembSP2=line[2:i]
		f2.write(str(AssembSP2)+'\tline: '+line)
	else:
		if AssembSP2 not in AssembToRemoveSP2:
			line=re.sub("^\s+","",line)
			line=re.sub("\s+","\t",line)
			row=line.split('\t')
			line = line.rstrip()
			AssembSP1 = row[0]
			myLen = int(row[3])
			if AssembSP1 not in AssembToRemoveSP1:
				if myLen in histo:
					histo[myLen] = histo[myLen]+1
				else:
					histo[myLen] = 1
			else:
				countRMSP2 = countRMSP2 + 1 
		else:
			countRMSP1 = countRMSP1 + 1

print("filtered SP1 : "+str(countRMSP1))
print("filtered SP2 : "+str(countRMSP2))

op.close()
outFile=mumFile+".h-filtered-mash"
print(outFile)

#histo = collections.OrderedDict(sorted(histo.items()))

sorted_histo = dict(sorted(histo.items()))
#+" "+str(counter)+" / "+str(len(os.listdir(myDirMum))/2))
f1 = open(outFile,"w")
for j in sorted_histo:
	f1.write(str(histo[j])+" "+str(j)+"\n")
	
