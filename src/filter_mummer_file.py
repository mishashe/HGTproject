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
parser.add_option( "--file","-f",dest="mumFile",action="store",
		  help="Query Species")



(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
#ori=defaultdict(list)
species1 = options.species1
species2 = options.species2
mumFile = options.mumFile


coordFileRegex = re.compile(".coo")
testFileCoord = lambda x: coordFileRegex.search(x)
#mumFileRegex = re.compile(".mum")
#testFileMum = lambda x: mumFileRegex.search(x)


myDir = "/cluster/CBIO/data1/data3/fmassip/HGT/ProjectMisha/HGTnew/data/processed/nucmer/"

fileSp1 = filter(testFileCoord,
os.listdir(myDir+species1))

fileSp2 = filter(testFileCoord,
os.listdir(myDir+species2))

#print(list(fileSp2))
AssembToRemoveSP2=dict()
for i in fileSp2:
	op = open(myDir+species2+"/"+i, 'r')
	for header in range(0,4):
		line = op.readline()
#		print("header "+line)
	for line in op.readlines():
		line = line.rstrip()
		row = line.split('\t')
		thisSp = row[7]
		if thisSp in AssembToRemoveSP2:
			AssembToRemoveSP2[thisSp]=AssembToRemoveSP2[thisSp]+"\n"+i
		else: 
			AssembToRemoveSP2[thisSp] = i  

		thisSp = row[8]
		if thisSp in AssembToRemoveSP2:
			AssembToRemoveSP2[thisSp]=AssembToRemoveSP2[thisSp]+"\n"+i
		else: 
			AssembToRemoveSP2[thisSp] = i  


AssembToRemoveSP1=dict()
for i in fileSp1:
	op = open(myDir+species1+"/"+i, 'r')
	for header in range(0,4):
		line = op.readline()
#		print("header "+line)
	for line in op.readlines():
		line = line.rstrip()
		row = line.split('\t')
		thisSp = row[7]
		if thisSp in AssembToRemoveSP1:
			AssembToRemoveSP1[thisSp]=AssembToRemoveSP1[thisSp]+"\n"+i
		else: 
			AssembToRemoveSP1[thisSp] = i  

		thisSp = row[8]
		if thisSp in AssembToRemoveSP1:
			AssembToRemoveSP1[thisSp]=AssembToRemoveSP1[thisSp]+"\n"+i
		else: 
			AssembToRemoveSP1[thisSp] = i  

#		AssembToRemoveSP1.append(row[7])
#		AssembToRemoveSP1.append(row[8])

#UniqAssembListToRemoveSP2 = list(set(AssembToRemoveSP2))
#UniqAssembListToRemoveSP1 = list(set(AssembToRemoveSP1))

f2=open("/cluster/CBIO/data1/data3/fmassip/HGT/ProjectMisha/HGTnew/HGTproject/src/test.txt","w")
for i in AssembToRemoveSP1:
	f2.write(str(i)+"\t"+str(AssembToRemoveSP1[i])+"\n")


#myDirMum = "/cluster/CBIO/data1/data3/fmassip/HGT/ProjectMisha/HGTnew/data/processed/mummer/"+species1+"_"+species2+"/"
#fileMum = filter(testFileMum,os.listdir(myDirMum))

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
		AssembSP1=line[2:i]
		f2.write(str(AssembSP1)+'\tline: '+line)
	else:
		if AssembSP1 not in AssembToRemoveSP2:
			line=re.sub("^\s+","",line)
			line=re.sub("\s+","\t",line)
			row=line.split('\t')
			line = line.rstrip()
			AssembSP2 = row[0]
			myLen = int(row[3])
			if AssembSP2 not in AssembToRemoveSP1:
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
outFile=mumFile+".h-filtered"
print(outFile)

#histo = collections.OrderedDict(sorted(histo.items()))

sorted_histo = dict(sorted(histo.items()))
#+" "+str(counter)+" / "+str(len(os.listdir(myDirMum))/2))
f1 = open(outFile,"w")
for j in sorted_histo:
	f1.write(str(histo[j])+" "+str(j)+"\n")
	
#with open (ori_file,'r') as ori: 
#	for line in ori.readlines():
#		line=line.rstrip()
#		row=line.split('\t')
#		ori.append(row)
#
