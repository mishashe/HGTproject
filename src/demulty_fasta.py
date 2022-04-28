# -*-coding:utf8 -*

import os
import sys
import pickle
#import numpy as np
import random
import re
from collections import defaultdict
from optparse import OptionParser
################################



parser = OptionParser()
parser.add_option("-o","--out", dest="out_dir",action="store",
		  help="directory to put output files")
parser.add_option("-i", "--input",dest="input_file",action="store",
		  help="input multifasta file")


(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
out_dir=options.out_dir
input_file=options.input_file

#ori=defaultdict(list)

listFile=[]

f1=open(out_dir+"/control_file.txt","w")
with open (input_file,'r') as op: 
	for line in op.readlines():
		if line[0] == ">":
			
			f1.write("\n")
			f1.close()
			i=1 
			while i < len(line) and line[i] != ' ':
				i=i+1
			Name=line[1:i]
			f1=open(out_dir+"/"+Name+".fasta","w")
			f1.write(line)
			listFile.append(Name)
		else:
			line=line.rstrip()
			f1.write(line)

f1.write("\n")
f1.close()

for myFile in listFile:
	comm='mash sketch -s 10000 '+out_dir+"/"+myFile+".fasta"
	os.system(comm)
	



