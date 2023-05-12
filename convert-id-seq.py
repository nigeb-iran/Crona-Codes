#from Bio import SeqIO
import string
import sys
import numpy
from pathlib import Path
import re
iddd=[]
handle = open("metadata2","r")
handle3 = open("sequence.fa","r")
seqq=[]
seq=[]
for lin in handle3 :
  while (lin[0]!=">"):
      seq.append(lin)
      break
for i in range(0, len(seq)):
	a=handle.readline()
	j=a.split("/")[5][0:].strip()
	name=j+".fas"
	out=open(name,"a")
	out.write(">")
	out.write(a)
	#out.write("\n")
	out.write(seq[i])
	#out.write("\n")
	


