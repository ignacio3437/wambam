#!/usr/bin/python


import re


handlenum=open('/Users/josec/Desktop/Datasets/arachnidtargets/Rix set/numbers.txt','rU')
numstring=handlenum.read()
handlenum.close()
handlepart=open('/Users/josec/Desktop/Datasets/arachnidtargets/Rix set/partitions.txt','rU')
partstring=handlepart.read()
handlepart.close()
handlegene=open('/Users/josec/Desktop/Datasets/arachnidtargets/Rix set/genespartitioned.fasta','rU')
teststring=handlegene.read()
handlegene.close()

nums=numstring.split('\n')
parts=partstring.split('\n')


for i,n in enumerate(nums):
	n= str(n)
	p=re.compile('\|%s\n'%n)
	m=('|%s\n'%parts[i])
	substring2=p.sub(m,teststring)
	teststring=substring2
	print n,m
# print nums 
# print parts
out= open('/Users/josec/Desktop/Datasets/arachnidtargets/Rix set/renamedgenes.fasta','w')
out.write(teststring)
out.close()
