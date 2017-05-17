#!/usr/bin/python
import os
import sys

"""
Makes a new fasta file of each group designated in the input csv file from a fasta file.
Name in sheet must match fasta file exatcly. No headers. First column will be name of new fasta file of all seqs specified if second column.

Eg input file will create 1.fasta and 2.fasta:
1,denovoMyti50draft_rep_c46044
1,denovoMyti50draft_dn_rep_c46045
1,denovoMyti50draft_rep_c46046
2,denovoMyti50draft_rep_c46047
2,denovoMyti50draft_rep_c46
"""

def fastadicter(input):
	inputf= open(input, 'Ur')
	fastdict={}
	faslines = inputf.readlines()
	for i,line in enumerate(faslines):
		if i ==0:
			nel=0
		if '>' in line:
			if nel:
				fastdict[fasname]=nel
			fasname = line.strip('>').rstrip('\n')
			nel=[]
		else:
			nel.append(line)
		if i ==len(faslines)-1:
			fastdict[fasname]=nel
		fastdict[fasname]=nel
	inputf.close()
	return fastdict



def inputdict(reglist):
	handle = open(reglist,'Ur')
	regdict={}
	relines=handle.readlines()
	handle.close()
	for line in relines:
		if ',' in line:
			line=line.strip('\n')
			linelist=line.split(',')
			regdict[linelist[1]]=linelist[0]
	return regdict

def fasparse(fastdict,reglist):
	newfasdict={}
	for r in reglist:
		if fastdict.get(r):
			newfasdict[r]=fastdict[r]
		else:
			print 'Warning: %s not in Fasta file' % (r)
	return newfasdict




def main():
	args = sys.argv[1:]
	usage='usage: FastaExtractor --fas input.fa[asta] --list RegNumList.csv'
	if not args:
		print usage
		sys.exit(1)

	if args[0] == '--fas':
		input = args[1]
		del args[0:2]
	if args[0] == '--list':
		reglist = args[1]
		del args[0:2]
	else:
		print usage

	fastdict=fastadicter(input)
	regdict=inputdict(reglist)
	outs= set(regdict.values())
	for out in outs:
		subreglist=[seq for seq in regdict.keys() if regdict[seq] == out]
		newfasdict=fasparse(fastdict,subreglist)
		outfilename=str('%s.fasta')%out
		outhandle=open(outfilename,'w')
		for x in newfasdict.items():
			outhandle.write('>'+x[0]+'\n')
			outhandle.write(''.join(x[1]))
		outhandle.close()



if __name__ == '__main__':
	main()
