#!/usr/bin/python


import os,sys

"""
Removes Fasta sequences that are not in RegList

"""


def fastadicter(input):
	inputf= open(input, 'ru')
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
	# print len(fastdict)
	return fastdict

def fasparse(fastdict,reglist):
	rfile= open(reglist, 'ru')
	rlines= rfile.readlines()
	newfasdict={}
	rlist=[]
	for line in rlines:
		rlist.append(line.rstrip('\n'))
	for r in rlist:
		if fastdict.get(r):
			newfasdict[r]=fastdict[r]
		else:
			print 'Warning: %s not in Fasta file' % (r)
	rfile.close()
	# print len(newfasdict)
	return newfasdict

def main():
	args = sys.argv[1:]
	if not args:
		print 'usage: FastaExtractor --fas input.fa[asta] --list RegNumList.txt --out output.fa[asta]'
		sys.exit(1)

	if args[0] == '--fas':
		input = args[1]
		del args[0:2]
	if args[0] == '--list':
		reglist = args[1]
		del args[0:2]
	if args[0] == '--out':
		output = args[1]
		del args[0:2]
	else:
		print 'ERROR. usage: FastaExtractor --fas input.fa[asta] --list RegNumList.txt --out output.fa[asta]'



	fastdict = fastadicter(input)
	newfasdict = fasparse(fastdict,reglist)
	out=open(output,'w')
	for x in newfasdict.items():
		out.write('>'+x[0]+'\n')
		out.write(''.join(x[1]))
	out.close()




if __name__ == '__main__':
	main()
