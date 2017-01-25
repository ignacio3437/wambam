import os,sys

"""
Given a fasta file. This will find duplicates and if there are any it will print a list of them.
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
			nel=''
		else:
			nel+=line.rstrip('\n')
		if i ==len(faslines)-1:
			fastdict[fasname]=nel
		fastdict[fasname]=nel
	inputf.close()
	return fastdict


def main():
	args = sys.argv[1:]
	if not args:
		print 'usage: BarcoderCOI.py --fas input.fa[asta] '
		sys.exit(1)
	else:
		input=args[1]
	if '.fas' in input:
		fastdict=fastadicter(input)
		countmult(fastdict)
	else:
		print 'Please use fasta file.'
		sys.exit(1)
if __name__ == '__main__':
	main()
