import os,sys

"""
Given a fasta file. This will find duplicates and if there are any it will print a list of them.
"""

def fasdict(input):
    inputf= open(input, 'ru')
    faslines = inputf.readlines()
    inputf.close()
    fastdict={}
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
    return fastdict


def main():
    args = sys.argv[1:]
    args = ['--fas','/Users/josec/Desktop/MelosaCOI.fasta','/Users/josec/Desktop/mel2.fasta','--out', '/Users/josec/Desktop/out.csv']
    if len(args)!=5:
    	print 'usage: demultiplexerCOI.py --fas name.fa[asta] seqs.fasta --out outfile.tsv'
    	sys.exit(1)
    else:
        namefile=args[1]
        seqfile=args[2]
        outfile=args[4]
    namedict=fasdict(namefile)
    seqdict=fasdict(seqfile)
    print namedict.values()

if __name__ == '__main__':
    main()
