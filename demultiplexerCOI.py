import os,sys

"""
Given a fasta file. This will find duplicates and if there are any it will print a list of them.
"""

def fasdict(input):
    inputf= open(input, 'ru')
    faslines = [line.rstrip('\n') for line in inputf.readlines()]
    inputf.close()
    fastdict={}
    for i,x in enumerate(faslines):
        if '>' in x:
            if i >0:
                if ''.join(seq) in fastdict:
                    fastdict[''.join(seq)]='dup'
                else:
                    fastdict[''.join(seq)]=name
            name = x.strip('>')
            seq=[]
        else:
            seq.append(x)
    if ''.join(seq) not in fastdict:
        fastdict[''.join(seq)]=name
    return fastdict


def merger(ndict,sdict,outfile):
    nseqs=[n for n in ndict.keys() if ndict[n] !='dup']
    with open (outfile,'w') as out:
        out.write("NEW_Name\tOLD_Name\tSeq\n")
        for n in nseqs:
            out.write('%s\t%s\t%s\n'%(ndict[n],sdict[n],n))
    return


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
    merger(namedict,seqdict,outfile)

    # print '\n'.join(namedict.keys())
    # print '\n'.join(namedict.values())

if __name__ == '__main__':
    main()
