#!/usr/bin/python


from Bio import SeqIO
import sys,os, commands
import shutil
import re

from collections import defaultdict

"""
 
"""


def blastsseq(query,db):
	#Parameter:
	blastcommand = "tblastn -query %s -task 'tblastn' -num_threads 7 -db %s -outfmt '6 qseqid sseqid sstart send' -max_target_seqs 1 -evalue 1e-70"%(query,db)
	print blastcommand
	output= commands.getoutput(blastcommand)
	outlines= output.split('\n')
	hits=[]
	md={}
	if '?' in outlines:
		print 'remove "?" from fastaqueryfile'
		sys.exit()
	for line in outlines:
		entrylist = line.split('\t')
		hits.append(entrylist)
	taxa= db.split('.')[0]
	taxaseq=taxa+'s'
	print hits
	if len(hits[0])>1:
		for hit in hits:
			if float(hit[2])<float(hit[3]):
				start=hit[2]
				end=hit[3]
			else:
				start=hit[3]
				end=hit[2]
			blastdbcommand= "blastdbcmd -db %s -entry \"%s\" -range %s-%s -outfmt %%s"%(db,hit[1],start,end)
			dbout= commands.getoutput(blastdbcommand)
			seq=dbout
			fasstring= '>'+taxa+'_'+hit[0]+'\n'+seq+'\n'
			md[hit[0]]=fasstring
	else:
		pass
		print md
	return md


def mafftalign(seqdir,aligndir):
	os.chdir(seqdir)
	toalign = [f for f in os.listdir(seqdir) if f.endswith('.fasta')]
	for fas in toalign:
		cmnd =  'mafft --maxiterate 1000 --localpair --adjustdirectionaccurately --reorder --quiet --thread 7 %s > %s'%(fas, aligndir+'/'+fas)
		print cmnd
		output = commands.getoutput(cmnd)
		print output
	return

def main():
	usage="BlastPipe.py --fas fastaoftargetstoblast --blast dir/blastdb --out dir/out --db dir/database/ --mafft dir/mafftout"
	args = sys.argv[1:]
	if len(args)==0:
		print usage
		sys.exit()
	if args[0] == '--fas':
		targetfas =os.path.abspath(args[1])
		del args[0:2]
	if args[0]=='--blast':
		blastdir= os.path.abspath(args[1])
		del args[0:2]
	if args[0]=='--out':
		tarout= os.path.abspath(args[1])
		if os.path.exists(tarout):
			shutil.rmtree(tarout)
			os.mkdir(tarout)
		else:
			os.mkdir(tarout)
		del args[0:2]
	if args[0]=='--db':
		os.chdir(args[1])
		fastadbs= [os.path.abspath(fas) for fas in os.listdir(args[1]) if fas.endswith('.fasta')]
		del args[0:2]
	if args[0] =='--mafft':
		mafftout=os.path.abspath(args[1])
		if os.path.exists(mafftout):
			shutil.rmtree(mafftout)
			os.mkdir(mafftout)
		else:
			os.mkdir(mafftout)
		del args[0:2]
	dbs = [n.replace('.nin','') for n in os.listdir(blastdir) if n.endswith('.nin')]
	os.chdir(blastdir)
	for db in dbs:
		md = blastsseq(targetfas,db)
		for x in md.keys():
			with open('%s/%s.fasta'%(tarout,x),'a+') as ahandle:
				ahandle.write(md[x])
	mafftalign(tarout,mafftout)
if __name__ == '__main__':
	main()
