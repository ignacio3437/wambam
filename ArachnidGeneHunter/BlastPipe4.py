#!/usr/bin/python


from Bio import SeqIO
import sys,os, commands
import shutil
import re

from collections import defaultdict

"""
#Parameter is something to change stringency of criteria
currently need more than 20 of 27 taxa to have hits
and length of hit must be more than half of query
and evalue is 1e-20
"""


def blastsseq(query,db,md):
	#Parameter:
	blastcommand = "blastn -query %s -task 'dc-megablast' -db %s -evalue '1e-20' -outfmt '6 qseqid sseqid qlen length sseq' -max_target_seqs 1"%(query,db)
	print blastcommand
	#returns dictionary with key = target and value = list of blast hits in transcriptome
	output= commands.getoutput(blastcommand)
	outlines= output.split('\n')
	hits=[]
	md2={}
	if '?' in outlines:
		print 'remove "?" from fastaqueryfile'
		sys.exit()
	for line in outlines:
		entrylist = line.split('\t')
		hits.append(entrylist)
	taxa= db.split('.')[0]
	taxaseq=taxa+'s'
	if len(hits[0])>1:
		for hit in hits:
			placeholder={}
			length= int(hit[3])
			qlen = int(hit[2])
			#Parameter:
			tooshort= qlen/2	
			if length>tooshort:
				placeholder[taxa]=hit[1]
				placeholder[taxaseq]=hit[4]
				md2[hit[0]]=placeholder
	for key in md.keys():
		if key in md2:
			md2[key].update(md[key])
		else:
			md2[key]=md[key]
	return md2
	
	
def fastafish(longstring,taxa):
	filestring= longstring
	myre=r">(\w*)"
	mysub=r">%s_\1"%taxa
	newstring2=re.sub(myre,mysub,filestring)
	myre2=r">(\w*)( .*)\n"
	mysub2=r">\1\n"
	newstring=re.sub(myre2,mysub2,newstring2)
	return newstring

def mafftalign(seqdir,aligndir):
	os.chdir(seqdir)
	toalign = [f for f in os.listdir(seqdir) if f.endswith('.fasta')]
	for fas in toalign:
		cmnd =  'mafft --auto %s > %s'%(fas, aligndir+'/'+fas)
		print cmnd
		output = commands.getoutput(cmnd)
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
		try:
			md
		except:
			md={}
		md = blastsseq(targetfas,db,md)
	#md is dict of with key for each target and then a sub dict of target[database]=seqnameofhit
	taxas=[db.split('.')[0] for db in dbs]
	targets = md.keys()
	writecount=0
	tooshortcount=0
	for target in targets:
		#Parameter
		if len(md[target].keys())>19:
			filename = os.path.join(tarout,'%s.fasta'%target)
			print 'writing to%s'%filename
			stringl=[]
			for taxa in taxas:
				staxa=taxa+'s'
				keyerror=0
				writecount2=0
				try:
					string = '>%s\n%s\n'%(taxa,md[target][staxa])
					stringl.append(string)
					writecount2+=1
				except KeyError:
					keyerror+=1
			print target
			print md[target].keys()
			print writecount2
			print keyerror
			longstring=''.join(stringl)
			taroutf=open(filename,'w')
			taroutf.write(longstring)
			taroutf.close()
			writecount+=1
		else:
			tooshortcount+=1
	print '%d files written, %s targets had too few hits' %(writecount,tooshortcount)
	try:
		mafftalign(tarout,mafftout)
	except:
		print 'nope'
		pass
	print '\a'
if __name__ == '__main__':
	main()
