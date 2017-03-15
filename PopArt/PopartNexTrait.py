#!/usr/bin/python

import os,sys,re
from collections import Counter
usage= """
usage: WAMsheetParseTrait.py --dir dir/w/tsvs --nex dir/w/nexus/gene/algn --pop popout/dir [--nexcleanup cleannex/dir]
This does the same as WAMsheetParseGeo.py but does not need lat lon data.

This script takes an input of a directory with tsv files and a directory of nexus alignments. The tsv files can have any info.
The name of the tsv must match name of nexus file.
Header line of tsvs MUST have the following columns:
 "Seq"
 "Pop"

The Seq cell in tsv must match sequence in nexus exactly. A waring will come up for sequences that don't match.
The --pop flag will create a new nexus file that is ready to load into PopArt. Sequences will be grouped by each area specified in the tsv under the  "Pop" column.

A count file will be written to the popout/dir. This file keeps a tally how many samples are in the nexus file compared to the tsv input file.

To clean up the nexus file enable the --nexcleanup option:
This will rename all of the sequences to only contain the first part of the name containing alpha-numerics and cut off everything else (like all of the stuff that Geneious appends to sequence names.)
The new cleaned nexus files will be saved to the cleannex/dir and used to make the new output files.
Example:
WAMS70039_consensussequence ----> WAMS70039
'WAMS98545(reversed)' ----> WAMS98545
'WAMS98546 ' ----> WAMS98546
"""

def grabdir(dir):
	#pull out the sheets in the directory and make a list of the files and taxa for each file
	filenames= os.listdir(dir)
	tsvs=[]
	taxalist=[]
	for file in filenames:
		if '.tsv' or '.txt' in file:
			tsvs.append(os.path.join(dir,file))
			file=file.replace('.tsv', '')
			file=file.replace('.txt', '')
			taxalist.append(file)
	return tsvs,taxalist

def tsvparser(tsv):
#returns wamdict with keys = subsamples
	f = open(tsv, 'rU')
	wamdict={}
	for line in f:
		if 'Seq' in line:
			#find index of important things in sheet
			headerlist= line.rstrip('\r\n').split('\t')
			wami=headerlist.index('Seq')
			areai= headerlist.index('Pop')
		try:
			if 'Seq' in line:
				continue
			#declare everything
			ll= line.rstrip('\r\n').split('\t')
			sarea=ll[areai]
			area=sarea.replace(' ','')
			#wamdict construction
			ldict={}
			ldict['area']= area
			wamdict[ll[wami]]=ldict
		except:
			pass
	return wamdict

def seqlistmaker(nex):
	#return the names of the sequences in nexus file
	f=open(nex,'rU')
	seqlist=[]
	NTAX_PATTERN = re.compile(r"""taxlabels""", re.IGNORECASE)
	NAME_PATTERN = re.compile(r"""[^\[\n\t]+""")
	lines=f.readlines()
	for i,l in enumerate(lines):
		ntax=re.search(NTAX_PATTERN,l)
		if ntax:
			while ';' not in lines[i+1]:
				nextline= lines[i+1]
				namere=re.search(NAME_PATTERN, nextline)
				name=namere.group()
				seqlist.append(name)
				i+=1
	return seqlist

def nexuscleanuper(nexfiles, cleanup):
	#pull out the sheets in the directory and make a list of the files and taxa for each file
	NTAX_PATTERN = re.compile(r"""taxlabels""", re.IGNORECASE)
	NAME_PATTERN = re.compile(r"""[A-Za-z0-9]+""")
	OLDNAME_PATTERN = re.compile(r"""[^\[\n\t]+""")
	for nex in nexfiles:
		f=open(nex,'rU')
		filestring=f.read()
		f=open(nex,'rU')
		lines=f.readlines()
		f.close()
		for i,l in enumerate(lines):
			ntax=re.search(NTAX_PATTERN,l)
			if ntax:
				while ';' not in lines[i+1]:
					nextline= lines[i+1]
					namere=re.search(NAME_PATTERN, nextline)
					name=namere.group()
					oldnamere=re.search(OLDNAME_PATTERN, nextline)
					oldname= oldnamere.group()
					filestring=filestring.replace(oldname,name)
					i+=1
		outfile=open(os.path.join(cleanup,os.path.basename(nex)),'w')
		outfile.write(filestring)
		outfile.close()
		newnexfiles=[os.path.join(cleanup,x) for x in os.listdir(cleanup) if '.nex' in x]
	print 'New nexus files written to: %s'%(cleanup)
	return newnexfiles

def subseqparser(subsamples,nex):
	#return list of subsamples that have sequence info
	rsubsam = []
	seqlist = seqlistmaker(nex)
	for sub in subsamples:
		if sub in seqlist:
			rsubsam.append(sub)
	for seq in seqlist:
		if seq not in rsubsam:
			print 'WARNING: %s is in %s but not in tsv sheet.'%(seq,os.path.basename(nex))
	return rsubsam

def nexfinder(nex):
	nexlist = nex.split('/')
	for x in nexlist:
		if '.nex' in x:
			return x.replace('.nex','')

def counter(rsubsamples,nsubsamples,wamdict,nexname, outc):
	locs=[]
	seqeddict={}
	notseqdict={}
	for rsub in rsubsamples:
		locs.append(wamdict[rsub]['area'])
	seqeddict= Counter(locs)
	seqlocations=seqeddict.keys()
	locations= locs
	locs=[]
	for nsub in nsubsamples:
		locs.append(wamdict[nsub]['area'])
		locations.append(wamdict[nsub]['area'])
	notseqdict= Counter(locs)
	locations= list(set(locations))
	for l in locations:
		outc.write(nexname+'\t'+l+'\t'+str(seqeddict[l])+'\t'+str(notseqdict[l])+'\n')
	return seqlocations

def traitmaker(seqlocations):
	numlocs= len(seqlocations)
	traitkey={}
	for i,sl in enumerate(seqlocations):
		olist=['0' for x in range(len(seqlocations))]
		olist.pop()
		olist.insert(i,'1')
		traitkey[sl]=olist
	return traitkey

def main():
	args = sys.argv[1:]
	if not args:
		print usage
		sys.exit(1)

	if args[0] == '--dir':
		dir = args[1]
		del args[0:2]
	if args[0] == '--nex':
		nexdir= args[1]
		nexfiles=[os.path.join(args[1],x) for x in os.listdir(args[1]) if '.nex' in x]
		del args[0:2]
	if '--nexcleanup' in args:
		try:
			cleanup = args[args.index("--nexcleanup")+1]
		except:
			print 'Please specify directory for cleaned nexus file.'
			sys.exit(1)
		cleanup = args[args.index("--nexcleanup")+1]
		if not os.path.exists(cleanup):
			os.mkdir(cleanup)
		nexfiles = nexuscleanuper(nexfiles, cleanup)
	else:
		cleanup=False
	tsvs,taxalist= grabdir(dir)
	if args[0] == '--pop':
		popout = os.path.abspath(args[1])
		print 'POP '+popout
		if not os.path.exists(popout):
			os.mkdir(popout)
		del args[0:2]
		#Make header(popblock) for nexus file:
		block1="""
BEGIN TRAITS;
Format labels=yes separator=Comma;
"""

		#get data and print to outfile
		countoutfile= popout+'/'+'counts.txt'
		outc= open(countoutfile, 'w')
		outc.write('Taxa'+'\t'+'Area'+'\t'+'Seqed'+'\t'+'notSeqed'+'\n')
		for nex in nexfiles:
			if '.nex' in nex:
				nexname= nexfinder(nex)
				tsvtaxa=tsvs[taxalist.index(nexname)]
				wamdict = tsvparser(tsvtaxa)
				subsamples=wamdict.keys()
				rsubsamples=subseqparser(subsamples,nex)
				nsubsamples=[]
				for sub in subsamples:
					if sub not in  rsubsamples:
						nsubsamples.append(sub)
				seqlocations = counter(rsubsamples,nsubsamples,wamdict,nexname,outc)
				traitkey=traitmaker(seqlocations)
				popoutfile= popout+'/'+os.path.basename(nex).rstrip('nex')
				block2="Dimensions NTRAITS=%d;\n"%(len(traitkey.keys()))
				block3="TraitLabels %s;\nMatrix\n"%(' '.join(seqlocations))
				popblock= block1+block2+block3
				outf= open(popoutfile+'nex', 'w')
				inf = open(nex,'rU')
				inlines= inf.read()
				outf.write(inlines)
				outf.write(popblock)
				for r in rsubsamples:
					ll=traitkey[wamdict[r]['area']]
					outf.write('%s %s\n'%(r, ','.join(ll)))
				outf.write(';\nEND;\n')
				outf.close()
if __name__ == '__main__':
	main()
