#!/usr/bin/python


import os,sys,re,collections

"""
Selects Primers from Geneious annotation output.csv

Each alignemnt must have sequences with the following name format:
Something_Something_Gene (ex: FamilyCode_species_gene, "CAACte_Cteniz_ABCF2")
Then create new primer (1 pair) in geneius> Save them to a folder> select all alignments and test each sequence to the primer pairs in folder> extract PCR products> extract annotations from PCR products
You can allow for degeneracy, and mismatches. Products should be 400-1200 bp. Should turn off mistmatches within 2 bp of 3' end

"""


def importer(filename):
#returns Annotation Dictionary by parsing csv
	ADict={}
	with open(filename,'rU') as incsv:
		linesh=incsv.readlines()
	#remove newline and spaces
	linesh2=[l.strip('\n').replace(' ','_') for l in linesh]
	#remove numerical commas ex: "1,200,345"
	lines=[re.sub(r'"(.*?)\,(.*?)"',r'\1\2',l) for l in linesh2]
	"""	
	with open('/Users/josec/Desktop/test.txt', 'w') as outf:
		for x in lines:
			outf.write(x+'\n')
	"""	
	header_line=lines.pop(0)
	headernames=header_line.split(',')
	position={}
	for i,x in enumerate(headernames):
		val=x
		position[val]=i
	testfor1=['Sequence_Name','Name','Sequence','Minimum','Maximum','Product_Size','Direction','Degeneracy','Mismatches']
# 	make sure correct attributes were exported from geneious
	for x in testfor1:
		try:
			position[x]
		except:
			print 'Make sure the following annoitations are exported from geneious: Sequence Name,Name,Sequence,Minimum,Maximum,Product Size,Direction,Degeneracy,Mismatches'
			sys.exit(1)
#Parses info in csv to ADict file
	for i,line in enumerate(lines):
		linedict={}
		for key in testfor1:
			loc=position[key]
			value=line.split(',')[loc]
			linedict[key]=value
		ADict[str(i)]=linedict
	return ADict

def findPairs(ADict,filename):
#Looks for the forward and reverse primer of each PCR product and stores info
	Fdict={}
	Rdict={}
	misprimed=[]
	for key in ADict.keys():
		if ADict[key]['Name'].split('_')[0]==ADict[key]['Sequence_Name'].split('_')[2]:
				PS= int(ADict[key]['Product_Size'])
				if PS > 399 and PS < 1200:
					if ADict[key]['Minimum']=='1':
						Fdict[ADict[key]['Sequence_Name']]=ADict[key]
					elif ADict[key]['Maximum']==ADict[key]['Product_Size']:
						Rdict[ADict[key]['Sequence_Name']]=ADict[key]
		else:
		# Will store misprimed info in misprimed list			
			misprimed.append(ADict[key]['Name'])
#Only keep pairs in dictionary.
	for key in Fdict.keys():
		try:
			Rdict[key]
		except:
			del Fdict[key]
	for key in Rdict.keys():
		try:
			Fdict[key]
		except:
			del Rdict[key]			
	PCRP= Fdict.keys()
	#PCRP is a list of PCR products
# Count how many time pair bound to a target in an alignment and how many times it misprimed. Store in Mispdict and Corpdict
	Mispdict={}
	Corpdict={}
	Mismdict={}
	cmisprimed=collections.Counter(misprimed)
	alignments=[]
	for pcrp in PCRP:
		pair=Fdict[pcrp]['Name'].split('_')[0]
		alignments.append(pair)
		if Fdict[pcrp]['Name'] in misprimed:
			Mispdict[pair]=cmisprimed[Fdict[pcrp]['Name']]
		else:
			Mispdict[pair]=0	
	calignments=collections.Counter(alignments)
	for pcrp in PCRP:
		pair=Fdict[pcrp]['Name'].split('_')[0]
		Corpdict[pair]=calignments[pair]	
#Sum mismatches 
		mismatches = int(Fdict[pcrp]['Mismatches'])+ int(Rdict[pcrp]['Mismatches'])
		try:
			Mismdict[pair]=mismatches+Mismdict[pair]
		except:
			Mismdict[pair]=mismatches
	Degdict={}
	Ampldict={}
	Seqdict={}
	genes=Mispdict.keys()
	for pcrp in PCRP:
		pair=Fdict[pcrp]['Name'].split('_')[0]	
		#Get sum of degeneracies in primer pair and add to Degdict
		try:
			Degdict[pair]
			continue
		except:
			if Fdict[pcrp]['Degeneracy'] or Rdict[pcrp]['Degeneracy']:
				try:
					Degdict[pair]=0
					Degdict[pair]=int(Fdict[pcrp]['Degeneracy'])+Degdict[pair]
					Degdict[pair]=int(Rdict[pcrp]['Degeneracy'])+Degdict[pair]
				except:
					pass
			else:
				Degdict[pair]=0
		#Get PCR product length and add to Ampldict
		try:
			Ampldict[pair]
			continue
		except:
			Ampldict[pair]=int(Fdict[pcrp]['Product_Size'])
		#Get primer forward and reverse sequence and add to Seqdict as tuple
		try:
			Seqdict[pair]
			continue
		except:
			Seqdict[pair]=(Fdict[pcrp]['Sequence'],Rdict[pcrp]['Sequence'])
			continue
	outdir=os.path.dirname(filename)
	outname="Primers."+os.path.basename(filename)
	outpath=os.path.join(outdir,outname)
	with open(outpath,'w') as outf:
		outf.write("Gene,Matches,Misprimes,Mismatches,DegenerateBP,Amplength,Fseq,Rseq\n")		
		for g in genes:
			glineh= [g,Corpdict[g],Mispdict[g],Mismdict[g],Degdict[g],Ampldict[g],Seqdict[g][0],Seqdict[g][1]]
			gline=[str(x) for x in glineh]
			outf.write(','.join(gline)+'\n')
	return	"%d pairs printed to %s"%(len(genes),outpath)	


	
	
	
	
def main():
	args = sys.argv[1:]
	if not args:
		print 'usage: PrimerSelector.py --infile geneiousannnotations.csv '
		sys.exit(1)
	if args[0] == '--infile':
		filename = args[1]
		del args[0:2]
	ADict=importer(filename)
	print findPairs(ADict,filename)
if __name__ == '__main__':
	main()
