#!/usr/bin/python


import os,sys
import numpy as np

"""
Problem: sumstats lists Snp stats not loci stats. Need average N per loci. 
Returns loci	average number of taxa for that loci. 
Also Prints number of Loci
"""


def read(file):
	with open(file,'rU') as f:
		data= f.readlines()
	f.closed
	return data
	
def counter(lines):
	lcount=[]
	loci=0
	taxoncovdict={}
	for l in lines[2:]:
		column=l.split('\t')
		n=int(column[8])
		if loci==0:
			loci=column[1]
			lcount=[99]
		elif loci==column[1]:
			lcount.append(n)
		else:
			taxoncov=int(np.mean(lcount))
			try:
				taxoncovdict[taxoncov]+=1
			except:
				taxoncovdict[taxoncov]=1
			loci=column[1]
			lcount=[]
			lcount.append(n)
	taxcovd=sorted(taxoncovdict.items())
			
		
	return taxcovd
	"""
	with open(os.path.join(os.path.dirname(file),'out.txt'), 'w') as o:
		o.write('Loci\#\n')
		for l,line in enumerate(data):
			vals=line.strip('\n').split('  ')
			for i,val in enumerate(vals):
				if float(val) > .9:
					o.write("%s\t%s\n"%(pops[l],int(i)+1))
	o.closed	
	return
	
"""


def main():
	file = os.path.abspath(sys.argv[1])
	lines=read(file)
	tc=counter(lines)
	for x in tc:
		print str(x[0])+'\t'+str(x[1])
	

if __name__ == '__main__':
	main()
