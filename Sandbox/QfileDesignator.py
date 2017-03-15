#!/usr/bin/python


import os,sys
"""
QfileDesignator.py qfile.meanq poporder.txt
outputs qout.csv with sample and qfile group designation

"""


def readpop(popfile):
	popdict={}
	with open(popfile, 'ru') as pop:
		poplist2 = pop.readlines()
		poplist = [p.strip('\n') for p in poplist2]
		for i,pop in enumerate(poplist):
			popdict[i]=pop
	return popdict
def readq(qfile):
	groupdict = {}
	with open(qfile, 'ru') as qf:
		qlines= [q.strip('\n').split() for q in qf.readlines()]
		for q,line in enumerate(qlines):
			for i,column in enumerate(line):
				if float(column)>float(0.7):
					group=i+1
					break
				else:
					group=len(column)+1
			groupdict[q]=group
	return groupdict


def main():
	args = sys.argv[1:]
	try:
		popfile= args[1]
	except:
		print """QfileDesignator.py qfile.meanq poporder.txt
outputs qout.csv with sample and qfile group designation

			"""
		sys.exit()
	qfile= args[0]
	popdict=readpop(popfile)
	groupdict= readq(qfile)
	outfile= 'order_'+str(qfile)+'.txt'
	with open(outfile,'w') as o:
		for x in range(len(popdict)):
			o.write("%s\t%s\n"%(popdict[x],groupdict[x]))

if __name__ == '__main__':
	main()
