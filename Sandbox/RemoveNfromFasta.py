#!/usr/bin/python


import os

"""
Insert description here
"""


def Nreplacer(filename):
	inputf= open(filename, 'ru')
	outfilename= 'out'+str(filename)
	output= open(outfilename, 'wu')
	lines=inputf.readlines()
	outstring=[]
	for line in lines:
		if '>' in line:
			lstr=line
		else:
			lstr=line.replace('N','-')
		outstring.append(lstr)
	for line in outstring:
		output.write(line)
	inputf.close()
	output.close()
	return 
  



def main():
	pwd=os.getcwd()
	dirs=os.listdir(pwd)
	for f in dirs:
		if '.DS' not in f and 'out' not in f:	
			Nreplacer(f)
if __name__ == '__main__':
  main()
