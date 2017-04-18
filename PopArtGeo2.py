#!/usr/bin/python

import os
import sys
import re
from collections import Counter
pwd="/Users/josec/Desktop/Pop2/"
#pwd=os.getcwd()

def read_param(param_txt):
    parameters=[]
    with open(param_txt, "rU") as param_file:
        params=param_file.readlines()
    for line in params:
        if "#" in line:
            parameters.append(line.split()[0])
    return parameters

def read_metadata():
    return

def read_popfile():
    return

def read_nexus():
    return

def nexuscleanuper(nexus_file, cleanup):
    #pull out the sheets in the directory and make a list of the files and taxa for each file
    NTAX_PATTERN = re.compile(r"""taxlabels""", re.IGNORECASE)
    NAME_PATTERN = re.compile(r"""[A-Za-z0-9]+""")
    OLDNAME_PATTERN = re.compile(r"""[^\[\n\t]+""")
    with open(nexus_file,'rU') as f:
        lines=f.readlines()
    with open(nexus_file,'rU') as f:
        filestring=f.read()
    for i,line in enumerate(lines):
        ntax=re.search(NTAX_PATTERN,line)
        if ntax:
            while ';' not in lines[i+1]:
                nextline= lines[i+1]
                namere=re.search(NAME_PATTERN, nextline)
                name=namere.group()
                oldnamere=re.search(OLDNAME_PATTERN, nextline)
                oldname= oldnamere.group()
                filestring=filestring.replace(oldname,name)
                i+=1
    outfile_name="%s_clean"%(os.path.basename(nexus_file))
    with open(os.path.join(pwd,outfile_name),"w") as outfile:
        outfile.write(filestring)
        outfile.close()
    print 'New nexus files written to: %s'%(cleanup)
    return

def write_out_nexus():
    return

def controller():
    nexus_file="/Users/josec/Desktop/Pop2/in.nex"
    cleanup="in_clean.nex"
    print nexuscleanuper(nexus_file, cleanup)
    return

def main():
    parameters=read_param("/Users/josec/Desktop/Pop2/param-in.txt")
    return


if __name__ == '__main__':
    main()
