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
            parameters.append(pwd+line.split()[0])
    return parameters

def read_popsort(pop_order):
    with open(pop_order, "rU") as pop_order_file:
        sorted_pops=[n.strip('\n') for n in pop_order_file.readlines()]
    return sorted_pops

def name_cleaner(line):
    NAME_PATTERN = re.compile(r"""[A-Za-z0-9]+""")
    OLDNAME_PATTERN = re.compile(r"""[^\[\n\t]+""")
    namere=re.search(NAME_PATTERN, line)
    name=namere.group()
    oldnamere=re.search(OLDNAME_PATTERN, line)
    oldname= oldnamere.group()
    return [oldname,name]


def read_popfile(pop_groups):
    pop_dict={}
    with open(pop_groups, "rU") as pop_file:
        lines=[n.strip('\n') for n in pop_file.readlines()]
    header=lines[0].split(',')
    Sample_index=header.index('Sample')
    Pop_index=header.index('Pop')
    for line in lines[1:]:
        sample=line.split(',')[Sample_index]
        [oldsample,renamed_sample]=name_cleaner(sample)
        pop_dict[renamed_sample]=line.split(',')[Pop_index]
    print pop_dict
    return pop_dict

def nexuscleanuper(nexus_file, cleanup):
    #pull out the sheets in the directory and make a list of the files and taxa for each file
    nex_seqnames=[]
    NTAX_PATTERN = re.compile(r"""taxlabels""", re.IGNORECASE)
    with open(nexus_file,'rU') as f:
        lines=f.readlines()
    with open(nexus_file,'rU') as f:
        filestring=f.read()
    for i,line in enumerate(lines):
        ntax=re.search(NTAX_PATTERN,line)
        if ntax:
            while ';' not in lines[i+1]:
                nextline= lines[i+1]
                [oldname,name]=name_cleaner(nextline)
                filestring=filestring.replace(oldname,name)
                nex_seqnames.append(name)
                i+=1
    outfile_name="%s_clean"%(os.path.basename(nexus_file))
    with open(os.path.join(pwd,outfile_name),"w") as outfile:
        outfile.write(filestring)
        outfile.close()
    return nex_seqnames



def main():
    parameters=read_param(os.path.join(pwd,"parameters.txt"))
    [nexus_file,pop_groups,pop_order]=parameters
    cleanup="clean_%s"%(os.path.basename(nexus_file))
    nex_seqnames=nexuscleanuper(nexus_file, cleanup)
    #Check if pop_order file
    if "#" in pop_order:
        pop_order=False
    else:
        sorted_pops=read_popsort(pop_order)
    pop_dict=read_popfile(pop_groups)
    return


if __name__ == '__main__':
    main()
