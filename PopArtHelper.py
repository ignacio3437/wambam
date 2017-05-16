#!/usr/bin/python
import os
import sys
import re

raw_pwd=raw_input("Please enter full path of parameter file: ")
clean_pwd=raw_pwd.replace('"','').strip()
clean_pwd=clean_pwd.replace('\\','').strip()
pwd=os.path.dirname(clean_pwd)
param_file=os.path.basename(clean_pwd)

print os.path.join(pwd,param_file)


def read_param(param_txt):
    parameters=[]
    with open(param_txt, "Ur") as param_file:
        params=param_file.readlines()
    for line in params:
        if "#" in line:
            parameters.append(os.path.join(pwd,line.split()[0]))
    return parameters

def read_popsort(pop_order):
    with open(pop_order, "Ur") as pop_order_file:
        sorted_pops=[n.strip('\n').replace(" ","_") for n in pop_order_file.readlines()]
    return sorted_pops

def name_cleaner(line):
    NAME_PATTERN = re.compile(r"[A-Za-z0-9]+")
    OLDNAME_PATTERN = re.compile(r"[^\[\n\t]+")
    namere=re.search(NAME_PATTERN, line)
    name=namere.group()
    oldnamere=re.search(OLDNAME_PATTERN, line)
    oldname= oldnamere.group()
    return [oldname,name]

def read_popfile(pop_groups):
    pop_dict={}
    with open(pop_groups, "Ur") as pop_file:
        lines=[n.strip('\n') for n in pop_file.readlines()if len(n)>1]
    header=lines[0].split(',')
    Sample_index=header.index('Sample')
    Pop_index=header.index('Pop')
    for line in lines[1:]:
        sample=line.split(',')[Sample_index]
        [oldsample,renamed_sample]=name_cleaner(sample)
        ipop=line.split(',')[Pop_index].replace(" ","_")
        pop_dict[renamed_sample]=ipop
    return pop_dict

def nexuscleanuper(nexus_file):
    #pull out the sheets in the directory and make a list of the files and taxa for each file
    nex_seqnames=[]
    NTAX_PATTERN = re.compile(r"matrix", re.IGNORECASE)
    with open(nexus_file,'Ur') as f:
        lines=f.readlines()
    with open(nexus_file,'Ur') as f:
        filestring=f.read()
    for i,line in enumerate(lines):
        ntax=re.search(NTAX_PATTERN,line)
        if ntax:
            while ';' not in lines[i+1]:
                nextline= lines[i+1]
                namepart=nextline.lstrip().split('\t')[0]
                if len(namepart)>1:
                    [oldname,name]=name_cleaner(namepart)
                    filestring=filestring.replace(oldname,name)
                    nex_seqnames.append(name)
                else:
                    pass
                i+=1
    nex_text=filestring
    return nex_seqnames,nex_text

def pop_binary(sorted_pops):
    binary_dict={}
    binary_code=[]
    for place in sorted_pops:
        binary_code.append("0")
    binary_code.pop()
    for i,place in enumerate(sorted_pops):
        binary_code.insert(i,"1")
        binary_dict[place]=",".join(binary_code)
        binary_code.pop(i)
    return binary_dict

def main():
    parameters=read_param(os.path.join(pwd,param_file))
    [nexus_file,pop_groups,pop_order]=parameters
    nex_seqnames,nex_text=nexuscleanuper(nexus_file)
    pop_dict=read_popfile(pop_groups)
    #Check if pop_order file
    if "#" in pop_order:
        sorted_pops=list(set(pop_dict.values()))
    else:
        sorted_pops=read_popsort(pop_order)
    binary_dict=pop_binary(sorted_pops)
    block1="BEGIN TRAITS;\nFormat labels=yes separator=Comma;\n"
    block2="Dimensions NTRAITS=%d;\n"%(len(sorted_pops))
    block3="TraitLabels %s;\nMatrix\n"%(' '.join(sorted_pops))
    popblock= block1+block2+block3
    nexus_file_base=os.path.basename(nexus_file.replace('.nex',''))
    with open(os.path.join(pwd,"%s_PopArt.nex"%(nexus_file_base)),'w') as outf:
        outf.write(nex_text)
        outf.write(popblock)
        for seq in nex_seqnames:
            seq_pop=pop_dict[seq]
            code= binary_dict[seq_pop]
            outf.write('%s %s\n'%(seq,code))
        outf.write(';\nEND;\n')
    return

if __name__ == '__main__':
    main()
