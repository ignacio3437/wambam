#!/usr/bin/python
import re
usage="""
Usage: VulgarityFilter.py --in exonerate_outfile.txt

Parses the output of exonerate wth --ryo '%qi\t%pi\t%qas\t%V\tEND\n'
output is the inputname.fa of a fasta with each exon listed as a seperate sequence
example exonerate command: exonerate --model est2genome Test44.fasta /Users/josec/Desktop/NudiSilicoTest/Exonerate/acl_ref_AplCal3.0_chrUn.fa -Q DNA -T DNA --showvulgar F --showalignment F --percent 90 --verbose 0 --ryo '%qi\t%pi\t%qas\t%V\tEND\n' --fsmmemory 20G --bestn 1 > exonerate_outfile.txt

"""
def openfile(infile):
    with open(infile,'Ur') as inhandle:
        target_list=inhandle.read().split('END')
    return target_list

def parser(target):
    target_dict={}
    #make list of items to be put in dictionary and remove new lines from strings
    tlist=[x.replace('\n','') for x in target.strip('\n').split('\t')]
    #name the parts to be put in dict from tlist
    target_id=tlist[0]
    pecentscore=tlist[1]
    CDS=tlist[2]
    vulgar_raw=tlist[3]
    #extract exon lengths from vulgar output.
    vlist=re.findall(r'M\s\d*',vulgar_raw)
    vulgar=[int(x.replace('M ','')) for x in vlist]
    #make dictionary of all info to recall later
    target_dict['Target']=target_id
    target_dict['Percent']=pecentscore
    target_dict['CDS']=CDS
    target_dict['Vulgar']=vulgar
    return target_dict

def splitter(string,vlist):
    #cuts up the CDS into exon sequences and stores them in a sequential list
    counter=0
    exonseq_list=[]
    for v in vlist:
        exonseq_list.append(string[counter:counter+v])
        counter+=v
    return exonseq_list

def writer(target_dict):
    #writes the exon sequences into the outfile.fa
    outstring=""
    vlist=target_dict['Vulgar']
    CDS=target_dict['CDS']
    exonseq_list=splitter(CDS,vlist)
    for exon_number,exonseq in enumerate(exonseq_list):
        outstring+=(">%s_%d\n%s\n")%(target_dict['Target'],exon_number,exonseq)
    return outstring

def main():
    #args = sys.argv[1:]
    args=['--in','/Users/josec/Desktop/NudiSilicoTest/Exonerate/Test44/test444out.txt']
    if not args:
        print usage
        sys.exit(1)
    if args[0] == '--in':
        infile = args[1]
        del args[0:2]
    outfile_name=infile.replace('.txt','.fa')
    outfile=open(outfile_name,'a')
    target_list= openfile(infile)
    for target in target_list:
        #check if target is blank
        if len(target)>2:
            target_dict=parser(target)
            ####Filter to remove low ID hits####
            if float(target_dict['Percent'])<95:
                continue
            outstring=writer(target_dict)
            outfile.write(outstring)
    outfile.close()

if __name__ == '__main__':
    main()
