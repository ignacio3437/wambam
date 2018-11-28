#!/usr/bin/python3
import os
import sys

"""
Converts the 4column output format of Baitfisher into a fasta file of bait seqs.
"""

def parser(infile):
    outlist = []
    with open(infile, 'r') as inhandle:
        lines=inhandle.readlines()
        for line in lines[1:]:
            linelist=line.split('\t')
            fasta_header = linelist[1].replace(".fasta","")
            fasta_seq = linelist[2]
            outstring = f">{fasta_header}\n{fasta_seq}\n"
            outlist.append(outstring)
    return outlist



def main():
    args = sys.argv[1:]
    usage='usage: python3 BaitFisherCleaner.py -in prefasta.txt -out baits.fasta'
    if not args:
        print(usage)
        sys.exit(1)

    if args[0] == '-in':
        infile = args[1]
        del args[0:2]
    if args[0] == '-out':
        outfile = args[1]
        del args[0:2]
    else:
        print(usage)

   
    with open(outfile,"w") as outhandle:
        for line in parser(infile):
            outhandle.write(line)
            
if __name__ == '__main__':
    main()
