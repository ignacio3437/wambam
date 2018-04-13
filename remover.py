#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import subprocess
from multiprocessing import Pool
import sys
from os.path import basename
from pathlib import Path


before_path = Path('/Users/josec/Desktop/Crinoid_capture/Nov_Update/Alignments/Cri_MarCut_AA_seqs')
badlist = ['ANNEI', 'ANTIDON', 'APORO', 'CENOL', 'CRINO', 'CYAT', 'DEMOC', 'DUMET', 'FLORO', 'GEPHY', 'ISOM', 'META', 'OLIGO', 'PROMA', 'PSATHY', 'PTILO', 'SACCO', 'STRONG', 'Spur']
outdir = "/Users/josec/Desktop/Crinoid_capture/Nov_Update/Alignments/Cri_MarCut_AA_hybseqs/"
for genefile in before_path.iterdir():
	if str(genefile.name).startswith('.'):
		pass
	else:
		print(genefile)
		records = [r for r in SeqIO.parse(str(genefile),"fasta") if r.id not in badlist]
		outpath = outdir + str(genefile.name)
		with open(outpath, 'w') as outhandle:
	            SeqIO.write(records, outhandle, "fasta")