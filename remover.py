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


before_path = Path('/Users/josec/Desktop/Crinoid_capture/Nov_Update/Alignments/HybTxt_gblocked')
badlist = ['PHRIX', 'STRONII', 'Endoxocrinus_macleraranus_C404j']
outdir = "/Users/josec/Desktop/Crinoid_capture/Nov_Update/Alignments/Cri_MarCut_AA_GBlo/"
for genefile in before_path.iterdir():
	if str(genefile.name).startswith('.'):
		pass
	else:
		print(genefile)
		records = [r for r in SeqIO.parse(str(genefile),"fasta") if r.id not in badlist]
		outpath = outdir + str(genefile.name)
		with open(outpath, 'w') as outhandle:
	            SeqIO.write(records, outhandle, "fasta")