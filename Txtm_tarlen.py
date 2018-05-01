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


def filelines_to_list(file):
    """Makes a list of each new line in a file. """
    with open(file, 'rU') as file_handle:
        file_list = [ind.rstrip() for ind in file_handle.readlines()]
    return file_list

def txtm_tarlen(loci_list_path, org_list_path, geneseq_path, out_path):
	"""Generates a txt file that can be appended to a HypPiper output. 
	Prints the lengths of the genes for each txtm. """
	loci_list = filelines_to_list(loci_list_path)
	org_list = filelines_to_list(org_list_path)
	lens_dict = {}
	for org in org_list:
		org_lens_list = []
		for loci in loci_list:
			loci_rec_dict = SeqIO.to_dict(SeqIO.parse(f"{geneseq_path}/{loci}.fa", "fasta"))
			# loci_rec_dict = SeqIO.to_dict(SeqIO.parse(f"{geneseq_path}/TRIM_{loci}.fasta", "fasta"))

			try:
				# org_lens_list.append(len(loci_rec_dict[f"{org}"].seq))
				org_lens_list.append(len(loci_rec_dict[f"{org}-{loci}"].seq))

			except KeyError:
				org_lens_list.append('0')
			lens_dict[org] = org_lens_list
	with open(out_path, 'w') as out_handle:
		for item in lens_dict:
			line_string = ('\t').join(str(x) for x in lens_dict[item])
			out_handle.write(f"{item}\t{line_string}\n")
	return



loci_list_path = "/Users/josec/Desktop/QuickHeatmap/genelist2.txt"
org_list_path = "/Users/josec/Desktop/QuickHeatmap/txtmlist.txt"
geneseq_path = "/Users/josec/Desktop/QuickHeatmap/GeneAln_14min"
out_path = "/Users/josec/Desktop/QuickHeatmap/test2.txt"


txtm_tarlen(loci_list_path, org_list_path, geneseq_path, out_path)