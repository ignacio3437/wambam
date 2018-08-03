#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
import re
import sys

def partitioner(part_file):
	#Parse Raxml style partition file
	#the partiotions list returned the start and stop coordinates of every gene in the supermatrix (concatenaged gene alignments).
	partition_strings = []
	partitions = []
	with open(part_file) as f:
		for line in f:
			line=line.strip ()
			line = line.split (" = ")
			partition_strings.append (line [1])
	for part in partition_strings:
		m = re.search("(.*)-(.*)",part)
		if m:
			start = int(m.groups()[0])
			end = int(m.groups()[1])
		partitions.append((start,end))
	return partitions

def alncutter(partitions,aln_file,genealn_outdir):
	#cuts matrix to create new alignment for each gene 
	#writes gene alignments to genealn_outdirectory
	for start,end in partitions:
		out_path = genealn_outdir/f"{start}.fas"
		with open(out_path,'w') as out:
			for record in SeqIO.parse (aln_file.name, "phylip-relaxed"):
				sequence = record.seq[start:end]
				#Skip seqs made only of gaps
				if len(set(sequence)) > 1: 
					out.write(f">{record.id}\n")
					out.write(f"{record.seq[start:end]}\n")
	return
	

def main():
	args = sys.argv[1:]
	usage = 'usage: m2g_parse.py --aln relaxedphy_alignment.phy --part partitions.txt --out outdirectory/'
	#parse inputs
	if not args:
		print(usage)
		sys.exit(1)
	if args[0] == '--aln':
		aln_file = Path(args[1])
	if args[2] == '--part':
		part_file = Path(args[3])
	if args[4] == '--out':
		outdir = Path(args[5])
	else:
		print(usage)
		sys.exit(1)
	# Run 
	partitions = partitioner(part_file)
	alncutter(partitions,aln_file,outdir)


if __name__ == "__main__":
	main()
