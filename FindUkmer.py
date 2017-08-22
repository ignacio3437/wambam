#!/usr/bin/python

"""
Returns a list of kmers the length of sliding_window that are unique in
your fasta file. Useful for finding genetic diagnostic markers.
"""


input='/Users/josec/Desktop/tester.fa'
sliding_window=5
# start_position=0
uniqe_Kmer_list=[]
def fastadicter(input):
	# Make a python dictionary out of fasta file
	# Key==Sequence_name 	Value==Sequence
	inputf= open(input, 'ru')
	fastdict={}
	faslines = inputf.readlines()
	for i,line in enumerate(faslines):
		if i ==0:
			nel=0
		if '>' in line:
			if nel:
				fastdict[fasname]=nel
			fasname = line.strip('>').rstrip('\n')
			nel=''
		else:
			nel+=line.rstrip('\n')
		if i ==len(faslines)-1:
			fastdict[fasname]=nel
		fastdict[fasname]=nel
	inputf.close()
	return fastdict

def generateKmerdict(fastdict,start_position):
	# Cuts up the sequences into sliding_window sized slices
	# kmerdict key ==  start index of slice 	value == kmer
	# Master dictionary is kmerdict for each sequence in fasta
	master_dict={}
	for key in fastdict.keys():
		seq=fastdict[key]
		kmer_list=[]
		start_index_list = range(len(seq))[start_position::sliding_window]
		for start in start_index_list:
			kmer_list.append(seq[start:start+sliding_window])
		kmerdict = dict(zip(range(len(kmer_list)),kmer_list))
		master_dict[key] = kmerdict
	return master_dict

def findUniqeKmers(master_dict):
	all_kmers=[]
	for kmer_dict in master_dict.values():
		for kmer in kmer_dict.values():
			all_kmers.append(kmer)
	for ukmer in list(set(all_kmers)):
		if all_kmers.count(ukmer)==1:
			uniqe_Kmer_list.append(ukmer)
	return

def main():
	for start_position in range(sliding_window):
		fastdict=fastadicter(input)
		master_dict=generateKmerdict(fastdict,start_position)
		findUniqeKmers(master_dict)
	print uniqe_Kmer_list
	for start_position in range(sliding_window):
		fastdict=fastadicter(input)
		master_dict=generateKmerdict(fastdict,start_position)
		for kmer in uniqe_Kmer_list:
			for name in master_dict.keys():
				if kmer in master_dict[name].values():
					print name,kmer
					pass
	return

if __name__ == '__main__':
	main()
