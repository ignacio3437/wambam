
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR838/SRR8386736/SRR8386736.sra
#for Trintiy assemblies
parallel -q fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files {} ::: *.sra
#as fasta not fastq
parallel -q fastq-dump.2 --defline-seq '@$sn[_$rn]/$ri' --fasta 0 --split-files {} ::: *.sra

# paired reads interleaved with /1/2 ending in name
parallel -q fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' -Z {} ::: *.sra
fastq-dump.2 --defline-seq '@$sn[_$rn]/$ri' --split-files -A SRR8386726 