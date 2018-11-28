


wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR1598945/SRR1598945.sra
#for Trintiy assemblies
SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1598945.sra
#as fasta not fastq
parallel -q fastq-dump.2 --defline-seq '@$sn[_$rn]/$ri' --fasta 0 --split-files {} ::: *.sra

# paired reads interleaved with /1/2 ending in name
fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' -Z SRR1598945.sra

