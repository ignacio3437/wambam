# trimmomatic PE -phred33 -threads 8 IngoB_R1.fastq IngoB_R2.fastq CIngoB_R1_paired.fastq CIngoB_R1_unpaired.fastq CIngoB_R2_paired.fastq CIngoB_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic PE -phred33 -threads 8 Kumel_R1.fastq Kumel_R2.fastq CKumel_R1_paired.fastq CKumel_R1_unpaired.fastq CKumel_R2_paired.fastq CKumel_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic PE -phred33 -threads 8 WAMS58123_R1.fastq WAMS58123_R2.fastq CWAMS58123_R1_paired.fastq CWAMS58123_R1_unpaired.fastq CWAMS58123_R2_paired.fastq CWAMS58123_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic PE -phred33 -threads 8 WAMS68185_R1.fastq WAMS68185_R2.fastq CWAMS68185_R1_paired.fastq CWAMS68185_R1_unpaired.fastq CWAMS68185_R2_paired.fastq CWAMS68185_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic PE -phred33 -threads 8 Pteraeolidia_1_final.fastq Pteraeolidia_2_final.fastq S96364_1_paired.fastq S96364_1_unpaired.fastq S96364_2_paired.fastq CPteraeolidia_2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic SE -phred33 -threads 8 C475907.fastq CC475907.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic SE -phred33 -threads 8 C475908.fastq CC475908.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic SE -phred33 -threads 8 M15820.fastq CM15820.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic SE -phred33 -threads 8 S96364.fastq CS96364.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
# trimmomatic SE -phred33 -threads 8 WAM103087.fastq CWAM103087.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50

# #Interleave Paired Reads with \1 \2 header
# parallel sed -e "'s/ 1.*/\/1/g'" {} ">" C_{} ::: *R1*
# parallel sed -e "'s/ 2.*/\/2/g'" {} ">" C_{} ::: *R2*
# for file in IngoB Kurnell S96364 WAMS58123 WAMS68185; do seqtk mergepe ${file}_R1_paired.fastq ${file}_R2_paired.fastq > ${file}_int.fastq ;done

# # #PE
# cd /Users/josec/Desktop/Pter/MITObim/IngoB
# # perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample IngoB -ref IngoB_ref -readpool /Users/josec/Desktop/Pter/Clean/IngoB_int.fastq --quick /Users/josec/Desktop/Pter/refs/IngoB_ref.fasta -end 50 --denovo --paired --clean &> IngoB_log
# perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample IngoB -ref IngoB_ref -readpool /Users/josec/Desktop/Pter/Clean/IngoB_int.fastq --quick /Users/josec/Desktop/Pter/refs/IngoB_ref.fasta -end 50 --denovo --clean &> IngoB_log

# cd /Users/josec/Desktop/Pter/MITObim/Kurnell
# perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample Kurnell -ref Kurnell_ref -readpool /Users/josec/Desktop/Pter/Clean/Kurnell_int.fastq --quick /Users/josec/Desktop/Pter/refs/Kurnell_ref.fasta -end 50 --denovo --paired --clean &> Kurnell_log

# cd /Users/josec/Desktop/Pter/MITObim/S96364
# # perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample S96364 -ref S96364_ref -readpool /Users/josec/Desktop/Pter/Clean/S96364_int.fastq --quick /Users/josec/Desktop/Pter/refs/S96364_ref.fasta -end 50 --denovo --paired --clean &> S96364_log
# perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample S96364 -ref S96364_ref -readpool /Users/josec/Desktop/Pter/Clean/S96364_int.fastq --quick /Users/josec/Desktop/Pter/refs/S96364_ref.fasta -end 50 --denovo --clean &> S96364_log

# cd /Users/josec/Desktop/Pter/MITObim/WAMS58123
# perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample WAMS58123 -ref WAMS58123_ref -readpool /Users/josec/Desktop/Pter/Clean/WAMS58123_int.fastq --quick /Users/josec/Desktop/Pter/refs/WAMS58123_ref.fasta -end 50 --denovo --paired --clean &> WAMS58123_log

# cd /Users/josec/Desktop/Pter/MITObim/WAMS68185
# perl /Users/josec/Desktop/Gitclones/MITObim/MITObim.pl -sample WAMS68185 -ref WAMS68185_ref -readpool /Users/josec/Desktop/Pter/Clean/WAMS68185_int.fastq --quick /Users/josec/Desktop/Pter/refs/WAMS68185_ref.fasta -end 50 --denovo --paired --clean &> WAMS68185_log

# #SE
# cd /Users/josec/Desktop/Pter/MITObim/C475908
# perl ~/Desktop/Gitclones/MITObim/MITObim.pl -sample C475908 -ref C475908_ref -readpool /Users/josec/Desktop/Pter/Clean/C475908.fastq --quick /Users/josec/Desktop/Pter/refs/C475908_ref.fasta -end 100 --clean &> C475908_log

# cd /Users/josec/Desktop/Pter/MITObim/C475907
# perl ~/Desktop/Gitclones/MITObim/MITObim.pl -sample C475907 -ref C475907_ref -readpool /Users/josec/Desktop/Pter/Clean/C475907.fastq --quick /Users/josec/Desktop/Pter/refs/C475907_ref.fasta -end 100 --clean &> C475907_log

# cd /Users/josec/Desktop/Pter/MITObim/M15820
# perl ~/Desktop/Gitclones/MITObim/MITObim.pl -sample M15820 -ref M15820_ref -readpool /Users/josec/Desktop/Pter/Clean/M15820.fastq --quick /Users/josec/Desktop/Pter/refs/M15820_ref.fasta -end 100 --clean &> M15820_log

# cd /Users/josec/Desktop/Pter/MITObim/M103087
# perl ~/Desktop/Gitclones/MITObim/MITObim.pl -sample M103087 -ref M103087_ref -readpool /Users/josec/Desktop/Pter/Clean/M103087.fastq --quick /Users/josec/Desktop/Pter/refs/M103087_ref.fasta -end 100 --clean &> M103087_log




# Map the reads to the circular mitochondrial genomes and remove those reads to create a new dataset to search for other mitochondrial signal.


# mkdir /Users/josec/Desktop/Pter/bb
cd /Users/josec/Desktop/Pter/
# Map reads to all genomes to check for anything obvious.
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/C475907.fastq  interleaved=f outm=bb/C475907_allmitcir.sam covstats=C475907_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/C475908.fastq interleaved=f outm=bb/C475908_allmitcir.sam covstats=C475908_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/IngoB_int.fastq interleaved=t outm=bb/IngoB_int_allmitcir.sam covstats=IngoB_int_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/Kurnell_int.fastq interleaved=t outm=bb/Kurnell_int_allmitcir.sam covstats=Kurnell_int_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/M15820.fastq interleaved=f outm=bb/M15820_allmitcir.sam covstats=M15820_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/M103087.fastq interleaved=f outm=bb/M103087_allmitcir.sam covstats=M103087_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/S96364_int.fastq interleaved=t outm=bb/S96364_int_allmitcir.sam covstats=S96364_int_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/S96364.fastq interleaved=t outm=bb/S96364_allmitcir.sam covstats=S96364_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/WAMS58123_int.fastq interleaved=t outm=bb/WAMS58123_int_allmitcir.sam covstats=WAMS58123_int_covstat.txt
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/WAMS68185_int.fastq interleaved=t outm=bb/WAMS68185_int_allmitcir.sam covstats=WAMS68185_int_covstat.txt

#Create file of unmapped reads to check more closely
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/C475907.fastq  interleaved=f outu=bb/C475907_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/C475908.fastq interleaved=f outu=bb/C475908_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/M15820.fastq interleaved=f outu=bb/M15820_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/M103087.fastq interleaved=f outu=bb/M103087_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/S96364_int.fastq interleaved=t outu=bb/S96364_int_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/S96364.fastq interleaved=t outu=bb/S96364_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/WAMS58123_int.fastq interleaved=t outu=bb/WAMS58123_int_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/WAMS68185_int.fastq interleaved=t outu=bb/WAMS68185_int_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/IngoB_int.fastq interleaved=t outu=bb/IngoB_int_unmapped.fastq
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/Kurnell_int.fastq interleaved=t outu=bb/Kurnell_int_unmapped.fastq
#Merge Paired end data
#MAybe not worth doing....
# for file in *int.fastq;do bbmerge.sh in=$file out=Merged_${file} outu=Unmerged_${file} ihist=${file}.txt;done
#Determine best kmer size for denovo assembly using tadpole
# tadwrapper.sh in=Clean/IngoB_int.fastq out=IngoB_%_contigs.fa k=31,42,62,93 quitearly=t bisect=t

#Craps out above 5M reads for SE datasets and 10M for PE datasets. Seems strange.... Cap with reads= flag. 
#Not keeping low coverage or short contigs. k=31 seems to work well according to tadwrapper.sh
cd /Users/josec/Desktop/Pter/Clean
tadpole.sh in=S96364_int.fastq out=../tadpole/contigs_S96364_int.fastq.fa k=31 mincoverage=3 mincontig=400 reads=7000000
tadpole.sh in=IngoB_int.fastq out=../tadpole/contigs_IngoB_int.fastq.fa k=31 mincoverage=3 mincontig=400 reads=9000000
tadpole.sh in=Kurnell_int.fastq out=../tadpole/contigs_Kurnell_int.fastq.fa k=31 mincoverage=3 mincontig=400 reads=9000000
tadpole.sh in=WAMS58123_int.fastq out=../tadpole/contigs_WAMS58123_int.fastq.fa k=31 mincoverage=3 mincontig=400 reads=9000000
tadpole.sh in=WAMS68185_int.fastq out=../tadpole/contigs_WAMS68185_int.fastq.fa k=31 mincoverage=3 mincontig=400 reads=9000000
tadpole.sh in=M103087.fastq out=../tadpole/contigs_M103087.fastq.fa k=31 mincoverage=3 mincontig=400 reads=5000000
tadpole.sh in=C475907.fastq out=../tadpole/contigs_C475907.fastq.fa k=31 mincoverage=3 mincontig=400 reads=5000000
tadpole.sh in=C475908.fastq out=../tadpole/contigs_C475908.fastq.fa k=31 mincoverage=3 mincontig=400 reads=5000000
tadpole.sh in=M15820.fastq out=../tadpole/contigs_M15820.fastq.fa k=31 mincoverage=3 mincontig=400 reads=5000000



#Do close mapping to other mitType
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/C475907_unmapped.fasta  interleaved=f oumu=bb/C475907_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/C475908_unmapped.fasta interleaved=f outm=bb/C475908_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Cir_genome.fasta in=Clean/Kurnell_int_unmapped.fasta interleaved=t outm=bb/Kurnell_int_vslow_othertypeMitMap.sam


# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type1_genomes.fasta in=Clean/IngoB_int_unmapped.fasta interleaved=t outm=bb/IngoB_int_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type1_genomes.fasta in=Clean/M15820_unmapped.fasta interleaved=f outm=bb/M15820_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type1_genomes.fasta in=Clean/M103087_unmapped.fasta interleaved=f outm=bb/M103087_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type1_genomes.fasta in=Clean/S96364_int_unmapped.fasta interleaved=t outm=bb/S96364_int_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type1_genomes.fasta in=Clean/S96364_unmapped.fasta interleaved=t outm=bb/S96364_vslow_othertypeMitMap.sam


# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type2_genomes.fasta in=Clean/WAMS58123_int_unmapped.fasta interleaved=t outm=bb/WAMS58123_int_vslow_othertypeMitMap.sam
# bbmap.sh vslow=t k=8 maxindel=200 minratio=0.1 local=t nodisk=t ref=CircularGenomes/Type2_genomes.fasta in=Clean/WAMS68185_int_unmapped.fasta interleaved=t outm=bb/WAMS68185_int_vslow_othertypeMitMap.sam





# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/SulawesiA.fasta in=Clean/WAMS58123_int_unmapped.fasta interleaved=t outm=bb/WAMS58123_int_fast_sualwesiA.sam
# bbmap.sh maxindel=100 semiperfectmode=f local=t nodisk=t ref=CircularGenomes/SulawesiA.fasta in=Clean/WAMS68185_int_unmapped.fasta interleaved=t outm=bb/WAMS68185_int_fast_sualwesiA.sam



# # Merge overlapping paired end data



