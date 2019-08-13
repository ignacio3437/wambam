trimmomatic PE -phred33 -threads 8 IngoB_R1.fastq IngoB_R2.fastq CIngoB_R1_paired.fastq CIngoB_R1_unpaired.fastq CIngoB_R2_paired.fastq CIngoB_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic PE -phred33 -threads 8 Kumel_R1.fastq Kumel_R2.fastq CKumel_R1_paired.fastq CKumel_R1_unpaired.fastq CKumel_R2_paired.fastq CKumel_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic PE -phred33 -threads 8 WAMS58123_R1.fastq WAMS58123_R2.fastq CWAMS58123_R1_paired.fastq CWAMS58123_R1_unpaired.fastq CWAMS58123_R2_paired.fastq CWAMS58123_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic PE -phred33 -threads 8 WAMS68185_R1.fastq WAMS68185_R2.fastq CWAMS68185_R1_paired.fastq CWAMS68185_R1_unpaired.fastq CWAMS68185_R2_paired.fastq CWAMS68185_R2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50

cd /Users/josec/Desktop/Pter/Raw
trimmomatic PE -phred33 -threads 8 Pteraeolidia_1_final.fastq Pteraeolidia_2_final.fastq S96364_1_paired.fastq S96364_1_unpaired.fastq S96364_2_paired.fastq CPteraeolidia_2_unpaired.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50


trimmomatic SE -phred33 -threads 8 C475907.fastq CC475907.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic SE -phred33 -threads 8 C475908.fastq CC475908.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic SE -phred33 -threads 8 M15820.fastq CM15820.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic SE -phred33 -threads 8 S96364.fastq CS96364.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
trimmomatic SE -phred33 -threads 8 WAM103087.fastq CWAM103087.fastq ILLUMINACLIP:/usr/local/opt/trimmomatic/share/trimmomatic/adapters/NexteraSE.fa:2:30:10 LEADING:7 TRAILING:7 SLIDINGWINDOW:4:15 MINLEN:50
echo '\b'