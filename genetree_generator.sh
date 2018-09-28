#remove alignments with less than x sequences
for f in *.fas;do a=($(wc $f));lines=${a[0]};if [ $lines -gt 30 ]; then mv $f test2; fi;done


#Run iqtree in parallel for each fasta alignment in a folder. Moves results to a treefile folder and concatenates them for astral


parallel --jobs 4 iqtree -s {} -nt 2 -m TEST -bb 1000 -wbt ::: *.fasta
mkdir treefiles
mkdir ufboot
mkdir iqtree_scratch
cat *.treefile > genetrees.tre
ls $PWD/*.ufboot > bootfiles.txt
mv *.ufboot > ufboot
ls $PWD/* > bootfiles.txt
mv *.treefile treefiles/
for f in *.fas.*; do mv $f iqtree_scratch/;done
# mv *.fas.* iqtree_scratch/
mv genetrees.tre ..




#java -jar ~/Desktop/Gitclones/ASTRAL/Astral/astral.5.5.6.jar -i /Users/josec/Desktop/exoncap/Mygal/genetrees.tre -o ~/Desktop/exoncap/Mygal/astral.tre &>astral.log
java -jar /home/gs_junior/ASTRAL-master/astral.5.5.9.jar -i genetrees.tre -o astral.tre &>astral.log
java -jar /Users/josec/Desktop/Gitclones/ASTRAL/astral.5.6.2.jar -i genetrees.tre -o astral.tre -b bs-files -r 1000 &>astral.log