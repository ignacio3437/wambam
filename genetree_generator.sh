#Run iqtree in parallel for each fasta alignment in a folder. Moves results to a treefile folder and concatenates them for astral


parallel --jobs 4 iqtree -s {} -nt 2 -m GTR20 ::: *.fas
mkdir treefiles
mkdir iqtree_scratch
cat *.treefile > genetrees.tre
mv *.treefile treefiles/
for f in *.fas.*; do mv $f iqtree_scratch/;done
# mv *.fas.* iqtree_scratch/
mv genetrees.tre ..




#java -jar ~/Desktop/Gitclones/ASTRAL/Astral/astral.5.5.6.jar -i /Users/josec/Desktop/exoncap/Mygal/genetrees.tre -o ~/Desktop/exoncap/Mygal/astral.tre &>astral.log
