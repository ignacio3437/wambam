#!/usr/bin/python


import os,sys,re,commands

# pwd=commands.getoutput('pwd')
pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis'
print "useful commands:"
basename="7pyrad5"
thin=basename+'_t'
os.chdir(pwd)
print commands.getoutput("vcftools --vcf %s.vcf --out %s_thin --thin 500000 --plink"%(basename,basename))
print commands.getoutput("awk '$1=1' %s_thin.map > %s_t.map"%(basename,basename))
print commands.getoutput("mv %s_thin.ped %s_t.ped"%(basename,basename))
with open(pwd+'/pca.par','w') as pca:
    pca.write("""
genotypename:    %s.ped
snpname:         %s.map
indivname:       %s.ped
evecoutname:     %s.evec
evaloutname:     %s.eval
altnormstyle:    NO
numoutevec:      2
familynames:     NO
grmoutname:      grmjunk
numthreads:   8
lsqproject:	YES
    """%(thin,thin,thin,thin,thin))
