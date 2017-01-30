#!/usr/bin/python
import os,sys,re,commands
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com




def vcftoPCA(arg):
    os.chdir(pwd)
    thin=basename+'_t'
    print commands.getoutput("vcftools --vcf %s.vcf --out %s_thin --thin 500000 --plink"%(basename,basename))
    print commands.getoutput("awk '$1=1' %s_thin.map > %s_t.map"%(basename,basename))
    print commands.getoutput("mv %s_thin.ped %s_t.ped"%(basename,basename))
    with open(pwd+'pca.par','w') as pca:
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
    pcafile='%s%s.evec'%(pwd,thin)
    pcaout=commands.getoutput("smartpca -p pca.par")
    with open(pwd+'pcaout.txt','w') as pcaoutfile:
        pcaoutfile.write(pcaout)
    return pcaoutfile


def RplotPCA(arg):
    pass


def main():
    # pwd=commands.getoutput('pwd')
    # basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis/'
    basename="7pyrad5"
    pcaoutfile=vcftoPCA(pwd,basename)

if __name__ == '__main__':
    main()
