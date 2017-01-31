#!/usr/bin/python
import os,sys,re,commands
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
grdevices = importr('grDevices')
rprint = ro.globalenv.get("print")
rplot=ro.r('plot')


def vcftoPCA(pwd,basename,thin):
    os.chdir(pwd)
    # print commands.getoutput("vcftools --vcf %s.vcf --out %s_thin --thin 500000 --plink"%(basename,basename))
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
    # pcaout=commands.getoutput("smartpca -p pca.par")
    # with open(pwd+'pcaout.txt','w') as pcaoutfile:
    #     pcaoutfile.write(pcaout)
    return pcafile


def RplotPCA(pcafile,pwd,basename):
    #Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
    ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
    ro.r('dat <- read.table(file="%sdata.txt",header=TRUE)'%(pwd))
    ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
    m2=ro.r['m2']
    print m2.colnames
    ##
    #Change color scheme by changing COI_cor
    ##
    grdevices.png(file="%s/%s_MAP_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("LON"),m2.rx2("LAT"),col=m2.rx2("COI_cor"),main="%s_MAP_COI"%(basename),ylab="Lattitude",xlab="Longitude",pch=20,cex=1.5)
    grdevices.png(file="%s/%s_PCA_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("V2"),m2.rx2("V3"),col=m2.rx2("COI_cor"),main="%s_PCA_COI"%(basename),ylab="V3",xlab="V2",pch=20,cex=1.5)
    grdevices.dev_off()
    return



def main():
    # pwd=commands.getoutput('pwd')
    # basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis/'
    basename="7pyrad5"
    thin=basename+'_t'
    pcafile=vcftoPCA(pwd,basename,thin)
    RplotPCA(pcafile,pwd,thin)

if __name__ == '__main__':
    main()
