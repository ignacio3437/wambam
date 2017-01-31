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


def thinner(pwd,basename):
    os.chdir(pwd)
    thin=basename+'_t'
    print commands.getoutput("vcftools --vcf %s.vcf --out %s --thin 500000 --recode"%(basename,thin))
    print commands.getoutput("mv %s.recode.vcf %s.vcf"%(thin,thin))
    return thin


def vcftoPCA(pwd,basename):
    os.chdir(pwd)
    baseout=basename+'_o'
    print commands.getoutput("vcftools --vcf %s.vcf --out %s --plink"%(basename,baseout))
    print commands.getoutput("awk '$1=1' %s.map > %s_temp.map"%(baseout,baseout))
    print commands.getoutput("mv %s_temp.map %s.map"%(baseout,baseout))
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
    outliersigmathresh: 7
        """%(baseout,baseout,baseout,baseout,baseout))
    pcafile='%s%s.evec'%(pwd,baseout)
    pcaout=commands.getoutput("smartpca -p pca.par")
    with open(pwd+'pcaout.txt','a') as pcaoutfile:
        pcaoutfile.write(pcaout)
    return pcafile


def RplotPCA(pcafile,pwd,basename):
    #Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
    #Datafile must be tsv with headers: Samples, LAT, LON, COI[_cor]
    ro.r('evec <- read.table(file="%s%s_o.evec")'%(pwd,basename))
    ro.r('dat <- read.table(file="%sdata.txt",header=TRUE)'%(pwd))
    ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
    m2=ro.r['m2']
    ####
    #Change color scheme by changing COI_cor
    ####
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
    pcafile=vcftoPCA(pwd,basename)
    RplotPCA(pcafile,pwd,basename)
    ####
    #Comment out for Thin Options
    ####
    thin=thinner(pwd,basename)
    thinpcafile=vcftoPCA(pwd,thin)
    RplotPCA(thinpcafile,pwd,thin)

if __name__ == '__main__':
    main()
