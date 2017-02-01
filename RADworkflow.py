#!/usr/bin/python
import os,sys,re,commands
import random
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
grdevices = importr('grDevices')
rprint = ro.globalenv.get("print")
rplot=ro.r('plot')

def controller(things):
    pass


def thinner(pwd,basename):
    os.chdir(pwd)
    thin=basename+'_t'
    print commands.getoutput("vcftools --vcf %s.vcf --out %s --thin 500000 --recode"%(basename,thin))
    print commands.getoutput("mv %s.recode.vcf %s.vcf"%(thin,thin))
    return thin


def vcftoevec(pwd,basename):
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
        """%(baseout,baseout,baseout,baseout,baseout))
    # add     outliersigmathresh: 7 if you want less outliers
    pcafile='%s%s.evec'%(pwd,baseout)
    pcaoutfile=pwd+basename+'_pcaout.txt'
    pcaout=commands.getoutput("smartpca -p pca.par")
    with open(pcaoutfile,'w') as pcaoutfile:
        pcaoutfile.write(pcaout)
    if 'REMOVED outlier' in pcaout:
        outlierlist=re.findall(r'REMOVED outlier (\w*)',pcaout)
    with open(pwd+baseout+'_outliers.txt', 'w') as outlierfile:
        for outlier in outlierlist:
            outlierfile.write(outlier+'\n')
    return pcafile


def RplotPCA(pcafile,pwd,basename):
    #Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
    #Datafile must be tsv with headers: Samples, LAT, LON, COI[_cor]
    ro.r('evec <- read.table(file="%s%s_o.evec")'%(pwd,basename))
    ro.r('dat <- read.table(file="%sdata.txt",header=TRUE)'%(pwd))
    ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
    m2=ro.r['m2']
    # print m2
    ####
    #Change color scheme by changing COI_cor
    ####
    grdevices.png(file="%s/%s_MAP_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("LON"),m2.rx2("LAT"),col=m2.rx2("Outlier"),main="%s_MAP_COI"%(basename),ylab="Lattitude",xlab="Longitude",pch=20,cex=1.5)
    grdevices.png(file="%s/%s_PCA_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("V2"),m2.rx2("V3"),col=m2.rx2("Outlier"),main="%s_PCA_COI"%(basename),ylab="V3",xlab="V2",pch=20,cex=1.5)
    grdevices.dev_off()
    return


def raxer(pwd,basename,bs):
    print commands.getoutput("raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),random.randint(0,999999),bs))
    tree="RAxML_bipartitionsBranchLabels."+basename
    return treefile

def main():
    #pwd=commands.getoutput('pwd')
    #basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis/'
    basename="7pyrad5"
    pcafile=vcftoevec(pwd,basename)
    RplotPCA(pcafile,pwd,basename)
    # Comment out for Thin Options
    thin=thinner(pwd,basename)
    thinpcafile=vcftoevec(pwd,thin)
    RplotPCA(thinpcafile,pwd,thin)

if __name__ == '__main__':
    main()
