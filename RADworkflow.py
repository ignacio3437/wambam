#!/usr/bin/python
import os,sys,re,commands,shutil
import random
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
grdevices = importr('grDevices')
rprint = ro.globalenv.get("print")
rplot=ro.r('plot')


def installtest():
    testlist=[]
    testlist.append(commands.getoutput("smartpca -v"))
    testlist.append(commands.getoutput("vcftools"))
    testlist.append(commands.getoutput("raxmlHPC-PTHREADS-SSE3 -v"))
    testlist.append(commands.getoutput("r --version"))
    testlist.append(commands.getoutput("ipyrad -v"))
    for x in testlist:
        if "command not found" in x:
            print x
            sys.exit()
    print "All programs installed"
    return


def thinner(pwd,basename):
    os.chdir(pwd)
    thin=basename+'_t'
    commands.getoutput("vcftools --vcf %s.vcf --out %s --thin 500000 --recode"%(basename,thin))
    print commands.getoutput("mv %s.recode.vcf %s.vcf"%(thin,thin))
    return thin


def vcftoevec(pwd,basename):
    os.chdir(pwd)
    baseout=basename+'_o'
    vcfoutput=commands.getoutput("vcftools --vcf %s.vcf --out %s --plink"%(basename,baseout))
    print commands.getoutput("awk '$1=1' %s.map > %s_temp.map"%(baseout,baseout))
    print commands.getoutput("mv %s_temp.map %s.map"%(baseout,baseout))

    ###Testing prune LD method
    # commands.getoutput("plink --noweb --file %s --indep-pairwise 50 10 0.1"%(baseout))
    # commands.getoutput("plink --noweb --file %s --extract plink.prune.in --recode --out %s"%(baseout,baseout))
    #

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
    ### add     outliersigmathresh: 7 if you want less outliers
    pcafile='%s%s.evec'%(pwd,baseout)
    pcaoutfile=pwd+basename+'_pcaout.txt'
    pcaout=commands.getoutput("smartpca -p pca.par")
    with open(pcaoutfile,'w') as pcaoutfile:
        pcaoutfile.write(pcaout)
    ### Pull out eigenvalues and P values. Calculate %Variation using Eigenvalues and put as tuple axes to be used for plotting PCA later.
    lines=re.findall(r"Tracy-Widom(.*)\s.*\s(.*)\s(.*)",pcaout)[0]
    eigsum=lines[0].split()[2]
    eigens=(lines[1].split()[1],lines[2].split()[1])
    pvals=(lines[1].split()[4],lines[2].split()[4])
    ipervar=[(float(eig)/float(eigsum))*100 for eig in eigens]
    pervar=["%.2f"%(i) for i in ipervar]
    axes= ("PCA1___PercentVar="+pervar[0]+"_Pval="+pvals[0],"PCA2___PercentVar="+pervar[1]+"_Pval="+pvals[1])
    ### Find outliers and print to outlierfile
    if 'REMOVED outlier' in pcaout:
        outliers=re.findall(r'REMOVED outlier (\w*)',pcaout)
    with open(pwd+baseout+'_outliers.txt', 'w') as outlierfile:
        for outlier in outliers:
            outlierfile.write(outlier+'\n')
    print "There are %d outliers in PCA analysis of %s"%(len(outliers),basename)
    ### Create keep.txt file for pyrad with outliers removed???
    return pcafile,axes


def RplotPCA(pcafile,pwd,basename,axes):
    ###Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
    ###Datafile must be tsv with headers: Samples, LAT, LON, COI[_cor]
    ro.r('evec <- read.table(file="%s%s_o.evec")'%(pwd,basename))
    ro.r('dat <- read.table(file="%sdata.txt",header=TRUE)'%(pwd))
    ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
    m2=ro.r['m2']
    # print m2
    ###Change color scheme by changing COI_col
    grdevices.png(file="%s/%s_MAP_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("LON"),m2.rx2("LAT"),col=m2.rx2("COI_cor"),main="%s_MAP_COI"%(basename),ylab="Lattitude",xlab="Longitude",pch=20,cex=1.5)
    grdevices.dev_off()
    grdevices.png(file="%s/%s_PCA_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("V2"),m2.rx2("V3"),col=m2.rx2("COI_cor"),main="%s_PCA_COI"%(basename),ylab=axes[1],xlab=axes[0],pch=20,cex=1.5)
    grdevices.dev_off()
    return


def raxer(pwd,basename,bs):
    print commands.getoutput("raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),bs,random.randint(0,999999)))
    tree="RAxML_bipartitionsBranchLabels."+basename
    return treefile


def cleanup(pwd):
    figpath=pwd+"figures/"
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    finpath=pwd+"final/"
    if not os.path.exists(finpath):
        os.mkdir(finpath)
    for file in os.listdir(pwd):
        if 'png' in file:
            shutil.move(pwd+file,figpath+file)
        else:
            pass
    return

def controller(things):
    pass

def main():
    #pwd=commands.getoutput('pwd')
    #basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis/'
    basename="7pyrad5"
    pcafile,axes=vcftoevec(pwd,basename)
    RplotPCA(pcafile,pwd,basename,axes)
    ### Comment out for Thin Options
    thin=thinner(pwd,basename)
    thinpcafile,axes2=vcftoevec(pwd,thin)
    RplotPCA(thinpcafile,pwd,thin,axes2)
    cleanup(pwd)

if __name__ == '__main__':
    main()
