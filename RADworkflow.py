#!/usr/bin/python
import os,sys,re,commands,shutil
import random
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
# from multiprocessing import Pool, TimeoutError
grdevices = importr('grDevices')
rprint = ro.globalenv.get("print")
rplot=ro.r('plot')


def installtest(pwd):
    testlist=[]
    os.chdir(pwd)
    testlist.append(commands.getoutput("smartpca -v"))
    testlist.append(commands.getoutput("vcftools"))
    testlist.append(commands.getoutput("raxmlHPC-PTHREADS-SSE3 -v"))
    testlist.append(commands.getoutput("r --version"))
    testlist.append(commands.getoutput("ipyrad -v"))
    testlist.append(commands.getoutput("admixture -v"))
    testlist.append(commands.getoutput("plink2 -v"))

    for x in testlist:
        if "command not found" in x:
            print x
            sys.exit()
    print "All programs installed"
    return


def vcftoevec(pwd,basename,prune):
    os.chdir(pwd)
    original=basename
    vcfoutput=commands.getoutput("vcftools --vcf %s.vcf --out %s_o --plink"%(basename,basename))
    basename=original+'_o'
    print commands.getoutput("awk '$1=1' %s.map > %s_temp.map"%(basename,basename))
    print commands.getoutput("mv %s_temp.map %s.map"%(basename,basename))
    if prune=="ON":
        baseprune=original+'_p'
        commands.getoutput("plink2 --file %s --threads 7 --indep-pairwise 50 10 0.1"%(basename))
        commands.getoutput("plink2 --file %s --threads 7 --extract plink.prune.in --recode --remove W92455.txt --out %s"%(basename,baseprune))
        basename=baseprune
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
        """%(basename,basename,basename,basename,basename))
    ### add     outliersigmathresh: 7 if you want less outliers
    pcafile='%s%s.evec'%(pwd,basename)
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
    else:
        outliers=[]
    with open(pwd+basename+'_outliers.txt', 'w') as outlierfile:
        for outlier in outliers:
            outlierfile.write(outlier+'\t'+outlier+'\n')
    print "There are %d outliers in PCA analysis of %s"%(len(outliers),basename)
    ### Create keep.txt file for pyrad with outliers removed???
    return pcafile,axes,basename


def RplotPCA(pcafile,pwd,basename,axes):
    ###Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
    ###Datafile must be tsv with headers: Samples, LAT, LON, COI[_cor]
    ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
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


def admixture(pwd,basename,k):
    os.chdir(pwd)
    admixout=[]
    cvd={}
    admixoutdir=pwd+"admixture/"
    if not os.path.exists(admixoutdir):
        os.mkdir(admixoutdir)
    shutil.copy(basename+".map",admixoutdir)
    shutil.copy(basename+".ped",admixoutdir)
    shutil.copy(basename+"_outliers.txt",admixoutdir)
    os.chdir(admixoutdir)
    ###Without removing outliers from the PCA analysis
    ###commands.getoutput("plink2 --threads 7 --file %s --make-bed --out %s"%(basename,basename))
    commands.getoutput("plink2 --threads 7 --file %s --make-bed --geno 0.99 --remove %s_outliers.txt --out %s"%(basename,basename,basename))
    for x in range(1,k):
        output = commands.getoutput("admixture -j7 -C=0.01 --cv %s.bed %d"%(basename,x))
        print "Admixture: Finished K=%d"%(x)
        admixout.append(output)
        cv=re.findall(r'CV .*K=(\d*).*: (\d*\.\d*)',output)
        cvd[x]=cv[0][1]
    with open(admixoutdir+"admixout.txt",'w') as admixoutfile:
        admixoutfile.write("\n\n".join(admixout))
    lowest=100
    lowestk="k"
    with open(admixoutdir+"cv.txt","w") as cvout:
        for key in cvd:
            cvout.write("%s\t%s\n"%(key,cvd[key]))
            if float(cvd[key])<lowest:
                lowest=float(cvd[key])
                lowestk=key
    ro.r('cvs <- read.table(file="%s",header=FALSE)'%(admixoutdir+"cv.txt"))
    cvs=ro.r['cvs']
    grdevices.png(file="%s/CVplot.png"%(admixoutdir), width=1000, height=1000)
    rplot(cvs.rx2("V1"),cvs.rx2("V2"),main="CVplot",ylab="CVeror",xlab="K",pch=20,cex=1.5)
    grdevices.dev_off()
    os.chdir(pwd)
    return lowestk,admixoutdir


def raxer(pwd,basename,bs):
    print commands.getoutput("raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),bs,random.randint(0,999999)))
    # print "raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),bs,random.randint(0,999999))
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
    for file in os.listdir(pwd+"/admixture/"):
        if 'png' in file:
            shutil.move(pwd+"/admixture/"+file,figpath+file)

    return

def controller(pwd,basename,prune,k):
    pcafile,axes,basename=vcftoevec(pwd,basename,prune)
    RplotPCA(pcafile,pwd,basename,axes)
    lowestk,admixoutdir=admixture(pwd,basename,k)
    cleanup(pwd)
    pass

def main():
    # pwd=commands.getoutput('pwd')
    #basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis/'
    installtest(pwd)
    basename="7pyrad5"
    prune="ON"
    k=30

    controller(pwd,basename,prune,k+1)
    # prune="OFF"
    # controller(pwd,basename,prune)
    return

if __name__ == '__main__':
    main()


    """
    # pool = Pool(processes=4)
    # pool.apply_async(
    # pool.close()
    # pool.join()


    # Pool(processes=4).map(admixture,range(10))"""
