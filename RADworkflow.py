#!/usr/bin/python
import os,sys,re,subprocess,shutil,shlex
import random
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
grdevices = importr('grDevices')
rprint = ro.globalenv.get("print")
rplot=ro.r('plot')


def installtest(pwd):
    os.chdir(pwd)
    subprocess.check_output(shlex.split("smartpca -v"))
    subprocess.check_output(shlex.split("vcftools"))
    subprocess.check_output(shlex.split("raxmlHPC-PTHREADS-SSE3 -v"))
    subprocess.check_output(shlex.split("r --version"))
    subprocess.check_output(shlex.split("ipyrad --version"))
    subprocess.check_output(shlex.split("admixture --help"))
    subprocess.check_output(shlex.split("plink2 --version"))
    print "All programs installed"
    return


def loaddata(pwd,datafile):
    ###Datafile must be tsv with headers: Samples, LAT, LON, COI[_cor]
    ro.r('dat <- read.table(file="%s%s",header=TRUE)'%(pwd,datafile))
    ### Friendly random color palette!
    ro.r("""palette(c('#8dd3c7','#ffffb3','#bebada',
    '#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5',
    '#d9d9d9','#bc80bd','#ccebc5','#ffed6f','black',
    'forestgreen','brown','deeppink'))""")
    return


def vcftoplink(pwd,basename,baseo,basep):
    os.chdir(pwd)
    vcfoutput=subprocess.check_output(shlex.split("vcftools --vcf %s.vcf --out %s --plink"%(basename,baseo)))
    ##changes the chromosome all of the SNP to chromosome 1 for forward compatibility with structure like programs.
    with open("%s_temp.map"%(basename),"w") as outfile:
        subprocess.call(shlex.split("awk '$1=1' %s.map"%(baseo)),stdout=outfile)
    print subprocess.check_output(shlex.split("mv %s_temp.map %s.map"%(basename,baseo)))
    ##Prunes the dataset
    subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --indep-pairwise 50 10 0.1"%(baseo)))
    subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --extract plink.prune.in --recode --geno 0.9 --mind 0.9 --maf 0.05 --out %s"%(baseo,basep)))
    return


def PCAer(pwd,basep):
    with open(pwd+'pca.par','w') as pca:
        pca.write("""
    genotypename:    %s.ped
    snpname:         %s.map
    indivname:       %s.ped
    evecoutname:     %s.evec
    evaloutname:     %s.eval
    altnormstyle:    YES
    numoutevec:      2
    familynames:     NO
    grmoutname:      grmjunk
    numthreads:   8
    lsqproject:    YES
        """%(basep,basep,basep,basep,basep))
    ### add     outliersigmathresh: 7 if you want more outliers
    pcaoutfile=pwd+basep+'_pcaout.txt'
    pcaout=subprocess.check_output(shlex.split("smartpca -p pca.par"))
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
    with open(pwd+basep+'_outliers.txt', 'w') as outlierfile:
        for outlier in outliers:
            outlierfile.write(outlier+'\t'+outlier+'\n')
    print "There are %d outliers in PCA analysis of %s"%(len(outliers),basep)
    return axes


def Rplot(pwd,basename,type):
    admixoutdir=pwd+"admixture/"
    if "PCA" in type[1]:
        ##Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
        ##Change coloring by searching from col=m2.rx2("XXXX")
        ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
        ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
        dat=ro.r['dat']
        m2=ro.r['m2']
        grdevices.png(file="%s%s_MAP_COI.png"%(pwd,basename), width=1000, height=1000)
        rplot(m2.rx2("LON"),m2.rx2("LAT"),col=m2.rx2("COI_13"),main="%s_MAP_COI"%(basename),ylab="Lattitude",xlab="Longitude",pch=20,cex=2)
        # ro.r('text(m2$LON,m2$LAT,m2$Sample,cex=0.8)')
        grdevices.dev_off()
        grdevices.png(file="%s%s_PCA_COI.png"%(pwd,basename), width=1000, height=1000)
        rplot(m2.rx2("V2"),m2.rx2("V3"),col=m2.rx2("COI_13"),main="%s_PCA_COI"%(basename),ylab=type[1],xlab=type[0],pch=20,cex=2,)
        grdevices.dev_off()
    elif "CV" in type:
        ##Print CV error plot to determine best K from admixture.
        ro.r('cvs <- read.table(file="cv.txt",header=FALSE)')
        cvs=ro.r['cvs']
        grdevices.png(file="%s/CVplot.png"%(pwd), width=1000, height=1000)
        rplot(cvs.rx2("V1"),cvs.rx2("V2"),main="CVplot",ylab="CVeror",xlab="K",pch=20,cex=1.5)
        grdevices.dev_off()
        os.chdir(pwd)
    elif "admix" in type:
        ro.r('palette(sample(rainbow(10)))')
        ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
        ro.r('am2 <- merge(am2,evec,all.x=TRUE,by.x="Sample",by.y="V1")')
        am2=ro.r['am2']
        grdevices.png(file="%s%s_MAP_AdmixGroup.png"%(admixoutdir,basename), width=1000, height=1000)
        rplot(am2.rx2("LON"),am2.rx2("LAT"),col=am2.rx2("Sgroup"),main="%s_MAP_AdmixGroup"%(basename),ylab="Lattitude",xlab="Longitude",pch=20,cex=2)
        grdevices.dev_off()
        grdevices.png(file="%s%s_PCA_AdmixGroup.png"%(admixoutdir,basename), width=1000, height=1000)
        rplot(am2.rx2("V2"),am2.rx2("V3"),col=am2.rx2("Sgroup"),main="%s_PCA_COI"%(basename),ylab=type[1],xlab=type[0],pch=20,cex=2,)
        grdevices.dev_off()
    return


def admixture(pwd,base,k):
    print "start admixture"
    os.chdir(pwd)
    cvd={}
    admixoutdir=pwd+"admixture/"
    if not os.path.exists(admixoutdir):
        os.mkdir(admixoutdir)
    shutil.copy(base+".map",admixoutdir)
    shutil.copy(base+".ped",admixoutdir)
    shutil.copy(base+"_outliers.txt",admixoutdir)
    shutil.copy("data.txt",admixoutdir)
    os.chdir(admixoutdir)
    subprocess.check_output(shlex.split("plink2 --threads 7 --file %s --make-bed --out %s"%(base,base)))
    ## make this a new function...?
    ## subprocess.check_output(shlex.split("plink2 --threads 7 --file %s --make-bed --remove %s_outliers.txt --out %s"%(base,base,base)))
    for x in range(1,k):
        def admixer():
            admixcommand=shlex.split("admixture -j7 -C=0.01 --cv %s.bed %d"%(base,x))
            popen = subprocess.Popen(admixcommand, stdout=subprocess.PIPE, universal_newlines=True)
            with open(admixoutdir+"admixout.txt",'a') as admixoutfile:
                admixoutfile.write("####################K=%s\n####################\n"%(x))
                for stdout_line in iter(popen.stdout.readline, ""):
                    # print stdout_line
                    admixoutfile.write(stdout_line)
                    if "nan" in stdout_line:
                        return
                    elif "CV" in stdout_line:
                        cv=re.findall(r'CV .*K=(\d*).*: (\d*\.\d*)',stdout_line)
                        cvd[x]=cv[0][1]
                popen.stdout.close()
                return_code = popen.wait()
                if return_code:
                    raise subprocess.CalledProcessError(return_code, admixcommand)
            print "Admixture: Finished K=%d"%(x)
            return
        admixer()
    lowest=100
    lowestk="k"
    with open(admixoutdir+"cv.txt","w") as cvout:
        for key in cvd:
            cvout.write("%s\t%s\n"%(key,cvd[key]))
            if float(cvd[key])<lowest:
                lowest=float(cvd[key])
                lowestk=key
    print "Best value for K is %s with a CV of %s"%(lowestk,lowest)
    return lowestk


def plotadmix(pwd,basename,lowestk):
    admixoutdir=pwd+"admixture/"
    os.chdir(admixoutdir)
    taxdict={}
    qdict = {}
    qfile="%s%s.%d.Q"%(admixoutdir,basename,lowestk)
    Rplot(admixoutdir,basename,"CV")
    with open("%s%s.fam"%(admixoutdir,basename),'ru') as ns:
        nslines=ns.readlines()
        tlist=[t.strip().split(" ")[0] for t in nslines]
        for i,taxa in enumerate(tlist):
            taxdict[i]=taxa
    with open(qfile, 'ru') as qf:
        qlines= [q.strip('\n').split() for q in qf.readlines()]
        for q,line in enumerate(qlines):
            for i,column in enumerate(line):
                if float(column)>float(0.5):
                    group=i+1
                    break
                else:
                    group=len(column)+1
            qdict[taxdict[q]]=group
    aoutname="%sassigned.%d.txt"%(admixoutdir,lowestk)
    with open(aoutname,'w') as aout:
        aout.write("Sample\tSgroup\n")
        for key in qdict.keys():
            aout.write("%s\t%d\n"%(key,qdict[key]))
    ro.r('aout <- read.table(file="%s",header=TRUE)'%(aoutname))
    ro.r('am2 <- merge(dat,aout,by.x="Sample",by.y="Sample")')
    Rplot(pwd,basename,"admix")
    return


def raxer(pwd,basename,bs):
    # print subprocess.check_output(shlex.split("raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),bs,random.randint(0,999999))))
    print "raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),bs,random.randint(0,999999))
    # tree="RAxML_bipartitionsBranchLabels."+basename
    return

def cleanup(pwd):
    figpath=pwd+"figures/"
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    finpath=pwd+"final/"
    if not os.path.exists(finpath):
        os.mkdir(finpath)
    for root, dirs, files in os.walk(pwd):
        for file in files:
            path=os.path.join(root,file)
            if file.endswith("png") and "figures" not in path:
                shutil.copy(path,figpath)
            elif "flag" in path and "final" not in paht:
                shutil.copy(path,finpath)
    return

def controller(pwd,basename,k):
    baseo='%s_o'%(basename)
    basep='%s_p'%(basename)
    pcafile='%s%s.evec'%(pwd,basename)
    ##Comment out steps to skip them. PCAer() and Rplot() need to be run together
    vcftoplink(pwd,basename,baseo,basep)
    axes=PCAer(pwd,basep)
    Rplot(pwd,basep,axes)
    lowestk=admixture(pwd,basep,k)
    # lowestk=2
    plotadmix(pwd,basep,lowestk)
    ro.r("write.table(am2, file='%smergeddata.txt', quote=FALSE, row.names=FALSE, sep='\t')"%(pwd))
    cleanup(pwd)
    return

def main():
    # pwd=subprocess.check_output(shlex.split('pwd'))
    #basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/freshstart/red/red_outfiles/'
    basename="red"
    datafile="data.txt"
    # installtest(pwd)
    loaddata(pwd,datafile)
    k=5
    bs=10
    controller(pwd,basename,k+1)
    # raxer(pwd,basename,bs)
    print "ALLDone"
    return

if __name__ == '__main__':
    main()
