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


def vcftoplink(pwd,basename,baseo,basep):
    os.chdir(pwd)
    vcfoutput=subprocess.check_output(shlex.split("vcftools --vcf %s.vcf --out %s --plink"%(basename,baseo)))
    ##changes the chromosome all of the SNP to chromosome 1 for forward compatibility with structure like programs.
    with open("%s_temp.map"%(basename),"w") as outfile:
        subprocess.call(shlex.split("awk '$1=1' %s.map"%(baseo)),stdout=outfile)
    print subprocess.check_output(shlex.split("mv %s_temp.map %s.map"%(basename,baseo)))
    ##Prunes the dataset
    subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --indep-pairwise 50 10 0.1"%(baseo)))
    # subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --extract plink.prune.in --recode --remove W92455.txt --out %s"%(baseo,basep)))
    subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --extract plink.prune.in --recode --geno 0.9 --mind 0.9 --out %s"%(baseo,basep)))
    return


def PCAer(pwd,basep):
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
        """%(basep,basep,basep,basep,basep))
    ### add     outliersigmathresh: 7 if you want less outliers
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


def RplotPCA(pcafile,pwd,basename,axes):
    ###Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
    ###Datafile must be tsv with headers: Samples, LAT, LON, COI[_cor]
    ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
    ro.r('dat <- read.table(file="%sdata.txt",header=TRUE)'%(pwd))
    ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
    ##Friendly random color palette!
    ro.r("""palette(c('#8dd3c7','#ffffb3','#bebada',
    '#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5',
    '#d9d9d9','#bc80bd','#ccebc5','#ffed6f','black','forestgreen','deeppink'))""")
    m2=ro.r['m2']
    # print m2
    ###Change color scheme by changing COI_col
    grdevices.png(file="%s/%s_MAP_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("LON"),m2.rx2("LAT"),col=m2.rx2("COI_5p"),main="%s_MAP_COI"%(basename),ylab="Lattitude",xlab="Longitude",pch=20,cex=2)
    grdevices.dev_off()
    grdevices.png(file="%s/%s_PCA_COI.png"%(pwd,basename), width=1000, height=1000)
    rplot(m2.rx2("V2"),m2.rx2("V3"),col=m2.rx2("COI_5p"),main="%s_PCA_COI"%(basename),ylab=axes[1],xlab=axes[0],pch=20,cex=2)
    grdevices.dev_off()
    return

def readpop(popfile):
	popdict={}
	with open(popfile, 'ru') as pop:
		poplist2 = pop.readlines()
		poplist = [p.strip('\n').strip() for p in poplist2]
		for i,pop in enumerate(poplist[0]):
			popdict[i]=pop
	return popdict
def readq(qfile):
	groupdict = {}
	with open(qfile, 'ru') as qf:
		qlines= [q.strip('\n').split() for q in qf.readlines()]
		for q,line in enumerate(qlines):
			for i,column in enumerate(line):
				if float(column)>float(0.7):
					group=i+1
					break
				else:
					group=len(column)+1
			groupdict[q]=group
	return groupdict

def admixture(pwd,basename,k):
    print "start admixture"
    os.chdir(pwd)
    cvd={}
    admixoutdir=pwd+"admixture/"
    if not os.path.exists(admixoutdir):
        os.mkdir(admixoutdir)
    shutil.copy(basename+".map",admixoutdir)
    shutil.copy(basename+".ped",admixoutdir)
    shutil.copy(basename+"_outliers.txt",admixoutdir)
    os.chdir(admixoutdir)
    ###Without removing outliers from the PCA analysis
    ###subprocess.check_output(shlex.split("plink2 --threads 7 --file %s --make-bed --out %s"%(basename,basename)))
    subprocess.check_output(shlex.split("plink2 --threads 7 --file %s --make-bed --remove %s_outliers.txt --out %s"%(basename,basename,basename)))
    for x in range(1,k):
        def admixer():
            admixcommand=shlex.split("admixture -j7 -C=0.01 --cv %s.bed %d"%(basename,x))
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
        admixer()
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
    print "Best value for K is %s with a CV of %s"%(lowestk,lowest)
    return lowestk


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
    for (dirpath, dirnames, filenames) in os.walk(pwd):
        for file in filenames:
            if file.endswith("png"):
                shutil.move(pwd+file,figpath+file)
    return

def controller(pwd,basename,k):
    baseo='%s_o'%(basename)
    basep='%s_p'%(basename)
    pcafile='%s%s.evec'%(pwd,basename)
    ##Comment out steps to skip them. PCAer() and RplotPCA() need to be run together
    vcftoplink(pwd,basename,baseo,basep)
    axes=PCAer(pwd,basep)
    RplotPCA(pcafile,pwd,basep,axes)
    # admixture(pwd,basename,k)
    # cleanup(pwd)
    pass

def main():
    # pwd=subprocess.check_output(shlex.split('pwd'))
    #basename = raw_input("Enter basename of VCF file:\n")
    pwd='/Users/josec/Desktop/Trapdoor/pyrad5/7pyrad5_outfiles/PopGenAnalysis/'
    basename="r_7pyrad5"

    installtest(pwd)
    k=10
    bs=50
    controller(pwd,basename,k+1)
    raxer(pwd,basename,bs)
    print "ALLDone"
    return

if __name__ == '__main__':
    main()
