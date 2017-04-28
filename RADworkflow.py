#!/usr/bin/python
import os,sys,re,subprocess,shutil,shlex
import random
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
ggmap = importr('ggmap')
ggplot = importr('ggplot2')
scatter3d= importr('scatterplot3d')
RColor=importr('RColorBrewer')
grdevices = importr('grDevices')
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
    ###Datafile must be tsv with headers: Samples, LAT, LON, COI_both
    ro.r('dat <- read.table(file="%s%s",header=TRUE,sep="\t")'%(pwd,datafile))
    return


def vcftoplink(pwd,basename,baseo,basep,taxon_cov):
    os.chdir(pwd)
    subprocess.check_output(shlex.split("vcftools --vcf %s.vcf --max-missing %.2f --out %s --plink"%(basename,taxon_cov,baseo)))
    ##change Remove FUNCTION --remove-indv WAMM57947
    with open("%s_temp.map"%(basename),"w") as outfile:
        subprocess.call(shlex.split("awk '$1=1' %s.map"%(baseo)),stdout=outfile)
    print subprocess.check_output(shlex.split("mv %s_temp.map %s.map"%(basename,baseo)))
    ##Prunes the dataset
    subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --indep-pairwise 50 10 0.1"%(baseo)))
    subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --extract plink.prune.in --recode --freq --maf 0.001 --out %s"%(baseo,basep)))
    # subprocess.check_output(shlex.split("plink2 --file %s --threads 7 --extract plink.prune.in --recode --out %s"%(baseo,basep)))
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
    numoutevec:      3
    familynames:     NO
    grmoutname:      grmjunk
    numthreads:   8
    lsqproject:    YES
        """%(basep,basep,basep,basep,basep))
    ### add     outliersigmathresh: 7 if you want more outliers
    pcaoutfile=pwd+basep+'_pcaoutI.txt'
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
    with open(pwd+basep+'_outliersI.txt', 'w') as outlierfile:
        for outlier in outliers:
            outlierfile.write(outlier+'\t'+outlier+'\n')
    print "There are %d outliers in PCA analysis of %s\n\n\n\n"%(len(outliers),basep)
    return axes


def Rplot(pwd,basename,type,lowestk,taxon_cov):
    admixoutdir=pwd+"admixture/"
    taxon_cov=str(taxon_cov).replace('0.','')
    if "PCA" in type[1]:
        ##Merges PCA output and Datafile in R then prints XY plots as PNG in pwd.
        ##Change coloring by searching from col=m2.rx2("XXXX")
        ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
        ro.r('m2 <- merge(dat,evec,by.x="Sample",by.y="V1")')
        dat=ro.r['dat']
        m2=ro.r['m2']
        grdevices.png(file="%s%s_MAP_Pop_%sTC.png"%(pwd,basename,taxon_cov), width=1000, height=1000)
        ro.r('lbound<-c(min(dat$LON),min(dat$LAT));ubound<-c(max(dat$LON),max(dat$LAT));bounds<-c(lbound,ubound);Category<-factor(dat$COI_both);options(warn=-1)')
        ro.r('map<-get_map(location=bounds,zoom=7,maptype="satellite");options(warn=0)')
        ro.r('plot(ggmap(map)+ggtitle("Map by Pop")+geom_point(data = dat,aes(LON,LAT,colour=Category,size=1),show.legend=T)+scale_color_brewer(palette="Set3"))')
        grdevices.dev_off()
        ro.r('png("%s%s_PCA_Pop_%sTC.png", width=1000, height=1000)'%(pwd,basename,taxon_cov))
        ro.r('myplot<-ggplot(m2,aes(x=V2,y=V3,color=factor(COI_both)))+ggtitle("PCA by Pop")+geom_point(size=4,show.legend=F)+scale_color_brewer(palette="Set3")+theme_dark(base_size=16)+labs(x="%s",y="%s")'%(type[0],type[1]))
        ro.r('print(myplot)')
        grdevices.dev_off()
        ro.r('palette(brewer.pal(12,"Set3"))')
        ro.r('png("%s%s_PCA3D_Pop_%sTC.png", width=1000, height=1000)'%(pwd,basename,taxon_cov))
        ro.r('myplot3D<-scatterplot3d(m2$V2,m2$V3,m2$V4,pch=16,color=as.numeric(m2$COI_both))')
        ro.r('print(myplot3D)')
        grdevices.dev_off()
    elif "CV" in type:
        ##Print CV error plot to determine best K from admixture.
        ro.r('cvs <- read.table(file="cv.txt",header=FALSE)')
        cvs=ro.r['cvs']
        ro.r('png("CV_plot_%sTC.png", width=1000, height=1000)'%(taxon_cov))
        ro.r('cvplot<-ggplot(cvs,aes(x=V1,y=V2))+geom_point(size=4)+geom_line(col="red")+theme_dark()+theme(text = element_text(size=20))+labs(x="k",y="CV error",size=4)')
        ro.r('print(cvplot)')
        grdevices.dev_off()
        os.chdir(pwd)
    elif "admix" in type:
        ro.r('evec <- read.table(file="%s%s.evec")'%(pwd,basename))
        ro.r('am2 <- merge(am2,evec,all.x=TRUE,by.x="Sample",by.y="V1")')
        am2=ro.r['am2']
        grdevices.png(file="%s%s_MAP_AdmixGroup_%d_%sTC.png"%(admixoutdir,basename,lowestk,taxon_cov), width=1000, height=1000)
        ro.r('plot(ggmap(map)+ggtitle("Map by Admixture Group")+geom_point(data = am2,aes(LON,LAT,colour=factor(Sgroup),size=1),show.legend = FALSE)+scale_color_brewer(palette="Accent"))')
        grdevices.dev_off()
        ro.r('png("%s%s_PCA_AdmixGroup_%d_%sTC.png", width=1000, height=1000)'%(admixoutdir,basename,lowestk,taxon_cov))
        ro.r('myAPCAplot<-ggplot(am2,aes(x=V2,y=V3,color=factor(Sgroup)))+ggtitle("PCA by Admixture K=%s")+geom_point(size=4,show.legend=F)+scale_color_brewer(palette="Accent")+theme_dark(base_size=16)+labs(x="PCA axis 1",y="PCA axis 2")'%(lowestk))
        ro.r('print(myAPCAplot)')
        grdevices.dev_off()
    return


def admixture(pwd,base,k,datafile):
    print "\n\n\n\nstart admixture"
    os.chdir(pwd)
    cvd={}
    admixoutdir=pwd+"admixture/"
    if not os.path.exists(admixoutdir):
        os.mkdir(admixoutdir)
    shutil.copy(base+".map",admixoutdir)
    shutil.copy(base+".ped",admixoutdir)
    shutil.copy(base+"_outliersI.txt",admixoutdir)
    shutil.copy(datafile,admixoutdir)
    os.chdir(admixoutdir)
    subprocess.check_output(shlex.split("plink2 --threads 7 --file %s --make-bed --out %s"%(base,base)))
    shutil.copy("%s.bed"%(base),"%sI.bed"%(base))
    shutil.copy("%s.bim"%(base),"%sI.bim"%(base))

    ## make this a new function...?
    ## subprocess.check_output(shlex.split("plink2 --threads 7 --file %s --make-bed --remove %s_outliersI.txt --out %s"%(base,base,base)))
    for x in range(1,k):
        def admixer():
            admixcommand=shlex.split("admixture -j7 -C=0.01 --cv %s.bed %d"%(base,x))
            popen = subprocess.Popen(admixcommand, stdout=subprocess.PIPE, universal_newlines=True)
            with open(admixoutdir+"admixoutI.txt",'a') as admixoutfile:
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


def plotadmix(pwd,basename,lowestk,taxon_cov):
    admixoutdir=pwd+"admixture/"
    os.chdir(admixoutdir)
    taxdict={}
    qdict = {}
    qfile="%s%s.%d.Q"%(admixoutdir,basename,lowestk)
    Rplot(admixoutdir,basename,"CV",lowestk,taxon_cov)
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
    Rplot(pwd,basename,"admix",lowestk,taxon_cov)
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
    finpath=pwd+"Analysis/"
    if not os.path.exists(finpath):
        os.mkdir(finpath)
    for root, dirs, files in os.walk(pwd):
        for file in files:
            path=os.path.join(root,file)
            if file.endswith("png") and "figures" not in path:
                shutil.copy(path,figpath)
            elif "I." in path and "Analysis" not in path:
                shutil.copy(path,finpath)
    return

def controller(pwd,basename,k,datafile,taxon_cov):
    baseo='%s_o'%(basename)
    basep='%s_p'%(basename)
    pcafile='%s%s.evec'%(pwd,basename)
    ##Comment out steps to skip them. PCAer() and Rplot() need to be run together
    vcftoplink(pwd,basename,baseo,basep,taxon_cov)
    axes=PCAer(pwd,basep)
    Rplot(pwd,basep,axes,1,taxon_cov)
    lowestk=admixture(pwd,basep,k,datafile)
    plotadmix(pwd,basep,lowestk,taxon_cov)
    # lowestk=6
    # plotadmix(pwd,basep,lowestk,taxon_cov)
    # lowestk=2
    # plotadmix(pwd,basep,lowestk,taxon_cov)
    # lowestk=3
    # plotadmix(pwd,basep,lowestk,taxon_cov)
    # lowestk=4
    # plotadmix(pwd,basep,lowestk,taxon_cov)
    ro.r("write.table(am2, file='%smergeddataI.txt', quote=FALSE, row.names=FALSE, sep='\t')"%(pwd))
    cleanup(pwd)
    return

def main():
    # pwd=subprocess.check_output(shlex.split('pwd'))
    #basename = raw_input("Enter basename of VCF file:\n")
    pwd="/Users/josec/Desktop/marzo/marzo_all/noout_50min_outfiles/"
    basename="noout_50min"
    datafile="MarzoSet2.txt"
    # installtest(pwd)
    loaddata(pwd,datafile)
    k=10
    bs=10
    taxon_cov=0.90
    controller(pwd,basename,k+1,datafile,taxon_cov)
    taxon_cov=0.75
    controller(pwd,basename,k+1,datafile,taxon_cov)
    taxon_cov=0.50
    controller(pwd,basename,k+1,datafile,taxon_cov)
    # raxer(pwd,basename,bs)
    print "\n\nALLDone"
    return

if __name__ == '__main__':
    main()
