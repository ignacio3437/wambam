#!/Users/josec/miniconda2/bin/python
import glob
import math
import os
import random
import re
import shlex
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.basemap import Basemap

"""
SNPViz, an automated workflow that:
- Reads a VCF file and associated metadata
- Prunes linked locci from the SNP dataset
- Removes locci with too much missing data
- Runs a PCA analysis to assess relatedness between samples
- Runs Admixture, which estimates the proportion of an individual's ancestries
- Picks a K value (number of ancestral populations) that best fits the
    data (i.e. lowest cross-validation error)
- Generates several figures including:
        - Map of sample locations and associated metadata (hypothetical groups)
        - Map of the samples and the estimated proportions of their ancestries
        - PCA plots colored by their metadata group and Admixture results
        - A CV plot showing the best value for K
"""


##############################Set Up##############################
pwd = "/Users/josec/Desktop/Testing/"  # Location of VCF file and Metadata File
basename = "Sy"  # Base name of VCF file. (Eg. So.vcf basename = "So")
datafile = "Sy_data.txt"  # Name of the Metadata.txt in TSV format
metagroup = "dummy"  # Name of column in Metadata.txt file to group sample by
k = 2  # Number of Ks to run the Admixture analysis
bs = 10  # Number of bootstraps for RaxML analysis
taxon_cov = 0.8  # If a locus percent missing data is below this number, it will be thrown out
threads = 7 # Number of cores to run the analysis
mapbuffer = 4  # Lat/Lon buffer to zoom out of box containing geographic distribution of samples
pretty_figures = False  # Set quality of output files 'high'=TRUE 'low'=FALSE
outlierbutton = ""  # Add "outliersigmathresh: 8" if you want outliers included, set to "" if you want outliers removed
##################################################################

# Set Global variables
datafile = pd.read_table(
    "%s%s" % (pwd, datafile), dtype={'Sample': object})
baseo = '%s_o' % (basename)
pcafile = '%s%s.evec' % (pwd, basename)
admixoutdir = pwd + "admixture/"
admix_colorpalette_list = sns.color_palette('Dark2', 8).as_hex()
colorpalette_list = sns.color_palette('Set1', 8).as_hex()
# Make high quality figures?
if pretty_figures:
    dpi, xpix, imgformat = 1200, 3000, '.svg'
else:
    dpi, xpix, imgformat = 300, 1000, '.png'


def VcfToPlink(baseo):
    # Changes from the VCF output from pyrad to a plink file
    # Linked locci are removed by passing --thin. This selects one SNP per locci
    # Assuming that locci are max 1000 bp long
    # Locci with missing data are removed with --max-missing
    subprocess.check_output(shlex.split("vcftools --vcf %s.vcf --thin 1000 --max-missing %f --out %s --plink" % (basename, taxon_cov, basename)))
    # Replace first column of map file to 1 to avoid errors concening too many Chromosomes in plink.
    # Plink treats each locci as a chromosome with the way the map file is encoded. This sets all locci on Chrom1
    # This does not affect anything because we are not doing any chromosomal analysis. Just analyzing a SNP dataset.
    subprocess.check_output(shlex.split("sed -ie 's/^0/1/' %s.map" % (basename)))
    return


def PCAer():
    # Sets up a pca.parameter file for eigenstrat and runs it, then sends it to be plotted in R.
    # Converts the eigenvalues to percent Variation for top axis. This axis will show percent variation for each PCA axis.
    with open(pwd + 'pca.par', 'w') as pca:
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
    numthreads:   %d
    lsqproject:    YES
    %s
        """ % (basename, basename, basename, basename, basename, threads,
               outlierbutton))
    pcaoutfile = pwd + basename + '_pcaoutI.txt'
    pcaout = subprocess.check_output(shlex.split("smartpca -p pca.par"))
    with open(pcaoutfile, 'w') as pcaoutfile:
        pcaoutfile.write(pcaout)
    # Pull out eigenvalues and P values from 4 lines starting with Tracy-Widom.
    lines = re.findall(r"Tracy-Widom(.*)\s.*\s(.*)\s(.*)", pcaout)[0]
    eigsum = lines[0].split()[2]
    eigens = (lines[1].split()[1], lines[2].split()[1])
    pvals = (lines[1].split()[4], lines[2].split()[4])
    # Calculate percent Variation using Eigenvalues and put as tuple to be used for plotting PCA later.
    ipervar = [(float(eig) / float(eigsum)) * 100 for eig in eigens]
    pervar = ["%.2f" % (i) for i in ipervar]
    # String tuple to be added to plots later
    pcaaxis = ("PCA1___PercentVar=" + pervar[0] + "_Pval=" + pvals[0],
               "PCA2___PercentVar=" + pervar[1] + "_Pval=" + pvals[1])
    # Find outliers and print to outlierfile
    if 'REMOVED outlier' in pcaout:
        outliers = re.findall(r'REMOVED outlier (\w*)', pcaout)
    else:
        outliers = []
    with open(pwd + basename + '_outliersI.txt', 'w') as outlierfile:
        for outlier in outliers:
            outlierfile.write(outlier + '\t' + outlier + '\n')
    print "There are %d outliers in PCA analysis of %s\n\n\n\n" % (
        len(outliers), basename)
    return pcaaxis


def PCApreper():
    # Store the PCA results in dataframe
    evectable = pd.read_table(
        "%s%s.evec" % (pwd, basename),
        skiprows=1,
        header=None,
        sep='\s*',
        names=["Sample", "Eigen1", "Eigen2", "Eigen3", "blank"],
        engine='python')
    # Combine PCA and metadata dataframes, merged by the sample name
    datatable = pd.merge(datafile, evectable, on='Sample', how='outer')
    # Remove entries with no metagroup assigned.
    datatable = datatable.dropna(subset=[metagroup], how='all')
    datatable.to_csv('%s%sM1.csv' % (pwd, basename))
    return


def PCAPlotter(df, colorby, colorpalette):
    # Plot PCA plot, colored by metagroup, axis are percent variation
    datatable = pd.read_csv('%s%s%s.csv' % (pwd, basename, df))
    datatable=datatable.sort_values(metagroup)
    sns.set_palette(colorpalette)
    pcapplot = sns.lmplot(
        'Eigen1',
        'Eigen2',
        markers="o",
        data=datatable,
        hue=colorby,
        fit_reg=False,
        legend=True,
        legend_out=False,
        scatter_kws={"linewidths": .1,
                     'edgecolors': 'black',
                     's': 25})
    pcapplot = (pcapplot.set_axis_labels(pcaaxis[0], pcaaxis[1]))
    plt.title('PCA plot', fontsize=8)
    pcapplot.savefig(
        '%s%s_PCA_%s_%s%s' % (pwd, basename, colorby, df, imgformat), dpi=dpi)
    plt.clf()
    return


def MapSetUp(datatable):
    # Function to download map given sampling locations and return as basemap object
    lats = datatable['LAT'].tolist()
    lons = datatable['LON'].tolist()
    # Bounding box for map
    lllat, lllon = (min(lats) - mapbuffer + 1), (min(lons) - mapbuffer)
    urlat, urlon = (max(lats) + mapbuffer - 1), (max(lons) + mapbuffer)
    # Generate the map with (australia=3577,mercator=3395) projection. Resolution set by pretty_figures variable
    # Flag this for further work. Why does epsg change the size of the map???? For now leave on 3577...
    m = Basemap(
        epsg=3577,
        llcrnrlat=lllat,
        urcrnrlat=urlat,
        llcrnrlon=lllon,
        urcrnrlon=urlon)
    m.arcgisimage(
        service='ESRI_Imagery_World_2D', xpixels=xpix, verbose=True, dpi=dpi)
    # If antarctic data use this. Not as nice....
    # m = Basemap(projection='spstere',boundinglat=-60,lon_0=45,resolution='h')
    # m.drawcoastlines(zorder=0,linewidth=0.5)
    return m


def SampleMapPlotter():
    # Plot Map with samples colored by metagroup
    datatable = pd.read_csv('%sM1.csv' % (basename))
    plt.title('Location of Samples', fontsize=12)
    m = MapSetUp(datatable)
    # Loop through each of the metagroups and plot points on map with different color
    # Sort to make the colors the same in all of the graphs
    datatable=datatable.sort_values(metagroup)
    metacategories = datatable[metagroup].tolist()
    noDupes = []
    [noDupes.append(i) for i in metacategories if not noDupes.count(i)]
    for i, mgroup in enumerate(noDupes):
        mtable = datatable.loc[datatable[metagroup] == mgroup]
        sublats = mtable['LAT'].tolist()
        sublons = mtable['LON'].tolist()
        x, y = m(sublons, sublats)
        m.scatter(
            x,
            y,
            linewidths=.1,
            edgecolors='black',
            c=colorpalette_list[i],
            marker='o',
            s=8)
    plt.savefig("%s%s_MAP_Pop%s" % (pwd, basename, imgformat), dpi=dpi)
    plt.clf()
    return


def AdmixtureRunner(k):
    # Runs Admixture program. Uses popout so that it will automatically terminate if admixutre
    # Hangs with "nan" at CV step. This then exits the program. Results are writen to admixout.txt
    # The CVD is a dictionary that stores each CV value to be called in AdmixtureSetUp.
    # CV is the cross validation error. The lowest CV value represents the best k.
    cvd = {}
    for x in range(1, k + 1):
        admixcommand = shlex.split("admixture -j%d -C=0.01 --cv=10 %s.bed %d" %
                                   (threads, basename, x))
        popen = subprocess.Popen(
            admixcommand, stdout=subprocess.PIPE, universal_newlines=True)
        with open(admixoutdir + "admixout.txt", 'a') as admixoutfile:
            admixoutfile.write(
                "####################K=%s\n####################\n" % (x))
            for stdout_line in iter(popen.stdout.readline, ""):
                admixoutfile.write(stdout_line)
                if "nan" in stdout_line:
                    return
                elif "CV" in stdout_line:
                    cv = re.findall(r'CV .*K=(\d*).*: (\d*\.\d*)', stdout_line)
                    cvd[x] = cv[0][1]
            popen.stdout.close()
            return_code = popen.wait()
            if return_code:
                raise subprocess.CalledProcessError(return_code, admixcommand)
        print "Admixture: Finished K=%d" % (x)
    return cvd


def AdmixtureSetUp(k):
    # Run admixture in new directory for range(1,k), keeps track of CV values to plot,
    # Set up the folder that admixture will run in and move input files to it
    print "\n\n\n\nstart admixture"
    os.chdir(pwd)
    # Make admixoutdir and add files to it
    if not os.path.exists(admixoutdir):
        os.mkdir(admixoutdir)
    for file in glob.glob(r"*%s*"%(basename)):
        shutil.copy(file, admixoutdir)
    os.chdir(admixoutdir)
    # Remove outliers from PCA from input files for Admixture analysis if present and make input files for Admixture.
    # Maybe make this a function in the future?
    outlierfile = basename + "_outliersI.txt"
    with open(outlierfile, 'r') as outlierhandle:
        outlierlist = outlierhandle.readlines()
    if len(outlierlist) > 0:
        subprocess.check_output(shlex.split("plink2 --file %s --remove %s --recode --make-bed --out %s" % (basename, outlierfile, basename)))
    else:
        subprocess.check_output(shlex.split("plink2 --file %s --recode --make-bed --out %s" % (basename, basename)))
    # Run admixture for each k 1->"k". cvd=dictionary for CV values for each K of admixture.
    cvd = AdmixtureRunner(k)
    # Pick lowest K based on CV error minimum
    lowestk = min(cvd, key=cvd.get)
    # Plot the CV results
    CVPloter(cvd)
    print "Best value for K is %s" % (lowestk)
    return lowestk


def draw_pie(ax, ratios, X, Y, size):
    # Helper function to plot a pie chart in a matplotlib.scatter() plot
    xy = []
    start = 0.0
    for ratio in ratios:
        x = [0] + np.cos(
            np.linspace(2 * math.pi * start, 2 * math.pi *
                        (start + ratio), 30)).tolist()
        y = [0] + np.sin(
            np.linspace(2 * math.pi * start, 2 * math.pi *
                        (start + ratio), 30)).tolist()
        xy1 = zip(x, y)
        xy.append(xy1)
        start += ratio
    for i, xyi in enumerate(xy):
        ax.scatter(
            [X], [Y],
            marker=(xyi, 0),
            s=size,
            facecolor=admix_colorpalette_list[i],
            alpha=.6,
            edgecolors='none')
    return


def MapAdmix(lowestk):
    # Plot the admix results as piecharts and put each sample on a map
    datadfM2 = pd.read_csv('%s%sM2.csv' % (pwd, basename), index_col=0)
    # Set up Map
    plt.title('Admixture Map k=%s' % (lowestk), fontsize=12)
    m = MapSetUp(datadfM2)
    ax = plt.subplot()
    # Plot each sample as piechard of admixture results on map.
    for i, sample in enumerate(datadfM2.index.tolist()):
        slat, slon = datadfM2.loc[[sample]]['LAT'].values[0], datadfM2.loc[[
            sample
        ]]['LON'].values[0]
        admix_result_list = datadfM2.loc[[sample]][[
            str(xx) for xx in range(lowestk)
        ]].values.tolist()[0]
        X, Y = m(slon, slat)
        draw_pie(ax, admix_result_list, X, Y, size=50)
        # Assign an each sample to a k population and save in M3_lowestk df
        admixgroup = str(admix_result_list.index(max(admix_result_list)) + 1)
        datadfM2.at[sample, "AdmixGroup"] = admixgroup
    plt.savefig(
        "%s%s_MAP_Admix%d%s" % (pwd, basename, lowestk, imgformat), dpi=dpi)
    plt.clf()
    datadfM2.to_csv('%s%sM3_%d.csv' % (pwd, basename, lowestk))
    return


def CVPloter(cvd):
    # Plot the admixture CV error results as a lineplot.
    plt.title('CV error for different K values', fontsize=8)
    xlist = [int(x) for x in cvd.keys()]
    ylist = [float(y) for y in cvd.values()]
    sns.pointplot(xlist, ylist)
    plt.savefig("%s%s_CVplot%s" % (pwd, basename, imgformat), dpi=dpi)
    plt.clf()
    return


def PlotAdmix(lowestk):
    # Plots admixture results as a map, 'structure plot', and CV plot.
    os.chdir(admixoutdir)
    # File Paths
    qfile = '%s%s.%d.Q' % (admixoutdir, basename, lowestk)
    sampleorderfile = '%s%s.fam' % (admixoutdir, basename)
    datafile = '%s%sM1.csv' % (pwd, basename)
    # Read the Plink file to get the sample names to assign to the admixture analysis
    with open(sampleorderfile, 'r') as sampleorder_handle:
        samplelist = [x.split(' ')[0] for x in sampleorder_handle.readlines()]
    # Combine with metadata into one table joining on the sample index
    qdf = pd.read_table(qfile, header=None, sep='\s', engine='python')
    ddf = pd.read_csv(datafile, index_col=1)
    ddf=ddf.sort_values(metagroup)
    qdf.index = samplelist
    # Set the index as a string to make sure they join correctly.Then save
    qdf.index = qdf.index.map(str)
    ddf.index = ddf.index.map(str)
    datadfM2 = ddf.join(qdf, how='inner')
    datadfM2.to_csv('%s%sM2.csv' % (pwd, basename))
    # Plot the admix results as piecharts and put each sample on a map
    MapAdmix(lowestk)
    # Plot the PCA colored by Admix results
    PCAPlotter('M3_%d' % (lowestk), 'AdmixGroup', admix_colorpalette_list)
    return


def Raxer(pwd, basename, bs):
    # Print subprocess.check_output(shlex.split("raxmlHPC-PTHREADS-SSE3 -T 8 -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d"%(basename,basename,random.randint(0,999999),bs,random.randint(0,999999))))
    print "raxmlHPC-PTHREADS-SSE3 -T %d -f a -n %s -s %s.phy -x %d -N %d -m GTRCAT -p %d" % (
        threads, basename, basename, random.randint(0, 999999), bs,
        random.randint(0, 999999))
    print "raxmlHPC-PTHREADS-SSE3 -T %d -f x -n %s -s %s.phy -m GTRCAT -p %d" % (
        threads, basename, basename, random.randint(0, 999999))
    return


def DirectoryCleaner(pwd):
    # Moves important files to Analysis directory and figures to figure Directory.
    figpath = pwd + "figures/"
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    for root, dirs, files in os.walk(pwd):
        for file in files:
            path = os.path.join(root, file)
            if file.endswith(imgformat) and "figures" not in path:
                shutil.copy(path, figpath)
    return


def Controller(k):
    # Engine for the workflow
    ################ Comment out steps to skip them.################
    os.chdir(pwd)
    VcfToPlink(baseo)
    global pcaaxis
    pcaaxis = PCAer()
    PCApreper()
    PCAPlotter('M1', metagroup, colorpalette_list)
    SampleMapPlotter()
    lowestk = AdmixtureSetUp(k)
    # PlotAdmix(lowestk)
    # Plot with lowestk then plot forcing k=2
    lowestk = 2
    PlotAdmix(lowestk)
    # Raxer(pwd,basename,bs)
    DirectoryCleaner(pwd)
    return


def main():
    Controller(k)
    print "\n\nAll done!"
    return


if __name__ == '__main__':
    main()
