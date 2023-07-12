import os, sys, re, csv, statistics
import argparse, logging, warnings, json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

## Plots numbers in the second column of a two-column dataframe, presumably insert lengths
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check(expected_ext, expected_ext3, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext) or filename.lower().endswith(expected_ext3)):
                        raise ValueError()
                return openner(filename)
        return extension

# Function: Intersection (default)

def Intersection(list1, list2):
    return sorted(set(list1).intersection(list2))

## Function: Complement

def Compl(list1, intersec):
    snpSet1 = set(list1)
    snpSet2 = set(intersec)
    return list(snpSet1 - snpSet2)


logger = logging.getLogger("vennDiagramPlotColumn.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='Computes intersection and discordant SNPs from two input lists of integers', usage="python vennDiagramPlotColumn.py filepath/snpsPositions1.tsv filepath/snpsPositions2.tsv --outputType [S|P]")

parser.add_argument('filename1', type=ext_check('.txt', '.tsv', argparse.FileType('r')))

parser.add_argument('filename2', type=ext_check('.txt', '.tsv', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='S', choices=['S', 'P'], help="--outputType S for simple, text output or P for venn diagram plot image.")

parser.add_argument('--title', '-t', default='SNPs Concordance', help='Title for plot and output file names')

args = parser.parse_args()

basePath = os.getcwd()
inPath = os.path.split(args.filename1.name)

logger.info("Reading file {}/{}".format(basePath, inPath[0]))

snpsFile1 = open(args.filename1.name, "r")
snpsFile2 = open(args.filename2.name, "r")
snpsName1 = os.path.basename(args.filename1.name)
snpsName1 = re.sub(r'\.(txt|tsv)', '', snpsName1)
snpsName2 = os.path.basename(args.filename2.name)
snpsName2 = re.sub(r'\.(txt|tsv)', '', snpsName2)

samples1 = []
positions1 = []
samples2 = []
positions2 = []

set1 = pd.DataFrame( { snpsName1+'_samples' : samples1 , snpsName1+'_SNPs' : positions1 } )
set2 = pd.DataFrame( { snpsName2+'_samples' : samples2 , snpsName2+'_SNPs' : positions2 } )

iter = 0

for line in snpsFile1:
        if(iter > 0):
                splitLine = line.split("\t")
                samples1.append(splitLine[0])
                positions1.append(int(splitLine[1]))
        iter = iter + 1

logger.info("Finished reading {}".format(args.filename1.name))

snpsFile1.close()

set1[snpsName1+"_samples"] = pd.Series(samples1)
set1[snpsName1+"_SNPs"] = pd.Series(positions1)

print(set1.describe())

iter = 0

for line in snpsFile2:
        if(iter > 0):
                splitLine = line.split("\t")
                samples2.append(splitLine[0])
                positions2.append(int(splitLine[1]))
        iter = iter + 1

logger.info("Finished reading {}".format(args.filename2.name))

snpsFile2.close()

set2[snpsName2+"_samples"] = pd.Series(samples2)
set2[snpsName2+"_SNPs"] = pd.Series(positions2)

print(set2.describe())

concordant = Intersection(positions1, positions2)
discordant1 = Compl(positions1, positions2)
discordant2 = Compl(positions2, positions1)

#fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(6,5))

#fig.text(0.04, 0.5, 'SNPs Count', va='center', rotation='vertical')
#fig.text(0.45, 0.04, 'Genome Position', va='center')

#axes.hist(positions1, bins = 30, color='green')
#axes.set_title("SARS-CoV-2 " + args.title)

venn2(subsets = (len(discordant1), len(discordant2), len(concordant)), set_labels = ('MiSeq', 'iSeq'))
plt.title("MiSeq to iSeq SNPs Concordance")

## remember to change filepath to your local installation of the Python virtual environment
plt.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/' + args.title + '_allSNPs_positions.png')


