import os, sys, re, csv, statistics
import argparse, logging, warnings, json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

## Plots numbers in the second column of a two-column dataframe, presumably insert lengths
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check(expected_ext, expected_ext3, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext) or filename.lower().endswith(expected_ext3)):
                        raise ValueError()
                return openner(filename)
        return extension


logger = logging.getLogger("histogramPlotDataframe.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq (Expects single fastq input)', usage="python fullPlotShuffledFastq.py filepath/filename.fastq(.gz) --outputType [F/J/N]")

parser.add_argument('filename', type=ext_check('.txt', '.tsv', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='F', choices=['F', 'J', 'N'], help="--outputType F for full output,... J for .json only, and N for no image.")

parser.add_argument('--title', '-t', default='Norovirus MiSeq', help='Title for plot and output file names')

args = parser.parse_args()

basePath = os.getcwd()
inPath = os.path.split(args.filename.name)

logger.info("Reading file {}/{}".format(basePath, inPath[0]))

insertsFile = open(args.filename.name, "r")
insertsName = os.path.basename(args.filename.name)
insertsName = re.sub(r'\.(txt|tsv)', '', insertsName)

samples = []
inserts = []

set1 = pd.DataFrame( { insertsName+'_samples' : samples , insertsName+'_inserts' : inserts } )

iter = 0

for line in insertsFile:
        if(iter > 0):
                splitLine = line.split("\t")
                samples.append(splitLine[0])
                inserts.append(int(splitLine[1]))
        iter = iter + 1

logger.info("Finished reading {}".format(args.filename.name))

set1[insertsName+"_samples"] = pd.Series(samples)
set1[insertsName+"_inserts"] = pd.Series(inserts)

print(set1.describe())

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(6,5))

#fig.text(0.04, 0.5, 'Insert Lengths', va='center', rotation='vertical')

axes.hist(inserts, bins = 40, color='green')
axes.set_title("Inserts of " + args.title)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/' + args.title + '_inserts.png')


