import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
from matplotlib import pyplot as plt

## Script for calculating sequence length and average PHRED score per read
## Outputs comma-delimited table to stdout
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.lower().endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq', usage="readShuffledFastq.py filepath/filename.fastq")

parser.add_argument('filename', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

args = parser.parse_args()

print(args.filename[0].name)

myTitle = re.split(r'\/', args.filename[0].name)
newTitle = re.sub('\.fastq', '', myTitle[len(myTitle) - 1])
print(newTitle)

csvRow = []
forwardLen = []
reverseLen = []
forwardAvg = []
reverseAvg = []

iter = 0

myFastq = open(args.filename[0].name, "r")

for record in SeqIO.parse(myFastq, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])), end="")
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

medLenFwd = statistics.median(forwardLen)
avgLenFwd = statistics.mean(forwardLen)
minLenFwd = min(forwardLen)
medFwdAvgPHRED = statistics.median(forwardAvg)
minFwdAvgPHRED = min(forwardAvg)
medLenRev = statistics.median(reverseLen)
avgLenRev = statistics.mean(reverseLen)
minLenRev = min(reverseLen)
medRevAvgPHRED = statistics.median(reverseAvg)
minRevAvgPHRED = min(reverseAvg)

#forwardAvg = np.random.normal(size = 500000) + 30
#reverseAvg = 1.5*np.random.normal(size = 500000) + 25

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))
fig.text(0.5, 0.04, 'Average Read Quality', ha='center')
fig.text(0.04, 0.5, 'Read Counts', va='center', rotation='vertical')

## Plot PHRED Quality for R1 reads as 1D histogram
axes[0].hist(forwardAvg, bins = 25, color='blue')
axes[0].set_title("R1 Quality " + newTitle)
#axes[0].xlabel('Average Read Quality')
#axes[0].ylabel('Read Counts')
#axes[0].savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/test_forwardPHRED.png')

## Plot PHRED Quality for R2 reads as 1D histogram
axes[1].hist(reverseAvg, bins = 25, color='red')
axes[1].set_title("R2 Quality " + newTitle)
#axes[1].xlabel('Average Read Quality')
#axes[1].ylabel('Read Counts')
#axes[1].savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/test_reversePHRED.png')


os.mkdir('QC_' + newTitle)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/fwd_and_revPHRED.png')

summary = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/ReadStatistics.README.txt', 'w')

summary.write("Summary Statistics for " + newTitle + "\n")
summary.write("R1 median length:\t0.2f" % (medLenFwd) + ",\tR2 median length:\t0.2f\n" % (medLenRev))
summary.write("R1 mean length:\t0.2f" % (avgLenFwd) + ",\tR2 mean length:\t0.2f\n" % (avgLenRev))
summary.write("R1 minimum length:\t" + str(minLenFwd) + ",\tR2 minimum length:\t\n" + str(minLenRev))
summary.write("R1 median PHRED:\t0.2f" % (medFwdAvgPHRED) + ",\tR2 median PHRED:\t0.2f\n" % (medRevAvgPHRED))
summary.write("R1 minimum PHRED:\t" + str(minFwdAvgPHRED) + ",\tR2 minimum PHRED:\t\n" + str(minRevAvgPHRED))

