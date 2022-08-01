import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

## Script for plotting average PHRED score per read and outputting data frames of PHRED averages
## Requires Biopython, Numpy, and Matplotlib
## Required Input: Two shuffled, paired-end .fastq files with equal numbers of forward (R1) and reverse (R2) reads
## Output: Two .png files and/or two .csv files
## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.lower().endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

## Two shuffled, paired-end .fastq files expected as input

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq', usage="plotDoubleShuffledFastq.py filepath/filename1.fastq filepath/filename2.fastq")

parser.add_argument('filename1', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))
parser.add_argument('filename2', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

## outputType enables suppression of dataframe (.csv) output files or suppression of histogram (.png) output files
parser.add_argument('--outputType', '-o', default='F', choices=['F', 'P', 'C'], help="--outputType F for full output (plots and .csv), P for plots only, and C for csv file with no plot.")

args = parser.parse_args()

print(args.filename1[0].name)
print(args.filename2[0].name)

## File names used in plot titles
myTitle1 = re.split(r'[\.\/]', args.filename1[0].name)
myTitle2 = re.split(r'[\.\/]', args.filename2[0].name)

print(myTitle1[len(myTitle1) - 2])
print(myTitle2[len(myTitle2) - 2])

csvRow1 = []
forwardLen1 = []
reverseLen1 = []
forwardAvg1 = []
reverseAvg1 = []

csvRow2 = []
forwardLen2 = []
reverseLen2 = []
forwardAvg2 = []
reverseAvg2 = []

iter = 0

myFastq1 = open(args.filename1[0].name, "r")
myFastq2 = open(args.filename2[0].name, "r")

r1Q30_1 = 0
r1Len_1 = 1
r2Q30_1 = 0
r2Len_1 = 1

for record in SeqIO.parse(myFastq1, "fastq"):
        if(iter % 2 == 0):
                j = 0
                r1Len_1 = r1Len_1 + len(record.seq)
                while( j < len(record.seq)):
                        if(record.letter_annotations["phred_quality"][j] >= 30):
                                r1Q30_1 = r1Q30_1 + 1
                        j = j + 1
                forwardAvg1.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                strand = record.description.split(" ")
                j = 0
                r2Len_1 = r2Len_1 + len(record.seq)
                while( j < len(record.seq) ):
                        if(record.letter_annotations["phred_quality"][j] >= 30 ):
                                r2Q30_1 = r2Q30_1 + 1
                        j = j + 1
                reverseAvg1.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1


print("%s, Forward: %0.2f, Reverse: %0.2f" % (myTitle1[len(myTitle1) - 2], r1Q30_1/r1Len_1, r2Q30_1/r2Len_1))


r1Q30_2 = 0
r1Len_2 = 1
r2Q30_2 = 0
r2Len_2 = 1

for record in SeqIO.parse(myFastq2, "fastq"):
        if(iter % 2 == 0):
                j = 0
                r1Len_2 = r1Len_2 + len(record.seq)
                while( j < len(record.seq)):
                        if(record.letter_annotations["phred_quality"][j] >= 30):
                                r1Q30_2 = r1Q30_2 + 1
                        j = j + 1
                forwardAvg2.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                strand = record.description.split(" ")
                j = 0
                r2Len_2 = r2Len_2 + len(record.seq)
                while( j < len(record.seq) ):
                        if(record.letter_annotations["phred_quality"][j] >= 30 ):
                                r2Q30_2 = r2Q30_2 + 1
                        j = j + 1
                reverseAvg2.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print("%s, Forward: %0.2f, Reverse: %0.2f" % (myTitle2[len(myTitle2) - 2], r1Q30_2/r1Len_2, r2Q30_2/r2Len_2))

if(args.outputType != 'C'):
        SMALL_SIZE = 20
        MEDIUM_SIZE = 24
        BIG_SIZE = 30
        fig1, axes1 = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(14,20))

        ## Plot PHRED Quality for R1 reads as 1D histogram
        axes1[0].xaxis.label.set_size(MEDIUM_SIZE)
        axes1[0].yaxis.label.set_size(MEDIUM_SIZE)
        axes1[0].tick_params(axis='x', labelsize=SMALL_SIZE)
        axes1[0].tick_params(axis='y', labelsize=SMALL_SIZE)
        axes1[0].hist(forwardAvg1, bins = 40, color='blue')
        axes1[0].set_title("R1 MiSeq " + myTitle1[len(myTitle1) - 2], fontsize = BIG_SIZE)
        axes1[0].set(ylabel='Read Counts')

        ## Plot PHRED Quality for R2 reads as 1D histogram
        axes1[1].xaxis.label.set_size(MEDIUM_SIZE)
        axes1[1].yaxis.label.set_size(MEDIUM_SIZE)
        axes1[1].tick_params(axis='x', labelsize=SMALL_SIZE)
        axes1[1].tick_params(axis='y', labelsize=SMALL_SIZE)
        axes1[1].hist(forwardAvg2, bins = 40, color='blue')
        axes1[1].set_title("R1 iSeq " + myTitle2[len(myTitle2) - 2], fontsize = BIG_SIZE)
        axes1[1].set(ylabel='Read Counts')
        axes1[1].set(xlabel='Average Read Quality')
        fig1.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/fwd_miSeq_and_iSeq_PHRED.png')
        fig2, axes2 = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(14,20))

        ## Plot PHRED Quality for R1 reads as 1D histogram
        axes2[0].xaxis.label.set_size(MEDIUM_SIZE)
        axes2[0].yaxis.label.set_size(MEDIUM_SIZE)
        axes2[0].tick_params(axis='x', labelsize=SMALL_SIZE)
        axes2[0].tick_params(axis='y', labelsize=SMALL_SIZE)
        axes2[0].hist(reverseAvg1, bins = 40, color='red')
        axes2[0].set_title("R2 MiSeq " + myTitle1[len(myTitle1) - 2], fontsize = BIG_SIZE)
        axes2[0].set(ylabel='Read Counts')

        ## Plot PHRED Quality for R2 reads as 1D histogram
        axes2[1].xaxis.label.set_size(MEDIUM_SIZE)
        axes2[1].yaxis.label.set_size(MEDIUM_SIZE)
        axes2[1].tick_params(axis='x', labelsize=SMALL_SIZE)
        axes2[1].tick_params(axis='y', labelsize=SMALL_SIZE)
        axes2[1].hist(reverseAvg2, bins = 40, color='red')
        axes2[1].set_title("R2 iSeq " + myTitle2[len(myTitle2) - 2], fontsize = BIG_SIZE)
        axes2[1].set(ylabel='Read Counts')
        axes2[1].set(xlabel='Average Read Quality')
        fig2.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/rev_miSeq_and_iSeq_PHRED.png')

if(args.outputType != 'P'):
        totalMiSeqPHRED = { "R1_MiSeq_PHRED" : forwardAvg1, "R2_MiSeq_PHRED" : reverseAvg1 }
        totaliSeqPHRED = { "R1_iSeq_PHRED" : forwardAvg2, "R2_iSeq_PHRED" : reverseAvg2 }
        dfMiSeqPHRED = pd.DataFrame(totalMiSeqPHRED)
        dfiSeqPHRED = pd.DataFrame(totaliSeqPHRED)
        dfMiSeqPHRED.to_csv('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/fwd_rev_miSeq_PHRED.csv')   
        dfiSeqPHRED.to_csv('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/fwd_rev_iSeq_PHRED.csv')   
