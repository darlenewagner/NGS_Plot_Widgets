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
parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq', usage="Q30singleFileShuffledFastq.py filepath/filename1.fastq")

parser.add_argument('filename1', type=ext_check('.fastq', argparse.FileType('r')), nargs='+')

## outputType enables suppression of dataframe (.csv) output files or suppression of histogram (.png) output files
parser.add_argument('--outputType', '-o', default='F', choices=['F', 'P', 'C', 'Q'], help="--outputType F for full output (plots and .csv), P for plots only, C for csv file with no plot, and Q for Q30 STDOUT summary only.")

parser.add_argument('--unpaired', '-u', default='F', choices=['T', 'F'], help="--unpaired F to calculate a proportion for forward and reverse reads separately and --unpaired T for combined forward and reverse or a file of single strand reads.")

args = parser.parse_args()

#print(args.filename1[0].name)
#print(args.filename2[0].name)

## File names used in plot titles
myTitle1 = re.split(r'[\.\/]', args.filename1[0].name)

csvRow1 = []
forwardName = []
reverseName = []
forwardLen1 = []
reverseLen1 = []
forwardAvg1 = []
reverseAvg1 = []

iter = 0

myFastq1 = open(args.filename1[0].name, "r")

r1Q30_1 = 0
r1Len_1 = 1
r2Q30_1 = 0
r2Len_1 = 1

for record in SeqIO.parse(myFastq1, "fastq"):
        if(iter % 2 == 0):
                forwardName.append(record.id)
                j = 0
                r1Len_1 = r1Len_1 + len(record.seq)
                while( j < len(record.seq)):
                        if(record.letter_annotations["phred_quality"][j] >= 30):
                                r1Q30_1 = r1Q30_1 + 1
                        j = j + 1
                forwardAvg1.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                reverseName.append(record.id)
                strand = record.description.split(" ")
                j = 0
                r2Len_1 = r2Len_1 + len(record.seq)
                while( j < len(record.seq) ):
                        if(record.letter_annotations["phred_quality"][j] >= 30 ):
                                r2Q30_1 = r2Q30_1 + 1
                        j = j + 1
                reverseAvg1.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

if(args.unpaired == 'F'):
        print("%s, Forward_Q30%%: %2.2f, Reverse_Q30%%: %2.2f" % (myTitle1[len(myTitle1) - 2], 100*r1Q30_1/r1Len_1, 100*r2Q30_1/r2Len_1))
else:
        print("%s, Paired_Q30%%: %2.2f" % (myTitle1[len(myTitle1) - 2], 100*(r1Q30_1 + r2Q30_1)/(r1Len_1 + r2Len_1)))


if(args.outputType == 'Q'):
        sys.exit()


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
        axes1[0].set_title("R1 " + myTitle1[len(myTitle1) - 2], fontsize = BIG_SIZE)
        axes1[0].set(ylabel='Read Counts')

        ## Plot PHRED Quality for R2 reads as 1D histogram
        axes1[1].xaxis.label.set_size(MEDIUM_SIZE)
        axes1[1].yaxis.label.set_size(MEDIUM_SIZE)
        axes1[1].tick_params(axis='x', labelsize=SMALL_SIZE)
        axes1[1].tick_params(axis='y', labelsize=SMALL_SIZE)
        axes1[1].hist(forwardAvg2, bins = 40, color='blue')
        axes1[1].set_title("R2 " + myTitle2[len(myTitle2) - 2], fontsize = BIG_SIZE)
        axes1[1].set(ylabel='Read Counts')
        axes1[1].set(xlabel='Average Read Quality')
        fig1.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/fwd_and_rev_PHRED.png')

if(args.outputType != 'P'):
        dfMiSeqPHRED = pd.DataFrame()
        if(args.unpaired == 'F'):
            pairedMiSeqPHRED = { "R1_Read_ID" : forwardName,  "R1_MiSeq_PHRED" : forwardAvg1, "R2_Read_ID" : reverseName, "R2_MiSeq_PHRED" : reverseAvg1 }
            dfMiSeqPHRED = pd.DataFrame(pairedMiSeqPHRED)
        else:
            forwardName.append(reverseName)
            forwardAvg1.append(reverseAvg1)
            singleMiSeqPHRED = { "R1_Read_ID" : forwardName,  "R1_MiSeq_PHRED" : forwardAvg1, "R2_Read_ID" : reverseName, "R2_MiSeq_PHRED" : reverseAvg1 }
            dfMiSeqPHRED = pd.DataFrame(singleMiSeqPHRED)
        
        dfMiSeqPHRED.to_csv('/scicomp/home-pure/ydn3/NGS_Plot_Widgets/fwd_and_rev_PHRED.csv')   

