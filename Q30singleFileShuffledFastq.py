import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tabulate import tabulate

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

## Check input file for zero-length
def check_for_empty(fastq_file):
        path1 = Path(fastq_file)
        if(path1.stat().st_size == 0):
            raise ArgumentTypeError(f'{fastq_file} cannot be empty')


def readable_dir(prospective_dir):
        if not os.path.isdir(prospective_dir):
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
                if( not prospective_dir.endswith("/") ):
                        prospective_dir = prospective_dir + "/"
                return prospective_dir
        else:
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

origWD = os.getcwd()


## Two shuffled, paired-end .fastq files expected as input
parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq', usage="Q30singleFileShuffledFastq.py filepath/filename1.fastq")

parser.add_argument('filename1', type=ext_check('.fastq', argparse.FileType('r')), nargs='+')

## outputType enables suppression of dataframe (.csv) output files or suppression of histogram (.png) output files
parser.add_argument('--outputType', '-o', default='F', choices=['F', 'P', 'C', 'Q'], help="--outputType F for full output (plots and .csv), P for plots only, C for csv file with no plot, and Q for Q30 STDOUT summary only.")

parser.add_argument('--paired', '-u', default='F', choices=['T', 'F'], help="--paired F to calculate a average PHRED for all reads together or --paired T for input file of forward and reverse shuffled together.")

## output folder
parser.add_argument('--outDir', '-D', type=readable_dir, required=True, action='store')

args = parser.parse_args()

outFolder = args.outDir

outFilePath = origWD + '/' + outFolder + '/'

## File names used in plot titles
myTitle1 = re.split(r'[\.\/]', args.filename1[0].name)

csvRow1 = []
forwardName = []
reverseName = []
allName = []
allAvg = []
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
                allName.append(record.id)
                j = 0
                r1Len_1 = r1Len_1 + len(record.seq)
                while( j < len(record.seq)):
                        if(record.letter_annotations["phred_quality"][j] >= 30):
                                r1Q30_1 = r1Q30_1 + 1
                        j = j + 1
                forwardAvg1.append(statistics.mean(record.letter_annotations["phred_quality"]))
                allAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                reverseName.append(record.id)
                allName.append(record.id)
                strand = record.description.split(" ")
                j = 0
                r2Len_1 = r2Len_1 + len(record.seq)
                while( j < len(record.seq) ):
                        if(record.letter_annotations["phred_quality"][j] >= 30 ):
                                r2Q30_1 = r2Q30_1 + 1
                        j = j + 1
                reverseAvg1.append(statistics.mean(record.letter_annotations["phred_quality"]))
                allAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
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
        fig1.savefig(outFilePath + 'fwd_and_rev_PHRED.png')

if((args.outputType != 'P') and (args.paired == 'T')):

        dfMiSeqPHRED = pd.DataFrame()
        stringList = []
        stringFwdList = []
        stringRevList = []
        for i in forwardAvg1:
            stringFwdList.append(str(round(float(i), 2)))
                
        for j in reverseAvg1:
            stringRevList.append(str(round(float(j), 2)))
            
        pairedMiSeqPHRED = { "R1_Read_ID" : forwardName,  "R1_PHRED" : stringFwdList, "R2_PHRED" : stringRevList }
        dfMiSeqPHRED = pd.DataFrame(pairedMiSeqPHRED)
            
        dfMiSeqPHRED.to_csv(outFilePath + 'fwd_and_rev_PHRED.csv', index=False)
            
if((args.outputType != 'P') and (args.paired == 'F')):
        dfMiSeqPHRED = pd.DataFrame()
        stringList = []
        
        for i in allAvg:
            stringList.append(str(round(float(i), 2)))
        
        singleMiSeqPHRED = {"Read_ID" : allName, "Read_PHRED" : stringList }
        dfMiSeqPHRED = pd.DataFrame(singleMiSeqPHRED)

        dfMiSeqPHRED.to_csv(outFilePath + 'all_PHRED.csv', index=False)



