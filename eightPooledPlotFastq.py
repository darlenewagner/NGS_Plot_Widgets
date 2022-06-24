import os, sys, re, csv, statistics
import argparse, logging, warnings, json
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

## Script for calculating sequence length and average PHRED score per read
## Outputs comma-delimited table to stdout
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check(expected_ext, expected_ext3, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext) or filename.lower().endswith(expected_ext3)):
                        raise ValueError()
                return openner(filename)
        return extension

logger = logging.getLogger("eightPooledPlotFastq.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq (Expects single fastq input)', usage="python eightPooledPlotFastq.py filepath/filename1.fastq(.gz) filepath/filename2.fastq(.gz) filepath/filename3.fastq(.gz) filepath/filename4.fastq(.gz) filepath/filename5.fastq(.gz) filepath/filename6.fastq(.gz) filepath/filename7.fastq(.gz) --outputType [F/J/N]")

parser.add_argument('filename1', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename2', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename3', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename4', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename5', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename6', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename7', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))
parser.add_argument('filename8', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='F', choices=['F', 'C', 'N'], help="--outputType F for full output including plot of total PHRED scores, C for pooled output and individual summary files with no plot, and N for pooled output only, no image.")

parser.add_argument('--title', '-t', default='SARS-CoV-2', help="Title for plot and output file names, if applicable")

args = parser.parse_args()

#print(args.filename1.name)print(args.filename2.name)print(args.filename3.name)print(args.filename4.name)print(args.filename5.name)print(args.filename6.name)print(args.filename7.name)

basePath = os.getcwd()
inPath1 = os.path.split(args.filename1.name)
inPath2 = os.path.split(args.filename2.name)
inPath3 = os.path.split(args.filename3.name)
inPath4 = os.path.split(args.filename4.name)
inPath5 = os.path.split(args.filename5.name)
inPath6 = os.path.split(args.filename6.name)
inPath7 = os.path.split(args.filename7.name)
inPath8 = os.path.split(args.filename8.name)

myTitle = re.split(r'\/', args.filename1.name)
newTitle = re.sub('\.f(ast)?q(\.gz)?', '', myTitle[len(myTitle) - 1])
#print(newTitle)

csvRow = []
forwardLen = []
reverseLen = []
forwardAvg = []
reverseAvg = []


nameCheck = ''
pairingPreserved = 'Y'
warned = 0

logger.info("Reading file {}/{}".format(basePath, inPath1[0]))
logger.info("Reading file {}/{}".format(basePath, inPath2[0]))
logger.info("Reading file {}/{}".format(basePath, inPath3[0]))
logger.info("Reading file {}/{}".format(basePath, inPath4[0]))
logger.info("Reading file {}/{}".format(basePath, inPath5[0]))
logger.info("Reading file {}/{}".format(basePath, inPath6[0]))
logger.info("Reading file {}/{}".format(basePath, inPath7[0]))
logger.info("Reading file {}/{}".format(basePath, inPath8[0]))

myFastq1 = ''
myFastq2 = ''
myFastq3 = ''
myFastq4 = ''
myFastq5 = ''
myFastq6 = ''
myFastq7 = ''
myFastq8 = ''

tempFas1 = re.sub(r'.gz', '', args.filename1.name)

if(args.filename1.name.endswith('.gz') == True ):
        os.system("gunzip -c {} > {}".format(args.filename1.name, tempFas))
        myFastq1 = open(tempFas1, "r")
else:
        myFastq1 = open(args.filename1.name, "r")
        fastqName1 = os.path.basename(args.filename1.name)
        fastqName1 = re.sub(r'\.f(ast)?q', '', fastqName1)
        logger.info("Read file {}".format(fastqName1))
        myFastq2 = open(args.filename2.name, "r")
        fastqName2 = os.path.basename(args.filename2.name)
        fastqName2 = re.sub(r'\.f(ast)?q', '', fastqName2)
        logger.info("Read file {}".format(fastqName2))
        myFastq3 = open(args.filename3.name, "r")
        fastqName3 = os.path.basename(args.filename3.name)
        fastqName3 = re.sub(r'\.f(ast)?q', '', fastqName3)
        logger.info("Read file {}".format(fastqName3))
        myFastq4 = open(args.filename4.name, "r")
        fastqName4 = os.path.basename(args.filename4.name)
        fastqName4 = re.sub(r'\.f(ast)?q', '', fastqName4)
        logger.info("Read file {}".format(fastqName4))
        myFastq5 = open(args.filename5.name, "r")
        fastqName5 = os.path.basename(args.filename5.name)
        fastqName5 = re.sub(r'\.f(ast)?q', '', fastqName5)
        logger.info("Read file {}".format(fastqName5))
        myFastq6 = open(args.filename6.name, "r")
        fastqName6 = os.path.basename(args.filename6.name)
        fastqName6 = re.sub(r'\.f(ast)?q', '', fastqName6)
        logger.info("Read file {}".format(fastqName6))
        myFastq7 = open(args.filename7.name, "r")
        fastqName7 = os.path.basename(args.filename7.name)
        fastqName7 = re.sub(r'\.f(ast)?q', '', fastqName7)
        logger.info("Read file {}".format(fastqName7))
        myFastq8 = open(args.filename8.name, "r")
        fastqName8 = os.path.basename(args.filename8.name)
        fastqName8 = re.sub(r'\.f(ast)?q', '', fastqName8)
        logger.info("Read file {}".format(fastqName8))

newTitle = args.title

set1 = pd.DataFrame( { fastqName1+'_R1_Len' : forwardLen, fastqName1+'_R2_Len' : reverseLen, fastqName1+'_R1_Qual' : forwardAvg, fastqName1+'_R2_Qual' : reverseAvg } )
set2 = pd.DataFrame( { fastqName2+"_R1_Len" : forwardLen, fastqName2+"_R2_Len" : reverseLen, fastqName2+"_R1_Qual" : forwardAvg, fastqName2+"_R2_Qual" : reverseAvg } )
set3 = pd.DataFrame( { fastqName3+"_R1_Len" : forwardLen, fastqName3+"_R2_Len" : reverseLen, fastqName3+"_R1_Qual" : forwardAvg, fastqName3+"_R2_Qual" : reverseAvg } )
set4 = pd.DataFrame( { fastqName4+"_R1_Len" : forwardLen, fastqName4+"_R2_Len" : reverseLen, fastqName4+"_R1_Qual" : forwardAvg, fastqName4+"_R2_Qual" : reverseAvg } )
set5 = pd.DataFrame( { fastqName5+"_R1_Len" : forwardLen, fastqName5+"_R2_Len" : reverseLen, fastqName5+"_R1_Qual" : forwardAvg, fastqName5+"_R2_Qual" : reverseAvg } )
set6 = pd.DataFrame( { fastqName6+"_R1_Len" : forwardLen, fastqName6+"_R2_Len" : reverseLen, fastqName6+"_R1_Qual" : forwardAvg, fastqName6+"_R2_Qual" : reverseAvg } )
set7 = pd.DataFrame( { fastqName7+"_R1_Len" : forwardLen, fastqName7+"_R2_Len" : reverseLen, fastqName7+"_R1_Qual" : forwardAvg, fastqName7+"_R2_Qual" : reverseAvg } )
set8 = pd.DataFrame( { fastqName8+"_R1_Len" : forwardLen, fastqName8+"_R2_Len" : reverseLen, fastqName8+"_R1_Qual" : forwardAvg, fastqName8+"_R2_Qual" : reverseAvg } )

iter = 0

for record in SeqIO.parse(myFastq1, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq1.close()
logger.info("Finished reading {}".format(args.filename1.name))

set1[fastqName1+"_R1_Len"] = pd.Series(forwardLen)
set1[fastqName1+"_R1_Qual"] = pd.Series(forwardAvg)
set1[fastqName1+"_R2_Len"] = pd.Series(reverseLen)
set1[fastqName1+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []

iter = 0

for record in SeqIO.parse(myFastq2, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq2.close()
logger.info("Finished reading {}".format(args.filename2.name))

set2[fastqName2+"_R1_Len"] = pd.Series(forwardLen)
set2[fastqName2+"_R1_Qual"] = pd.Series(forwardAvg)
set2[fastqName2+"_R2_Len"] = pd.Series(reverseLen)
set2[fastqName2+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []

iter = 0

for record in SeqIO.parse(myFastq3, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq3.close()
logger.info("Finished reading {}".format(args.filename3.name))

set3[fastqName3+"_R1_Len"] = pd.Series(forwardLen)
set3[fastqName3+"_R1_Qual"] = pd.Series(forwardAvg)
set3[fastqName3+"_R2_Len"] = pd.Series(reverseLen)
set3[fastqName3+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []

iter = 0

for record in SeqIO.parse(myFastq4, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq4.close()
logger.info("Finished reading {}".format(args.filename4.name))

set4[fastqName4+"_R1_Len"] = pd.Series(forwardLen)
set4[fastqName4+"_R1_Qual"] = pd.Series(forwardAvg)
set4[fastqName4+"_R2_Len"] = pd.Series(reverseLen)
set4[fastqName4+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []

iter = 0

for record in SeqIO.parse(myFastq5, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq5.close()
logger.info("Finished reading {}".format(args.filename5.name))

set5[fastqName5+"_R1_Len"] = pd.Series(forwardLen)
set5[fastqName5+"_R1_Qual"] = pd.Series(forwardAvg)
set5[fastqName5+"_R2_Len"] = pd.Series(reverseLen)
set5[fastqName5+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []

iter = 0

for record in SeqIO.parse(myFastq6, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq6.close()
logger.info("Finished reading {}".format(args.filename6.name))

set6[fastqName6+"_R1_Len"] = pd.Series(forwardLen)
set6[fastqName6+"_R1_Qual"] = pd.Series(forwardAvg)
set6[fastqName6+"_R2_Len"] = pd.Series(reverseLen)
set6[fastqName6+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []


iter = 0

for record in SeqIO.parse(myFastq7, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq7.close()
logger.info("Finished reading {}".format(args.filename7.name))

set7[fastqName7+"_R1_Len"] = pd.Series(forwardLen)
set7[fastqName7+"_R1_Qual"] = pd.Series(forwardAvg)
set7[fastqName7+"_R2_Len"] = pd.Series(reverseLen)
set7[fastqName7+"_R2_Qual"] = pd.Series(reverseAvg)

forwardLen = []
forwardAvg = []
reverseLen = []
reverseAvg = []

iter = 0

for record in SeqIO.parse(myFastq8, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])),end="")
                nameCheck = record.name
                forwardLen.append(len(record.seq))
                forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        elif(iter % 2 == 1):
                if((nameCheck != record.name) and (warned == 0)):   ## Confirm that next read has same name, i.e., is reverse read of previous
                        logger.warning("Read pairing for input file {} is broken!".format(myTitle[len(myTitle) - 1]))
                        pairingPreserved = 'N'
                        warned = 1
                reverseLen.append(len(record.seq))
                strand = record.description.split(" ")
                #print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
        iter = iter + 1

print(iter)
myFastq8.close()
logger.info("Finished reading {}".format(args.filename8.name))

set8[fastqName8+"_R1_Len"] = pd.Series(forwardLen)
set8[fastqName8+"_R1_Qual"] = pd.Series(forwardAvg)
set8[fastqName8+"_R2_Len"] = pd.Series(reverseLen)
set8[fastqName8+"_R2_Qual"] = pd.Series(reverseAvg)


print(set1.describe())
print(set2.describe())
print(set3.describe())
print(set4.describe())
print(set5.describe())
print(set6.describe())
print(set7.describe())
print(set8.describe())

lisAvgLenFwd = set1[fastqName1+'_R1_Len'].to_numpy(np.float64).tolist() + set2[fastqName2+'_R1_Len'].to_numpy(np.float64).tolist() + set3[fastqName3+'_R1_Len'].to_numpy(np.float64).tolist() + set4[fastqName4+"_R1_Len"].to_numpy(np.float64).tolist() + set5[fastqName5+"_R1_Len"].to_numpy(np.float64).tolist() + set6[fastqName6+"_R1_Len"].to_numpy(np.float64).tolist() + set7[fastqName7+"_R1_Len"].to_numpy(np.float64).tolist() + set8[fastqName8+"_R1_Len"].to_numpy(np.float64).tolist() 
medLenFwd = np.percentile(lisAvgLenFwd, 50)
print(medLenFwd)
#print(statistics.median(medLenFwd.values()))
avgLenFwd = np.mean(lisAvgLenFwd)
print(avgLenFwd)
#lisMinLenFwd = set1[fastqName1+'_R1_Len'].to_numpy(np.float64) 
minLenFwd = min(lisAvgLenFwd)
print(minLenFwd)
lisFwdAvgPHRED = set1[fastqName1+'_R1_Qual'].to_numpy(np.float64).tolist() + set2[fastqName2+'_R1_Qual'].to_numpy(np.float64).tolist() + set3[fastqName3+'_R1_Qual'].to_numpy(np.float64).tolist() + set4[fastqName4+"_R1_Qual"].to_numpy(np.float64).tolist() + set5[fastqName5+"_R1_Qual"].to_numpy(np.float64).tolist() + set6[fastqName6+"_R1_Qual"].to_numpy(np.float64).tolist() + set7[fastqName7+"_R1_Qual"].to_numpy(np.float64).tolist() + set8[fastqName8+"_R1_Qual"].to_numpy(np.float64).tolist()
medFwdAvgPHRED = np.percentile(lisFwdAvgPHRED, 50)
print(medFwdAvgPHRED)
#lqFwdAvgPHRED = np.percentile(set1[fastqName1+'_R1_Qual', fastqName2+'_R1_Qual'].values.tolist(), 25)
#uqFwdAvgPHRED = np.percentile(set1[fastqName1+'_R1_Qual', fastqName2+'_R1_Qual'].values.tolist(), 75)
avgFwdAvgPHRED = np.mean(lisFwdAvgPHRED)
print(avgFwdAvgPHRED)
minFwdAvgPHRED = min(lisFwdAvgPHRED)
print(minFwdAvgPHRED)

lisMedLenRev = set1[fastqName1+'_R2_Len'].to_numpy(np.float64).tolist() + set2[fastqName2+'_R2_Len'].to_numpy(np.float64).tolist() + set3[fastqName3+'_R2_Len'].to_numpy(np.float64).tolist() + set4[fastqName4+"_R2_Len"].to_numpy(np.float64).tolist() + set5[fastqName5+"_R2_Len"].to_numpy(np.float64).tolist() + set6[fastqName6+"_R2_Len"].to_numpy(np.float64).tolist() + set7[fastqName7+"_R2_Len"].to_numpy(np.float64).tolist() + set8[fastqName8+"_R2_Len"].to_numpy(np.float64).tolist()
medLenRev = np.percentile(lisMedLenRev, 50)
print(medLenRev)
avgLenRev = np.mean(lisMedLenRev)
print(avgLenRev)
minLenRev = min(lisMedLenRev)
print(minLenRev)
lisRevAvgPHRED = set1[fastqName1+'_R2_Qual'].to_numpy(np.float64).tolist() + set2[fastqName2+'_R2_Qual'].to_numpy(np.float64).tolist() + set3[fastqName3+'_R2_Qual'].to_numpy(np.float64).tolist() + set4[fastqName4+"_R2_Qual"].to_numpy(np.float64).tolist() + set5[fastqName5+"_R2_Qual"].to_numpy(np.float64).tolist() + set6[fastqName6+"_R2_Qual"].to_numpy(np.float64).tolist() + set7[fastqName7+"_R2_Qual"].to_numpy(np.float64).tolist() + set8[fastqName8+"_R2_Qual"].to_numpy(np.float64).tolist()
medRevAvgPHRED = np.percentile(lisRevAvgPHRED, 50)
print(medRevAvgPHRED)
#lqRevAvgPHRED = np.percentile(set1[fastqName1+'_R2_Qual', fastqName2+'_R2_Qual',].values.tolist(), 25)
#uqRevAvgPHRED = np.percentile(set1[fastqName1+'_R2_Qual', fastqName2+'_R2_Qual',].values.tolist(), 75)
avgRevAvgPHRED = np.mean(lisRevAvgPHRED)
print(avgRevAvgPHRED)
minRevAvgPHRED = min(lisRevAvgPHRED)
print(minRevAvgPHRED)

#forwardAvg = np.random.normal(size = 500000) + 30
#reverseAvg = 1.5*np.random.normal(size = 500000) + 25

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))

if(pairingPreserved == 'Y'):
        fig.text(0.5, 0.04, 'Average Read Quality', ha='center')
else:
        fig.text(0.5, 0.04, 'Read Quality - Pairs Broken!', ha='center')

fig.text(0.04, 0.5, 'Read Counts', va='center', rotation='vertical')

## Plot PHRED Quality for R1 reads as 1D histogram
axes[0].hist(lisFwdAvgPHRED, bins = 25, color='blue')
axes[0].set_title("R1 Quality " + newTitle)

## Plot PHRED Quality for R2 reads as 1D histogram
axes[1].hist(lisRevAvgPHRED, bins = 25, color='red')
axes[1].set_title("R2 Quality " + newTitle)

## create output folder
if(args.outputType != 'J'):
        os.mkdir('QC_' + newTitle)

## write figure file as .png
if(args.outputType == 'F'):
        fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/fwd_and_revPHRED_histo.png')

## output average per-read PHRED scores into a .csv file
        totalPHRED = { "R1_Reads_Avg_PHRED" : lisFwdAvgPHRED, "R2_Reads_Avg_PHRED" : lisRevAvgPHRED }
        dfTotalPHRED = pd.DataFrame(totalPHRED)
        dfTotalPHRED.to_csv('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/fwd_and_revPHRED_allEight.csv')

## write summary statistics into human readable text file
if(args.outputType != 'N'):
        summary = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/ReadStatistics.README.txt', 'w')
        summary.write("Summary Statistics for " + newTitle + "\n")
        summary.write("%s\t%0.2f" % ("R1 median length:", medLenFwd) + ",\t%s\t%0.2f\n" % ("R2 median length:", medLenRev))
        summary.write("%s\t%0.2f" % ("R1 mean length:", avgLenFwd) + ",\t%s\t%0.2f\n" % ("R1 mean length:", avgLenRev))
        summary.write("%s\t%s" % ("R1 minimum length:", str(minLenFwd)) + ",\t%s\t%s\n" % ("R2 minimum length:", str(minLenRev)))
        summary.write("%s\t%0.2f" % ("R1 median PHRED:", medFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 median PHRED:", medRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 mean PHRED:", avgFwdAvgPHRED)   + ",\t%s\t%0.2f\n" % ("R2 mean PHRED:", avgRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 minimum PHRED:", minFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 minimum PHRED:", minRevAvgPHRED))
        summary.close()

        summary1 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName1 + 'Statistics.README.txt', 'a')
        summary1.write("Summary Statistics for " + fastqName1 + "\n")
        set1.describe().to_string(summary1)
        summary1.close()

        summary2 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName2 + 'Statistics.README.txt', 'a')
        summary2.write("Summary Statistics for " + fastqName2 + "\n")
        set2.describe().to_string(summary2)
        summary2.close()

        summary3 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName3 + 'Statistics.README.txt', 'a')
        summary3.write("Summary Statistics for " + fastqName3 + "\n")
        set3.describe().to_string(summary3)
        summary3.close()

        summary4 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName4 + 'Statistics.README.txt', 'a')
        summary4.write("Summary Statistics for " + fastqName4 + "\n")
        set4.describe().to_string(summary4)
        summary4.close()

        summary5 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName5 + 'Statistics.README.txt', 'a')
        summary5.write("Summary Statistics for " + fastqName5 + "\n")
        set5.describe().to_string(summary5)
        summary5.close()

        summary6 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName6 + 'Statistics.README.txt', 'a')
        summary6.write("Summary Statistics for " + fastqName6 + "\n")
        set6.describe().to_string(summary6)
        summary6.close()

        summary7 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName7 + 'Statistics.README.txt', 'a')
        summary7.write("Summary Statistics for " + fastqName7 + "\n")
        set7.describe().to_string(summary7)
        summary7.close()

        summary8 = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/' + fastqName8 + 'Statistics.README.txt', 'a')
        summary8.write("Summary Statistics for " + fastqName8 + "\n")
        set8.describe().to_string(summary8)
        summary8.close()

if(args.outputType == 'N'):
        summary = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_' + newTitle + '/AllSevenReadStatistics.README.txt', 'w')
        summary.write("Summary Statistics for " + newTitle + "\n")
        summary.write("%s\t%0.2f" % ("R1 median length:", medLenFwd) + ",\t%s\t%0.2f\n" % ("R2 median length:", medLenRev))
        summary.write("%s\t%0.2f" % ("R1 mean length:", avgLenFwd) + ",\t%s\t%0.2f\n" % ("R1 mean length:", avgLenRev))
        summary.write("%s\t%s" % ("R1 minimum length:", str(minLenFwd)) + ",\t%s\t%s\n" % ("R2 minimum length:", str(minLenRev)))
        summary.write("%s\t%0.2f" % ("R1 median PHRED:", medFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 median PHRED:", medRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 mean PHRED:", avgFwdAvgPHRED)   + ",\t%s\t%0.2f\n" % ("R2 mean PHRED:", avgRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 minimum PHRED:", minFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 minimum PHRED:", minRevAvgPHRED))
        summary.close()



## write summary statistics into json file
#jsonDict = {
#        newTitle : {"read length" : { "median" : { "R1" : medLenFwd, "R2" : medLenRev }, "average" : { "R1" : avgLenFwd, "R2" : avgLenRev }, "minimum" : { "R1" : minLenFwd, "R2" : minLenRev } }, 
#                    "read PHRED quality" : { "median" : { "R1" : medFwdAvgPHRED, "R2" : medRevAvgPHRED }, "average" : {"R1" : avgFwdAvgPHRED, "R2" : avgRevAvgPHRED}, "minimum" : {"R1" : minFwdAvgPHRED, "R2" : minRevAvgPHRED} }, "Read Pairing Preserved" : pairingPreserved } }
#if(args.outputType == 'J'):
#        with open(basePath + "/" + inPath1[0] + "/QC_" + newTitle + ".json", "w") as jsummary:
#                json.dump(jsonDict, jsummary)
#else:
#        with open("/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/QC_" + newTitle + "/ReadStatistics.json", "w") as jsummary:
#                json.dump(jsonDict, jsummary)

