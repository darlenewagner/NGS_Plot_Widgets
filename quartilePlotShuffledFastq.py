import os, sys, re, csv, statistics, time
import argparse, logging, warnings, json
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
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


logger = logging.getLogger("quartilePlotShuffledFastq.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq (Expects single fastq input)', usage="python fullPlotShuffledFastq.py filepath/filename.fastq(.gz) --outputType [F/J/N]")

parser.add_argument('filename', type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='F', choices=['F', 'J', 'N'], help="--outputType F for full output,... J for .json only, and N for no image.")

parser.add_argument('--force', '-f', default='N', choices=['Y', 'N'], help="--force Y to remove previous folder with same name and remake it.")

args = parser.parse_args()

print(args.filename.name)

basePath = os.getcwd()
inPath = os.path.split(args.filename.name)

myTitle = re.split(r'\/', args.filename.name)
newTitle = re.sub('\.f(ast)?q(\.gz)', '', myTitle[len(myTitle) - 1])
#print(newTitle)

csvRow = []
forwardLen = []
reverseLen = []
forwardAvg = []
reverseAvg = []

iter = 0
nameCheck = ''
pairingPreserved = 'Y'
warned = 0

logger.info("Reading file {}/{}".format(basePath, inPath[0]))

myFastq = ''
tempFas = re.sub(r'.gz', '', args.filename.name)

if(args.filename.name.endswith('.gz') == True ):
        os.system("gunzip -c {} > {}".format(args.filename.name, tempFas))
        myFastq = open(tempFas, "r")
else:
        myFastq = open(args.filename.name, "r")

for record in SeqIO.parse(myFastq, "fastq"):
        if(iter % 2 == 0):
                #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])), end="")
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

logger.info("Finished reading {}".format(args.filename.name))


## Compute medians, averages, and quartiles
medLenFwd = statistics.median(forwardLen)
avgLenFwd = statistics.mean(forwardLen)
minLenFwd = min(forwardLen)
medFwdAvgPHRED = statistics.median(forwardAvg)
lqFwdAvgPHRED = np.percentile(forwardAvg, 25)
uqFwdAvgPHRED = np.percentile(forwardAvg, 75)
avgFwdAvgPHRED = statistics.mean(forwardAvg)
minFwdAvgPHRED = min(forwardAvg)
medLenRev = statistics.median(reverseLen)
avgLenRev = statistics.mean(reverseLen)
minLenRev = min(reverseLen)
medRevAvgPHRED = statistics.median(reverseAvg)
lqRevAvgPHRED = np.percentile(reverseAvg, 25)
uqRevAvgPHRED = np.percentile(reverseAvg, 75)
avgRevAvgPHRED = statistics.mean(reverseAvg)
minRevAvgPHRED = min(reverseAvg)

#forwardAvg = np.random.normal(size = 500000) + 30
#reverseAvg = 1.5*np.random.normal(size = 500000) + 25

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))

if(pairingPreserved == 'Y'):
        fig.text(0.5, 0.04, 'Average Read Quality', ha='center')
else:
        fig.text(0.5, 0.04, 'Read Quality - Pairs Broken!', ha='center')

fig.text(0.04, 0.5, 'Read Counts', va='center', rotation='vertical')

## Plot PHRED Quality for R1 reads as 1D histogram
axes[0].hist(forwardAvg, bins = 25, color='blue')
axes[0].set_title("R1 Quality " + newTitle)

## Plot PHRED Quality for R2 reads as 1D histogram
axes[1].hist(reverseAvg, bins = 25, color='red')
axes[1].set_title("R2 Quality " + newTitle)

## create output folder
if((args.outputType != 'J') and (args.force == 'N') ):
        try:
                os.mkdir(basePath + "/" + inPath[0] + "/QC_" + newTitle)
        except OSError:
                print("Output folder {} already exists!".format(basePath + "/" + inPath[0] + "/QC_" + newTitle))
                sys.exit("Command os.mkdir exception - terminating script.")
elif(args.outputType != 'J'):
        os.system("rm -vr {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle))
        os.mkdir(basePath + "/" + inPath[0] + "/QC_" + newTitle)
        logger.info("Created new {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle))

## write figure file as .png
if(args.outputType == 'F'):
        fig.savefig(basePath + "/" + inPath[0] + "/QC_" + newTitle + "/fwd_and_revPHRED_histo.png")
        logger.info("Finished plotting {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle + "/fwd_and_revPHRED_histo.png"))

## write summary statistics into human readable text file
if(args.outputType != 'J'):
        summary = open(basePath + "/" + inPath[0] + "/QC_" + newTitle + "/ReadStatistics.README.txt", "w")
        summary.write("Summary Statistics for " + newTitle + "\n")
        summary.write("%s\t%0.2f" % ("R1 median length:", medLenFwd) + ",\t%s\t%0.2f\n" % ("R2 median length:", medLenRev))
        summary.write("%s\t%0.2f" % ("R1 mean length:", avgLenFwd) + ",\t%s\t%0.2f\n" % ("R1 mean length:", avgLenRev))
        summary.write("%s\t%s" % ("R1 minimum length:", str(minLenFwd)) + ",\t%s\t%s\n" % ("R2 minimum length:", str(minLenRev)))
        summary.write("%s\t%0.2f" % ("R1 PHRED Upper Quartile:", uqFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 PHRED Upper Quartile:", uqRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 median PHRED:", medFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 median PHRED:", medRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 PHRED Lower Quartile:", lqFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 PHRED Lower Quartile:", lqRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 average PHRED:", avgFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 average PHRED:", avgRevAvgPHRED))
        summary.write("%s\t%0.2f" % ("R1 minimum PHRED:", minFwdAvgPHRED) + ",\t%s\t%0.2f\n" % ("R2 minimum PHRED:", minRevAvgPHRED))
        summary.close()
        time.sleep(0.25)
        logger.info("Finished writing {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle + "/ReadStatistics.README.txt"))

## write summary statistics into json file

jsonDict = {
        newTitle : {"read length" : { "median" : { "R1" : medLenFwd, "R2" : medLenRev }, "average" : { "R1" : avgLenFwd, "R2" : avgLenRev }, "minimum" : { "R1" : minLenFwd, "R2" : minLenRev } }, 
                    "read PHRED quality" : {"upperQuartile" : {"R1" : uqFwdAvgPHRED, "R2" : uqRevAvgPHRED }, "median" : { "R1" : medFwdAvgPHRED, "R2" : medRevAvgPHRED }, "lowerQuartile" : {"R1" : lqFwdAvgPHRED, "R2" : lqRevAvgPHRED }, "average" : {"R1" : avgFwdAvgPHRED, "R2" : avgRevAvgPHRED}, "minimum" : {"R1" : minFwdAvgPHRED, "R2" : minRevAvgPHRED} }, "Read Pairing Preserved" : pairingPreserved } }

if((args.outputType == 'J') and (args.force == 'N')):
        try:
                with open(basePath + "/" + inPath[0] + "/QC_" + newTitle + "_quartile.json", "w") as jsummary:
                        json.dump(jsonDict, jsummary)
        except OSError:
                print("Output file {} already exists!".format(basePath + "/" + inPath[0] + "/QC_" + newTitle + "_quartile.json"))
                sys.exit("File open exception - terminating script.")
elif(args.outputType == 'J'):
        os.system("rm -v {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle + "_quartile.json"))
        with open(basePath + "/" + inPath[0] + "/QC_" + newTitle + "_quartile.json", "w") as jsummary:
                json.dump(jsonDict, jsummary)
        logger.info("Created new {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle + "_quartile.json"))
else:
        with open(basePath + "/" + inPath[0] + "/QC_" + newTitle + "/ReadStatistics.json", "w") as jsummary:
                json.dump(jsonDict, jsummary)

time.sleep(0.25)
logger.info("Finished writing {}".format(basePath + "/" + inPath[0] + "/QC_" + newTitle + "/ReadStatistics.json"))
