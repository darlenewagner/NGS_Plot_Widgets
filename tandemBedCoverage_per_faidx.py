import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt

## Script for calculating sequence length and average PHRED score per read
## Requires a .bed or .bedGraph file as input followed by a corresponding faidx file
## Outputs log messages to stdout and .png graphic to a file via Matplotlib
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check1(expected_ext, another_ext, openner):
        def extension(filename):
                if not (filename.endswith(expected_ext) or filename.endswith(another_ext)):
                        raise ValueError()
                return openner(filename)
        return extension

def ext_check2(fai_ext, other_fai_ext, openner):
        def extension(filename):
                if not (filename.endswith(fai_ext) or filename.endswith(other_fai_ext)):
                        raise ValueError()
                return openner(filename)
        return extension


parser = argparse.ArgumentParser(description='Plots coverage from a .bed file and its .faidx file as a line plot', usage="plotBedCoverage.py filepath/filename1.bed filepath/filename2.bed filepath/filename.fasta.faidx")

parser.add_argument('filename1', type=ext_check1('.bedGraph', '.bed', argparse.FileType('r')))

parser.add_argument('filename2', type=ext_check1('.bedGraph', '.bed', argparse.FileType('r')))

parser.add_argument('filename3', type=ext_check2('.fasta.fai', '.fa.fai', argparse.FileType('r')))

parser.add_argument('--window', '-w', default='10', type=int)

args = parser.parse_args()

# echo input file name
print(args.filename1.name)


myTitle = re.split(r'[\/.]', args.filename1.name)

print(myTitle[len(myTitle) - 2])

shortTitle = myTitle[1]

csvRow = []
forwardLen = []
reverseLen = []
coordinates = []
coverage = []
referenceName = ''

iter = 0
# avoid plotting every point, plot average of every 5th, 10th, or 20th point
smooth = []

window = int(args.window)

## Coerce original args.window values to be divisible by 5
if window < 6:
        window = 5
elif window > 14:
        window = 20

print("window = " + str(window))

first_line = "";
myCoverage1 = open(args.filename1.name, "r")
myCoverage2 = open(args.filename2.name, "r")

with open(args.filename3.name, "r") as myLength:
        first_line = myLength.readline().rstrip('\n')

genomeLength = re.split(r'\s+', first_line)

myInitialCoord = 1
myEndCoord = int(genomeLength[1])
print(myEndCoord)

for line in myCoverage1:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        coordinates.append(int(lineData[1]))
        coverage.append(int(lineData[2]))
        iter = iter + 1

idx = 0
adjustedCoordinate1 = []
hiAdjustedCoordinates1 = []
hiAdjustedCoverage1 = []
loAdjustedCoordinates1 = []
loAdjustedCoverage1 = []


while myInitialCoord < myEndCoord - 1:
        if(( myInitialCoord == coordinates[idx] ) and (idx < len(coordinates) - 1) ):
                if(coverage[idx] < 21):
                        loAdjustedCoverage1.append(coverage[idx])
                        loAdjustedCoordinates1.append( myInitialCoord )
                        adjustedCoordinate1.append( myInitialCoord )
                else:
                        hiAdjustedCoverage1.append(coverage[idx])
                        hiAdjustedCoordinates1.append( myInitialCoord )
                        adjustedCoordinate1.append( myInitialCoord )
                idx = idx + 1
                myInitialCoord = myInitialCoord + 1
        else:
                #hiAdjustedCoverage.append( 0 )
                #hiAdjustedCoordinates.append( myInitialCoord )
                loAdjustedCoverage1.append( 0 )
                loAdjustedCoordinates1.append( myInitialCoord )
                adjustedCoordinate1.append( myInitialCoord )
                myInitialCoord = myInitialCoord + 1

iter = 0

coordinates2 = []
coverage2 = []

for line in myCoverage2:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        coordinates2.append(int(lineData[1]))
        coverage2.append(int(lineData[2]))
        iter = iter + 1

idx = 0
adjustedCoordinate2 = []
hiAdjustedCoordinates2 = []
hiAdjustedCoverage2 = []
loAdjustedCoordinates2 = []
loAdjustedCoverage2 = []


myInitialCoord = 1

while myInitialCoord < myEndCoord - 1:
        if(( myInitialCoord == coordinates2[idx] ) and (idx < len(coordinates2) - 1) ):
                if(coverage2[idx] < 21):
                        loAdjustedCoverage2.append(coverage2[idx])
                        loAdjustedCoordinates2.append( myInitialCoord )
                        adjustedCoordinate2.append( myInitialCoord )
                else:
                        hiAdjustedCoverage2.append(coverage2[idx])
                        hiAdjustedCoordinates2.append( myInitialCoord )
                        adjustedCoordinate2.append( myInitialCoord )
                idx = idx + 1
                myInitialCoord = myInitialCoord + 1
        else:
                #hiAdjustedCoverage.append( 0 )
                #hiAdjustedCoordinates.append( myInitialCoord )
                loAdjustedCoverage2.append( 0 )
                loAdjustedCoordinates2.append( myInitialCoord )
                adjustedCoordinate2.append( myInitialCoord )
                myInitialCoord = myInitialCoord + 1




  #for t in loAdjustedCoverage1:
  #     print(t)

myXticks = []
iter = 0
# Generate X-axis labels
for count in adjustedCoordinate1:
        if(iter % 500 == 0):
                myXticks.append(count)
        iter = iter + 1


## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(10,7), gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Plot 1 Genome Coverage as line plot
axes[0].plot(hiAdjustedCoordinates1, hiAdjustedCoverage1, linestyle='-', marker='.')
axes[0].plot(loAdjustedCoordinates1, loAdjustedCoverage1, linestyle='None', color='red', marker='.')
axes[0].set_xticks(myXticks, myXticks, rotation='vertical')
axes[0].set_title("Sample " + shortTitle + " Reads")
axes[0].set_xlabel('Reference Genome, ' + referenceName + ', Coordinates')
axes[0].set_ylabel('Coverage (X) at Position')
#axes.margins(0.2)


myXticks2 = []
iter = 0
# Generate X-axis labels
for count in adjustedCoordinate2:
        if(iter % 500 == 0):
                myXticks2.append(count)
        iter = iter + 1


## Plot 2 Genome Coverage as line plot
axes[1].plot(hiAdjustedCoordinates2, hiAdjustedCoverage2, linestyle='-', marker='.')
axes[1].plot(loAdjustedCoordinates2, loAdjustedCoverage2, linestyle='None', color='red', marker='.')
axes[1].set_xticks(myXticks2, myXticks2, rotation='vertical')
axes[1].set_title("Sample " + shortTitle + " Reads")
axes[1].set_xlabel('Reference Genome, ' + referenceName + ', Coordinates')
axes[1].set_ylabel('Coverage (X) at Position')


fig.savefig('/scicomp/home-pure/ydn3/NGS_Plot_Widgets/Hidden_Files/tandem_win' + str(window) + '_' + shortTitle + '_to_' + referenceName + '.png')

