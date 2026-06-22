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


parser = argparse.ArgumentParser(description='Creates a line plot of coverage from a .bed file with its corresponding .faidx file defining the x-axis', usage="plotBedCoverage.py filepath/filename.bed filepath/filename.fasta.faidx")

parser.add_argument('filename1', type=ext_check1('.bedGraph', '.bed', argparse.FileType('r')))

parser.add_argument('filename2', type=ext_check2('.fasta.fai', '.fa.fai', argparse.FileType('r')))

## The parameter, --rotavirus, determines whether or not to annotate segments in the line plot.
parser.add_argument('--rotavirus', '-r', default='N', help="Add --rotavirus Y to command line for annotation of rotavirus segments, otherwise, omit --rotavirus parameter.")

parser.add_argument('--virusType', '-v', default='miscellaneous', help="Add --virusType virus_name to create annotated output file name.")

parser.add_argument('--ignoreLow', '-i', default='60', type=int, help="For not plotting coverage between 20X and an upper limit given by --ignoreLow")

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

rotavirus = args.rotavirus

first_line = "";
myCoverage = open(args.filename1.name, "r")

with open(args.filename2.name, "r") as myLength:
        first_line = myLength.readline().rstrip('\n')

genomeLength = re.split(r'\s+', first_line)

myInitialCoord = 1
myEndCoord = int(genomeLength[1])
print(myEndCoord)

for line in myCoverage:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        coordinates.append(int(lineData[1]))
        coverage.append(int(lineData[2]))
        iter = iter + 1

idx = 0
adjustedCoordinate = []
hiAdjustedCoordinates = []
hiAdjustedCoverage = []
loAdjustedCoordinates = []
loAdjustedCoverage = []


while myInitialCoord < myEndCoord - 1:
        if(( myInitialCoord == coordinates[idx] ) and (idx < len(coordinates) - 1) ):
                if(coverage[idx] < 20):
                        loAdjustedCoverage.append(coverage[idx])
                        loAdjustedCoordinates.append( myInitialCoord )
                        adjustedCoordinate.append( myInitialCoord )
                else:
                        hiAdjustedCoverage.append(coverage[idx])
                        hiAdjustedCoordinates.append( myInitialCoord )
                        adjustedCoordinate.append( myInitialCoord )
                idx = idx + 1
                myInitialCoord = myInitialCoord + 1
        else:
                #hiAdjustedCoverage.append( 0 )
                #hiAdjustedCoordinates.append( myInitialCoord )
                loAdjustedCoverage.append( 0 )
                loAdjustedCoordinates.append( myInitialCoord )
                adjustedCoordinate.append( myInitialCoord )
                myInitialCoord = myInitialCoord + 1

for t in loAdjustedCoverage:
     print(t)   
                
myXticks = []
iter = 0
# Generate X-axis labels
for count in adjustedCoordinate:
        if(iter % 500 == 0):
                myXticks.append(count)
        iter = iter + 1


## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(10,7), gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Mask out hiAdjustedCoverage between 20X and 60X
npAdjustedCoverage = np.array(hiAdjustedCoverage)
y_filtered = np.where(npAdjustedCoverage < int(args.ignoreLow), np.nan, npAdjustedCoverage)

## Plot Genome Coverage as line plot
axes.plot(hiAdjustedCoordinates, y_filtered, linestyle='-', marker='.')
axes.plot(loAdjustedCoordinates, loAdjustedCoverage, linestyle='None', color='red', marker='.')
if((rotavirus == 'Y') or (rotavirus == 'y')):
        axes.text(0.078, 0.8, 'Segment 1', transform=axes.transAxes)
        axes.text(0.22, 0.75, 'Segment 2', transform=axes.transAxes)
        axes.text(0.34, 0.7, 'Segment 3', transform=axes.transAxes)
        axes.text(0.46, 0.65, 'Segment 4', transform=axes.transAxes)
        axes.text(0.58, 0.6, 'Seg. 5', transform=axes.transAxes)
        axes.text(0.65, 0.55, 'Seg. 6', transform=axes.transAxes)
        axes.text(0.72, 0.5, 'Seg. 7', transform=axes.transAxes)
        axes.text(0.775, 0.45, 'Seg. 8', transform=axes.transAxes)
        axes.text(0.825, 0.4, 'Seg. 9', transform=axes.transAxes)
        axes.text(0.875, 0.35, 'Seg. 10', transform=axes.transAxes)
        axes.text(0.915, 0.3, 'Seg. 11', transform=axes.transAxes)
        args.virusType = 'rotavirus'
        
axes.set_xticks(myXticks, myXticks, rotation='vertical')
axes.set_title( shortTitle + " Read Mapping Depth")
axes.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates')
axes.set_ylabel('Coverage (X) at Position')
#axes.margins(0.2)

fig.savefig('/scicomp/home-pure/ydn3/NGS_Plot_Widgets/Hidden_Files/' + args.virusType + '_' + shortTitle + '_to_' + referenceName + '.png')

