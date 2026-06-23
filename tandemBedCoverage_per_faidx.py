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

## tandemBedCoverage_per_faidx.py not recommended for rotavirus or other segmented viral genomes.
## parser.add_argument('--rotavirus', '-r', default='N', help="Add --rotavirus Y to command line for annotation of rotavirus segments, otherwise, omit --rotavirus parameter.")

parser.add_argument('--virusType', '-v', default='miscellaneous', help="Add --virusType virus_name to create annotated output file name.")

parser.add_argument('--plotTitle1', '-p1', default='', help="Provide --plotTitle1 string if input filename 1 is not the desired title.")

parser.add_argument('--ignoreLow1', '-i1', default='60', type=int, help="For not plotting coverage between 20X and an upper limit given by --ignoreLow")

parser.add_argument('--plotTitle2', '-p2', default='', help="Provide --plotTitle2 string if input filename 2 is not the desired title.")

parser.add_argument('--ignoreLow2', '-i2', default='60', type=int, help="For not plotting coverage between 20X and an upper limit given by --ignoreLow")

args = parser.parse_args()

# echo input file name
print(args.filename1.name)


myTitle1 = re.split(r'[\/.]', args.filename1.name)
print(myTitle1[len(myTitle1) - 2])
shortTitle1 = myTitle1[1]

myTitle2 = re.split(r'[\/.]', args.filename2.name)
print(myTitle2[len(myTitle2) - 2])
shortTitle2 = myTitle2[1]

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

maxCoverage = 0;

while myInitialCoord < myEndCoord - 1:
        if(( myInitialCoord == coordinates[idx] ) and (idx < len(coordinates) - 1) ):
                if(coverage[idx] < 21):
                        loAdjustedCoverage1.append(coverage[idx])
                        loAdjustedCoordinates1.append( myInitialCoord )
                        adjustedCoordinate1.append( myInitialCoord )
                else:
                        hiAdjustedCoverage1.append(coverage[idx])
                        if(coverage[idx] > maxCoverage):
                                maxCoverage = coverage[idx]
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
                        if(coverage2[idx] > maxCoverage):
                                maxCoverage = coverage2[idx]
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

yLimit = maxCoverage + (maxCoverage % 1000)

myYtick = []
ii = 0
if(yLimit > 2500):
        while( ii < yLimit ):
                myYtick.append(ii)
                ii = ii + 1000
else:
        while( ii < yLimit ):
                myYtick.append(ii)
                ii = ii + 500

## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig = plt.figure( figsize=(10,7) )
##gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Mask out hiAdjustedCoverage1 between 20X and 60X
npAdjustedCoverage1 = np.array(hiAdjustedCoverage1)
y_filtered_1 = np.where(npAdjustedCoverage1 < int(args.ignoreLow1), np.nan, npAdjustedCoverage1)

axe1 = fig.add_axes([0.08, 0.6, 0.89, 0.27])

## Plot 1 Genome Coverage as line plot
axe1.plot(hiAdjustedCoordinates1, y_filtered_1, linestyle='-', marker='.')
axe1.plot(loAdjustedCoordinates1, loAdjustedCoverage1, linestyle='None', color='red', marker='.')
axe1.set_xticks(myXticks, myXticks, rotation='vertical')
axe1.set_yticks(myYtick, myYtick)
axe1.set_title( shortTitle1 + " Read Mapping Depth")
axe1.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates')
axe1.set_ylabel('Coverage (X) at Position')
#axes.margins(0.2)


myXticks2 = []
iter = 0
# Generate X-axis labels
for count in adjustedCoordinate2:
        if(iter % 500 == 0):
                myXticks2.append(count)
        iter = iter + 1

## Mask out hiAdjustedCoverage1 between 20X and 60X
npAdjustedCoverage2 = np.array(hiAdjustedCoverage2)
y_filtered_2 = np.where(npAdjustedCoverage2 < int(args.ignoreLow2), np.nan, npAdjustedCoverage2)

axe2 = fig.add_axes([0.08, 0.15, 0.89, 0.27])

## Plot 2 Genome Coverage as line plot
axe2.plot(hiAdjustedCoordinates2, y_filtered_2, linestyle='-', marker='.')
axe2.plot(loAdjustedCoordinates2, loAdjustedCoverage2, linestyle='None', color='red', marker='.')
axe2.set_xticks(myXticks2, myXticks2, rotation='vertical')
axe2.set_yticks(myYtick, myYtick)
axe2.set_title(shortTitle2 + " Read Mapping Depth")
axe2.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates')
axe2.set_ylabel('Coverage (X) at Position')


fig.savefig('/scicomp/home-pure/ydn3/NGS_Plot_Widgets/Hidden_Files/' + args.virusType + '_' + shortTitle1 + '_' + shortTitle2 + '_to_' + referenceName + '.png')

