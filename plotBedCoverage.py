import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt

## Script for calculating sequence length and average PHRED score per read
## Outputs comma-delimited table to stdout
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

parser = argparse.ArgumentParser(description='Plots coverage from a .bed file as a line plot', usage="plotBedCoverage.py filepath/filename.bedGraph")

parser.add_argument('filename', nargs='+', type=ext_check('.bedGraph', argparse.FileType('r')))

parser.add_argument('--window', '-w', default='10', type=int)

args = parser.parse_args()

# echo input file name
print(args.filename[0].name)


myTitle = re.split(r'[\/]', args.filename[0].name)

print(myTitle[len(myTitle) - 1])

shortTitle = re.split(r'[\.]', myTitle[len(myTitle) - 1])

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

## Coerce original args.window values to divisible by 5
if window < 6:
        window = 5
elif window > 14:
        window = 20

print("window = " + str(window))

myCoverage = open(args.filename[0].name, "r")

for line in myCoverage:
        lineData = re.split(r'\s+', line)
        if(iter % window == 0):
                coordinates.append(lineData[1])
                # add plotting point to smoothing window
                smooth.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage.append(statistics.mean(smooth))
                smooth = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                referenceName = lineData[0]
                smooth.append(int(lineData[2]))
        else:
                smooth.append(int(lineData[2]))
        iter = iter + 1

myXticks = []
iter = 0
# Generate X-axis labels
for count in coordinates:
        if(iter % 100 == 0):
                myXticks.append(count)
        iter = iter + 1


## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(10,7), gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Plot Genome Coverage as line plot
axes.plot(coordinates, coverage)
axes.set_xticks(myXticks, myXticks, rotation='vertical')
axes.set_title("Sample " + shortTitle[0] + " Reads")
axes.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates')
axes.set_ylabel('Coverage (X) at Position')
#axes.margins(0.2)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/w' + str(window) + '_' + shortTitle[0] + '_to_' + referenceName + '.png')

