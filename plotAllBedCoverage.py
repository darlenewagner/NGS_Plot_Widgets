import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

## Script for plotting genome coverages and SNPs position coverages
## Outputs .png file of line plot 
## Requires Biopython, Numpy, and Matplotlib

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

parser = argparse.ArgumentParser(description='Plots coverage from two .bed files line plots', usage="plotBedCoverage.py filepath/filename1.bed filepath/filename2.bed filepath/filename3.bed filepath/filename4.bed filepath/filename5.bed filepath/filename6.bed")

parser.add_argument('filename1', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename2', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('snpFile1', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile2', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('filename3', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename4', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('snpFile3', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile4', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('filename5', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename6', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('snpFile5', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile6', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('--window', '-w', default='1', type=int)

args = parser.parse_args()

# echo input file name
print(args.filename1.name)
print(args.filename3.name)
print(args.filename5.name)

myTitle1 = re.split(r'[\/]', args.filename1.name)
myTitle2 = re.split(r'[\/]', args.filename3.name)
myTitle3 = re.split(r'[\/]', args.filename5.name)

#print(myTitle[len(myTitle) - 1])

shortTitle1 = re.split(r'[\.]', myTitle1[len(myTitle1) - 1])
shortTitle1t = shortTitle1[0].split(sep='_S')

shortTitle2 = re.split(r'[\.]', myTitle2[len(myTitle2) - 1])
shortTitle2t = shortTitle2[0].split(sep='_S')

shortTitle3 = re.split(r'[\.]', myTitle3[len(myTitle3) - 1])
shortTitle3t = shortTitle3[0].split(sep='_S')


csvRow = []
forwardLen = []
reverseLen = []
coordinates1 = []
coordinates2 = []

coverage1 = []
coverage2 = []

coordinates3 = []
coordinates4 = []

coverage3 = []
coverage4 = []

coordinates5 = []
coordinates6 = []

coverage5 = []
coverage6 = []

referenceName = ''

iter = 0
# avoid plotting every point, plot average of every 5th, 10th, or 20th point
smooth1 = []
smooth2 = []
smooth3 = []
smooth4 = []
smooth5 = []
smooth6 = []

window = int(args.window)

## Coerce original args.window values to divisible by 5
if window < 1:
        window = 1
elif window > 14:
        window = 20

print("window = " + str(window))

myCoverage1 = open(args.filename1.name, "r")

for line in myCoverage1:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates1.append(lineData[1])
                # add plotting point to smoothing window
                smooth1.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage1.append(statistics.mean(smooth1))
                smooth1 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth1.append(int(lineData[2]))
        else:
                smooth1.append(int(lineData[2]))
        iter = iter + 1

myCoverage1.close()

myCoverage2 = open(args.filename2.name, "r")

line = 0
iter = 0

for line in myCoverage2:
        lineData = re.split(r'\s+', line)
        if(iter % window == 0):
                coordinates2.append(lineData[1])
                # add plotting point to smoothing window
                smooth2.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage2.append(statistics.mean(smooth2))
                smooth2 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                #referenceName = lineData[0]
                smooth2.append(int(lineData[2]))
        else:
                smooth2.append(int(lineData[2]))
        iter = iter + 1

myCoverage2.close()

iter = 0

myCoverage3 = open(args.filename3.name, "r")

for line in myCoverage3:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates3.append(lineData[1])
                # add plotting point to smoothing window
                smooth3.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage3.append(statistics.mean(smooth3))
                smooth3 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth3.append(int(lineData[2]))
        else:
                smooth3.append(int(lineData[2]))
        iter = iter + 1

myCoverage3.close()

iter = 0

myCoverage4 = open(args.filename4.name, "r")

for line in myCoverage4:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates4.append(lineData[1])
                # add plotting point to smoothing window
                smooth4.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage4.append(statistics.mean(smooth4))
                smooth4 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth4.append(int(lineData[2]))
        else:
                smooth4.append(int(lineData[2]))
        iter = iter + 1

myCoverage4.close()

iter = 0

myCoverage5 = open(args.filename5.name, "r")

for line in myCoverage5:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates5.append(lineData[1])
                # add plotting point to smoothing window
                smooth5.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage5.append(statistics.mean(smooth5))
                smooth5 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth5.append(int(lineData[2]))
        else:
                smooth5.append(int(lineData[2]))
        iter = iter + 1

myCoverage5.close()

iter = 0

myCoverage6 = open(args.filename6.name, "r")

for line in myCoverage6:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates6.append(lineData[1])
                # add plotting point to smoothing window
                smooth6.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage6.append(statistics.mean(smooth6))
                smooth6 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth6.append(int(lineData[2]))
        else:
                smooth6.append(int(lineData[2]))
        iter = iter + 1

myCoverage6.close()

iter = 0

snpCoords1 = []
snpDepth1 = []

mySnps1 = open(args.snpFile1.name, "r")

for line in mySnps1:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords1.append(int(lineData[1]))
               snpDepth1.append(int(lineData[5]))
        iter = iter + 1

mySnps1.close()

iter = 0

snpCoords2 = []
snpDepth2 = []

mySnps2 = open(args.snpFile2.name, "r")

for line in mySnps2:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords2.append(int(lineData[1]))
               snpDepth2.append(int(lineData[5]))
        iter = iter + 1

mySnps2.close()

iter = 0

snpCoords3 = []
snpDepth3 = []

mySnps3 = open(args.snpFile3.name, "r")

for line in mySnps3:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords3.append(int(lineData[1]))
               snpDepth3.append(int(lineData[5]))
        iter = iter + 1

mySnps3.close()

iter = 0

snpCoords4 = []
snpDepth4 = []

mySnps4 = open(args.snpFile4.name, "r")

for line in mySnps4:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords4.append(int(lineData[1]))
               snpDepth4.append(int(lineData[5]))
        iter = iter + 1

mySnps4.close()

iter = 0

snpCoords5 = []
snpDepth5 = []

mySnps5 = open(args.snpFile5.name, "r")

for line in mySnps5:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords5.append(int(lineData[1]))
               snpDepth5.append(int(lineData[5]))
        iter = iter + 1

mySnps5.close()

iter = 0

snpCoords6 = []
snpDepth6 = []

mySnps6 = open(args.snpFile6.name, "r")

for line in mySnps6:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords6.append(int(lineData[1]))
               snpDepth6.append(int(lineData[5]))
        iter = iter + 1

mySnps6.close()


myXticks = []
iter = 0
# Generate X-axis labels
for count in coordinates1:
        if(iter % 200 == 0):
                myXticks.append(int(count) - 1)
        iter = iter + 1

myYticks = []

iter = 0
for count in range(0, max(coverage1) + 1000):
        if(iter % 2000 == 0):
                myYticks.append(count)
        iter = iter + 1

iter = 0

tabColors = mcolors.TABLEAU_COLORS

colors=[tabColors['tab:blue'], tabColors['tab:orange'], tabColors['tab:cyan'], tabColors['tab:red'], tabColors['tab:green'], tabColors['tab:pink']]

##print(tabColors)

label1 = shortTitle1t[0] + " MiSeq"
label2 = shortTitle1t[0] + " iSeq"
label3 = shortTitle2t[0] + " MiSeq"
label4 = shortTitle2t[0] + " iSeq"
label5 = shortTitle3t[0] + " MiSeq"
label6 = shortTitle3t[0] + " iSeq"


## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(14,8.5), gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Plot Genome Coverage as line plot
axes.plot(coordinates1, coverage1, label=label1, c=colors[0], linewidth=2)
axes.plot(coordinates2, coverage2, label=label2, c=colors[1])
axes.plot(coordinates3, coverage3, label=label3, c=colors[2])
axes.plot(coordinates4, coverage4, label=label4, c=colors[3])
axes.plot(coordinates5, coverage5, label=label5, c=colors[4])
axes.plot(coordinates6, coverage6, label=label6, c=colors[5])
axes.scatter(snpCoords1, snpDepth1, marker='.', c=colors[0])
axes.scatter(snpCoords2, snpDepth2, marker='.', c=colors[1])
axes.scatter(snpCoords3, snpDepth3, marker='.', c=colors[2])
axes.scatter(snpCoords4, snpDepth4, marker='.', c=colors[3])
axes.scatter(snpCoords5, snpDepth5, marker='.', c=colors[4])
axes.scatter(snpCoords6, snpDepth6, marker='.', c=colors[5])
axes.legend()
axes.set_xticks(myXticks, myXticks, rotation='vertical', size=12)
axes.set_yticks(myYticks, myYticks, size=12)
axes.set_title("Norovirus Samples " + shortTitle1t[0] + ", " + shortTitle2t[0] + ", and " + shortTitle3t[0] + " Coverage", size=16)
axes.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates', size=14)
axes.set_ylabel('Coverage (X) at Position', size=14)
#axes.margins(0.2)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/tabColSNPs_w' + str(window) + '_' + shortTitle1t[0] + '_' + shortTitle2t[0] + '_' + shortTitle3t[0] + '_to_' + referenceName + '.png')

