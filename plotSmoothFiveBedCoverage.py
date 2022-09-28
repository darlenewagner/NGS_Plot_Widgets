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

parser.add_argument('filename3', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename4', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename5', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename6', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename7', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename8', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename9', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename10', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('--window', '-w', default='1', type=int)

args = parser.parse_args()

# echo input file name
print(args.filename1.name)
print(args.filename3.name)
print(args.filename5.name)
print(args.filename7.name)
print(args.filename9.name)

myTitle1 = re.split(r'[\/]', args.filename1.name)
myTitle2 = re.split(r'[\/]', args.filename3.name)
myTitle3 = re.split(r'[\/]', args.filename5.name)
myTitle4 = re.split(r'[\/]', args.filename7.name)
myTitle5 = re.split(r'[\/]', args.filename9.name)

#print(myTitle[len(myTitle) - 1])

shortTitle1 = re.split(r'[\.]', myTitle1[len(myTitle1) - 1])
shortTitle1t = shortTitle1[0].split(sep='_S')

shortTitle2 = re.split(r'[\.]', myTitle2[len(myTitle2) - 1])
shortTitle2t = shortTitle2[0].split(sep='_S')

shortTitle3 = re.split(r'[\.]', myTitle3[len(myTitle3) - 1])
shortTitle3t = shortTitle3[0].split(sep='_S')

shortTitle4 = re.split(r'[\.]', myTitle4[len(myTitle4) - 1])
shortTitle4t = shortTitle4[0].split(sep='_S')

shortTitle5 = re.split(r'[\.]', myTitle5[len(myTitle5) - 1])
shortTitle5t = shortTitle5[0].split(sep='_S')


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

coordinates7 = []
coordinates8 = []

coverage7 = []
coverage8 = []

coordinates9 = []
coordinates10 = []

coverage9 = []
coverage10 = []

referenceName = ''

iter = 0
# avoid plotting every point, plot average of every 5th, 10th, or 20th point
smooth1 = []
smooth2 = []
smooth3 = []
smooth4 = []
smooth5 = []
smooth6 = []
smooth7 = []
smooth8 = []
smooth9 = []
smooth10 = []

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
                coordinates1.append(iter)
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

#print(coordinates1)

myCoverage1.close()

myCoverage2 = open(args.filename2.name, "r")

line = 0
iter = 0

for line in myCoverage2:
        lineData = re.split(r'\s+', line)
        if(iter % window == 0):
                coordinates2.append(iter)
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
                coordinates3.append(iter)
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
                coordinates4.append(iter)
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
                coordinates5.append(iter)
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
                coordinates6.append(iter)
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

myCoverage7 = open(args.filename7.name, "r")

for line in myCoverage7:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates7.append(iter)
                # add plotting point to smoothing window
                smooth7.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage7.append(statistics.mean(smooth7))
                smooth7 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth7.append(int(lineData[2]))
        else:
                smooth7.append(int(lineData[2]))
        iter = iter + 1

myCoverage7.close()

iter = 0

myCoverage8 = open(args.filename8.name, "r")

for line in myCoverage8:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates8.append(iter)
                # add plotting point to smoothing window
                smooth8.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage8.append(statistics.mean(smooth8))
                smooth8 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth8.append(int(lineData[2]))
        else:
                smooth8.append(int(lineData[2]))
        iter = iter + 1

myCoverage8.close()

iter = 0

myCoverage9 = open(args.filename9.name, "r")

for line in myCoverage9:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates9.append(iter)
                # add plotting point to smoothing window
                smooth9.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage9.append(statistics.mean(smooth9))
                smooth9 = []
        elif(iter == 1):
                #collect non-plotted points for smoothing average
                smooth9.append(int(lineData[2]))
        else:
                smooth9.append(int(lineData[2]))
        iter = iter + 1

myCoverage9.close()

iter = 0

myCoverage10 = open(args.filename10.name, "r")

for line in myCoverage10:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates10.append(iter)
                # add plotting point to smoothing window
                smooth10.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage10.append(statistics.mean(smooth10))
                smooth10 = []
        else:
                smooth10.append(int(lineData[2]))
        iter = iter + 1

myCoverage10.close()

iter = 0


myXticks = []
iter = 0
# Generate X-axis labels
for count in range(0, max(coordinates1)):
        if(iter % 200 == 0):
                myXticks.append(int(count))
        iter = iter + 1

#print(myXticks)

myYticks = []

topCoverage = coverage1

if(max(topCoverage) < max(coverage2)):
        topCoverage = coverage2
if(max(topCoverage) < max(coverage3)):
        topCoverage = coverage3
if(max(topCoverage) < max(coverage4)):
        topCoverage = coverage4
if(max(topCoverage) < max(coverage5)):
        topCoverage = coverage5
if(max(topCoverage) < max(coverage6)):
        topCoverage = coverage6
if(max(topCoverage) < max(coverage7)):
        topCoverage = coverage7
if(max(topCoverage) < max(coverage8)):
        topCoverage = coverage8
if(max(topCoverage) < max(coverage9)):
        topCoverage = coverage9
if(max(topCoverage) < max(coverage10)):
        topCoverage = coverage10


print(max(topCoverage))

iter = 0
for count in range(0, int(max(topCoverage)) + 1000):
        if(iter % 2000 == 0):
                myYticks.append(int(count))
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
label7 = shortTitle4t[0] + " MiSeq"
label8 = shortTitle4t[0] + " iSeq"
label9 = shortTitle5t[0] + " MiSeq"
label10 = shortTitle5t[0] + " iSeq"

## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(14,8.5), gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Plot Genome Coverage as line plot
axes.plot(coordinates1, coverage1, label=label1, c=colors[0])
axes.plot(coordinates2, coverage2, label=label2, c=colors[0], linestyle='--', linewidth=2)
axes.plot(coordinates3, coverage3, label=label3, c=colors[1])
axes.plot(coordinates4, coverage4, label=label4, c=colors[1], linestyle='--', linewidth=2)
axes.plot(coordinates5, coverage5, label=label5, c=colors[2])
axes.plot(coordinates6, coverage6, label=label6, c=colors[2], linestyle='--', linewidth=2)
axes.plot(coordinates7, coverage7, label=label7, c=colors[3])
axes.plot(coordinates8, coverage8, label=label8, c=colors[3], linestyle='--', linewidth=2)
axes.plot(coordinates9, coverage9, label=label9, c=colors[4])
axes.plot(coordinates10, coverage10, label=label10, c=colors[4], linestyle='--', linewidth=2)

axes.legend()
axes.set_xticks(myXticks, myXticks, rotation='vertical', size=12)
axes.set_yticks(myYticks, myYticks, size=12)
axes.set_title("Raw Polio " + shortTitle1t[0] + ", " + shortTitle2t[0] + ", " + shortTitle3t[0] + ", " + shortTitle4t[0] + ", and " + shortTitle5t[0], size=15)
axes.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates', size=14)
axes.set_ylabel('Coverage (X) at Position', size=14)
#axes.margins(0.2)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/rawSmooth_w' + str(window) + '_' + shortTitle1t[0] + '_' + shortTitle2t[0] + '_' + shortTitle3t[0] + '_' + shortTitle4t[0] + '_and_' + shortTitle5t[0] + '_to_' + referenceName + '.png')

