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

parser.add_argument('filename11', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename12', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename13', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('filename14', type=ext_check('.bed', argparse.FileType('r')))

#parser.add_argument('filename15', type=ext_check('.bed', argparse.FileType('r')))

#parser.add_argument('filename16', type=ext_check('.bed', argparse.FileType('r')))

parser.add_argument('--window', '-w', default='1', type=int)

args = parser.parse_args()

# echo input file name
print(args.filename1.name)
print(args.filename3.name)
print(args.filename5.name)
print(args.filename7.name)
print(args.filename9.name)
print(args.filename11.name)
print(args.filename13.name)
#print(args.filename15.name)

myTitle1 = re.split(r'[\/]', args.filename1.name)
myTitle2 = re.split(r'[\/]', args.filename3.name)
myTitle3 = re.split(r'[\/]', args.filename5.name)
myTitle4 = re.split(r'[\/]', args.filename7.name)
myTitle5 = re.split(r'[\/]', args.filename9.name)
myTitle6 = re.split(r'[\/]', args.filename11.name)
myTitle7 = re.split(r'[\/]', args.filename13.name)
#myTitle8 = re.split(r'[\/]', args.filename15.name)

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

shortTitle6 = re.split(r'[\.]', myTitle6[len(myTitle6) - 1])
shortTitle6t = shortTitle6[0].split(sep='_S')

shortTitle7 = re.split(r'[\.]', myTitle7[len(myTitle7) - 1])
shortTitle7t = shortTitle7[0].split(sep='_S')

#shortTitle8 = re.split(r'[\.]', myTitle8[len(myTitle8) - 1])
#shortTitle8t = shortTitle8[0].split(sep='_S')


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

coordinates11 = []
coordinates12 = []

coverage11 = []
coverage12 = []

coordinates13 = []
coordinates14 = []

coverage13 = []
coverage14 = []

#coordinates15 = []
#coordinates16 = []

#coverage15 = []
#coverage16 = []

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
smooth11 = []
smooth12 = []
smooth13 = []
smooth14 = []
#smooth15 = []
#smooth16 = []

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

iter = 0

myCoverage11 = open(args.filename11.name, "r")

for line in myCoverage11:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates11.append(iter)
                # add plotting point to smoothing window
                smooth11.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage11.append(statistics.mean(smooth11))
                smooth11 = []
        else:
                smooth11.append(int(lineData[2]))
        iter = iter + 1

myCoverage11.close()

iter = 0

myCoverage12 = open(args.filename12.name, "r")

for line in myCoverage12:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates12.append(iter)
                # add plotting point to smoothing window
                smooth12.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                coverage12.append(statistics.mean(smooth12))
                smooth12 = []
        else:
                smooth12.append(int(lineData[2]))
        iter = iter + 1

myCoverage12.close()

iter = 0

myCoverage13 = open(args.filename13.name, "r")

for line in myCoverage13:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates13.append(iter)
                # add plotting point to smoothing window
                smooth13.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                if(statistics.mean(smooth13) > 30000):
                        coverage13.append(30000)
                else:
                        coverage13.append(statistics.mean(smooth13))
                smooth13 = []
        else:
                smooth13.append(int(lineData[2]))
        iter = iter + 1

myCoverage13.close()

iter = 0

myCoverage14 = open(args.filename14.name, "r")

for line in myCoverage14:
        lineData = re.split(r'\s+', line)
        referenceName = lineData[0]
        if(iter % window == 0):
                coordinates14.append(iter)
                # add plotting point to smoothing window
                smooth14.append(int(lineData[2]))
                # take average of smoothing window, 5, 10, or 20 points
                if(statistics.mean(smooth14) > 30000):
                        coverage14.append(30000)
                else:
                        coverage14.append(statistics.mean(smooth14))
                smooth14 = []
        else:
                smooth14.append(int(lineData[2]))
        iter = iter + 1

myCoverage14.close()

iter = 0


myXticks = []
iter = 0
# Generate X-axis labels
for count in range(0, max(coordinates1)):
        if(iter % 1000 == 0):
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
if(max(topCoverage) < max(coverage11)):
        topCoverage = coverage11
if(max(topCoverage) < max(coverage12)):
        topCoverage = coverage12
if(max(topCoverage) < max(coverage13)):
        topCoverage = coverage13
if(max(topCoverage) < max(coverage14)):
        topCoverage = coverage14
#if(max(topCoverage) < max(coverage15)):
#        topCoverage = coverage15
#if(max(topCoverage) < max(coverage16)):
#        topCoverage = coverage16

#topCoverage = 30000

print(max(topCoverage))

iter = 0
for count in range(0, int(max(topCoverage)) + 1000):
        if(iter % 2000 == 0):
                myYticks.append(int(count))
        iter = iter + 1

iter = 0

tabColors = mcolors.TABLEAU_COLORS

colors=[tabColors['tab:blue'], tabColors['tab:orange'], tabColors['tab:cyan'], tabColors['tab:red'], tabColors['tab:green'], tabColors['tab:pink'], tabColors['tab:olive'], tabColors['tab:purple']]

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
label11 = shortTitle6t[0] + " MiSeq"
label12 = shortTitle6t[0] + " iSeq"
label13 = shortTitle7t[0] + " MiSeq"
label14 = shortTitle7t[0] + " iSeq"
#label15 = shortTitle8t[0] + " MiSeq"
#label16 = shortTitle8t[0] + " iSeq"

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
axes.plot(coordinates11, coverage11, label=label11, c=colors[5])
axes.plot(coordinates12, coverage12, label=label12, c=colors[5], linestyle='--', linewidth=2)
axes.plot(coordinates13, coverage13, label=label13, c=colors[6])
axes.plot(coordinates14, coverage14, label=label14, c=colors[6], linestyle='--', linewidth=2)
#axes.plot(coordinates15, coverage15, label=label15, c=colors[7])
#axes.plot(coordinates16, coverage16, label=label16, c=colors[7], linestyle='--', linewidth=2)

axes.legend()
axes.set_xticks(myXticks, myXticks, rotation='vertical', size=12)
axes.set_yticks(myYticks, myYticks, size=12)
axes.set_title("Raw Norovirus " + shortTitle1t[0] + ", " + shortTitle2t[0] + ", " + shortTitle3t[0] + ", " + shortTitle4t[0] + ", " + shortTitle5t[0] + ", " + shortTitle6t[0] + ", and " + shortTitle7t[0], size=14)
axes.set_xlabel('Reference Genome, ' + referenceName + ', Coordinates', size=14)
axes.set_ylabel('Coverage (X) at Position', size=14)
#axes.margins(0.2)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/rawSmooth_w' + str(window) + '_' + shortTitle1t[0] + '_' + shortTitle2t[0] + '_' + shortTitle3t[0] + '_' + shortTitle4t[0] + '_' + shortTitle5t[0] + '_' + shortTitle6t[0] + '_' + shortTitle7t[0] + '_to_' + referenceName + '.png')

