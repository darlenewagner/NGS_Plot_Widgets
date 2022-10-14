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

parser.add_argument('snpFile1', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile2', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile3', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile4', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile5', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile6', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile7', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile8', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile9', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile10', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile11', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('snpFile12', type=ext_check('.tab', argparse.FileType('r')))

#parser.add_argument('snpFile13', type=ext_check('.tab', argparse.FileType('r')))

#parser.add_argument('snpFile14', type=ext_check('.tab', argparse.FileType('r')))

#parser.add_argument('snpFile15', type=ext_check('.tab', argparse.FileType('r')))

#parser.add_argument('snpFile16', type=ext_check('.tab', argparse.FileType('r')))

parser.add_argument('--window', '-w', default='1', type=int)

args = parser.parse_args()

# echo input file name
print(args.snpFile1.name)
print(args.snpFile3.name)
print(args.snpFile5.name)
print(args.snpFile7.name)
print(args.snpFile9.name)
print(args.snpFile11.name)
#print(args.snpFile13.name)
#print(args.snpFile15.name)

myTitle1 = re.split(r'[\/]', args.snpFile1.name)
myTitle2 = re.split(r'[\/]', args.snpFile3.name)
myTitle3 = re.split(r'[\/]', args.snpFile5.name)
myTitle4 = re.split(r'[\/]', args.snpFile7.name)
myTitle5 = re.split(r'[\/]', args.snpFile9.name)
myTitle6 = re.split(r'[\/]', args.snpFile11.name)
#myTitle7 = re.split(r'[\/]', args.snpFile13.name)
#myTitle8 = re.split(r'[\/]', args.snpFile15.name)

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

#shortTitle7 = re.split(r'[\.]', myTitle7[len(myTitle7) - 1])
#shortTitle7t = shortTitle7[0].split(sep='_S')

#shortTitle8 = re.split(r'[\.]', myTitle8[len(myTitle8) - 1])
#shortTitle8t = shortTitle8[0].split(sep='_S')


csvRow = []
forwardLen = []
reverseLen = []


referenceName = ''

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

iter = 0

snpCoords7 = []
snpDepth7 = []

mySnps7 = open(args.snpFile7.name, "r")

for line in mySnps7:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords7.append(int(lineData[1]))
               snpDepth7.append(int(lineData[5]))
        iter = iter + 1

mySnps7.close()

iter = 0

snpCoords8 = []
snpDepth8 = []

mySnps8 = open(args.snpFile8.name, "r")

for line in mySnps8:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords8.append(int(lineData[1]))
               snpDepth8.append(int(lineData[5]))
        iter = iter + 1

mySnps8.close()

iter = 0

snpCoords9 = []
snpDepth9 = []

mySnps9 = open(args.snpFile9.name, "r")

for line in mySnps9:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords9.append(int(lineData[1]))
               snpDepth9.append(int(lineData[5]))
        iter = iter + 1

mySnps9.close()

iter = 0

snpCoords10 = []
snpDepth10 = []

mySnps10 = open(args.snpFile10.name, "r")

for line in mySnps10:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords10.append(int(lineData[1]))
               snpDepth10.append(int(lineData[5]))
        iter = iter + 1

mySnps10.close()

iter = 0

snpCoords11 = []
snpDepth11 = []

mySnps11 = open(args.snpFile11.name, "r")

for line in mySnps11:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords11.append(int(lineData[1]))
               snpDepth11.append(int(lineData[5]))
        iter = iter + 1

mySnps11.close()

iter = 0

snpCoords12 = []
snpDepth12 = []

mySnps12 = open(args.snpFile12.name, "r")

for line in mySnps12:
        lineData = re.split(r'\s+', line)
        if(iter > 1):
               snpCoords12.append(int(lineData[1]))
               snpDepth12.append(int(lineData[5]))
        iter = iter + 1

mySnps12.close()


myXticks = []
iter = 0
# Generate X-axis labels
for count in range(0, max(snpCoords1)):
        if(iter % 2000 == 0):
                myXticks.append(int(count))
        iter = iter + 1

#print(myXticks)

myYticks = []

topCoverage = snpDepth1

if(max(topCoverage) < max(snpDepth2)):
        topCoverage = snpDepth2
if(max(topCoverage) < max(snpDepth3)):
        topCoverage = snpDepth3
if(max(topCoverage) < max(snpDepth4)):
        topCoverage = snpDepth4
if(max(topCoverage) < max(snpDepth5)):
        topCoverage = snpDepth5
if(max(topCoverage) < max(snpDepth6)):
        topCoverage = snpDepth6
if(max(topCoverage) < max(snpDepth7)):
        topCoverage = snpDepth7
if(max(topCoverage) < max(snpDepth8)):
        topCoverage = snpDepth8
if(max(topCoverage) < max(snpDepth9)):
        topCoverage = snpDepth9
if(max(topCoverage) < max(snpDepth10)):
        topCoverage = snpDepth10
if(max(topCoverage) < max(snpDepth11)):
        topCoverage = snpDepth11
if(max(topCoverage) < max(snpDepth12)):
        topCoverage = snpDepth12
#if(max(topCoverage) < max(coverage13)):
#        topCoverage = coverage13
#if(max(topCoverage) < max(coverage14)):
#        topCoverage = coverage14
#if(max(topCoverage) < max(coverage15)):
#        topCoverage = coverage15
#if(max(topCoverage) < max(coverage16)):
#        topCoverage = coverage16

print(max(topCoverage))

iter = 0
for count in range(0, int(max(topCoverage)) + 3000):
        if(iter % 2000 == 0):
                myYticks.append(int(count))
        iter = iter + 1

iter = 0

tabColors = mcolors.TABLEAU_COLORS

colors=[tabColors['tab:blue'], tabColors['tab:orange'], tabColors['tab:cyan'], tabColors['tab:red'], tabColors['tab:green'], tabColors['tab:pink'], tabColors['tab:olive'], tabColors['tab:purple'] ]

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
#label13 = shortTitle7t[0] + " MiSeq"
#label14 = shortTitle7t[0] + " iSeq"
#label15 = shortTitle8t[0] + " MiSeq"
#label16 = shortTitle8t[0] + " iSeq"

## A single plot in the subplots.  Padding of 15% on bottom margin and 10% for the other three margins.
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(14,8.5), gridspec_kw=dict(left=0.1, right=0.9, bottom=0.15, top=0.9))

## Plot Genome Coverage as line plot

axes.scatter(snpCoords1, snpDepth1, label=label1, marker='.', c=colors[0], s=100)
axes.scatter(snpCoords2, snpDepth2, label=label2, marker='+', c=colors[0], s=85)
axes.scatter(snpCoords3, snpDepth3, label=label3, marker='.', c=colors[1], s=100)
axes.scatter(snpCoords4, snpDepth4, label=label4, marker='+', c=colors[1], s=85)
axes.scatter(snpCoords5, snpDepth5, label=label5, marker='.', c=colors[2], s=100)
axes.scatter(snpCoords6, snpDepth6, label=label6, marker='+', c=colors[2], s=85)
axes.scatter(snpCoords7, snpDepth7, label=label7, marker='.', c=colors[3], s=100)
axes.scatter(snpCoords8, snpDepth8, label=label8, marker='+', c=colors[3], s=85)
axes.scatter(snpCoords9, snpDepth9, label=label9, marker='.', c=colors[4], s=100)
axes.scatter(snpCoords10, snpDepth10, label=label10, marker='+', c=colors[4], s=85)
axes.scatter(snpCoords11, snpDepth11, label=label11, marker='.', c=colors[5], s=100)
axes.scatter(snpCoords12, snpDepth12, label=label12, marker='+', c=colors[5], s=85)

axes.legend()
axes.set_xticks(myXticks, myXticks, rotation='vertical', size=14)
axes.set_yticks(myYticks, myYticks, size=14)

#axes.set_title("SNPs of " + shortTitle1t[0] + ", " + shortTitle2t[0] + ", " + shortTitle3t[0] + ", " + shortTitle4t[0] + ", " + shortTitle5t[0] + ", and " + shortTitle6t[0], size=16)
axes.set_title("SNPs of MiSeq and iSeq SARS-CoV-2 Samples", size=28)
axes.set_xlabel('Reference Genome, MN908947, Coordinates', size=15)
axes.set_ylabel('Coverage (X) at Position', size=15)
#axes.margins(0.2)

fig.savefig('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/rawSNP_w' + '_' + shortTitle1t[0] + '_' + shortTitle2t[0] + '_' + shortTitle3t[0] + '_' + shortTitle4t[0] + '_' + shortTitle5t[0] + '_and_' + shortTitle6t[0] + '_to_MN908947.png')

