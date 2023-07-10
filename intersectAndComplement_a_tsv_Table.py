#!/usr/bin/python

import sys, os.path, argparse, re, logging, warnings, csv, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

## Function: A closure for .tsv or .csv extension checking

def tsv_check(expected_ext1, expected_ext2, expected_ext3, openner):
    def extension(filename):
        if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2)):
            raise ValueError()
        return openner(filename)
    return extension

## Function: Intersection (default)

def Intersection(list1, list2):
    return sorted(set(list1).intersection(list2))

## Function: Union 

def Union(list1, list2):
    return list(set().union(list1, list2))

## Function: Complement

def Compl(list1, intersec):
    snpSet1 = set(list1)
    snpSet2 = set(intersec)
    return list(snpSet1 - snpSet2)

logger = logging.getLogger("intersectAndComplement_a_tsv_Table.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Find intersection of SNPs positions from two plain text, single-column files', usage="intersectAndComplement_a_tsv_Table.py tableFile.tsv.txt")

parser.add_argument("tableFile1", type=tsv_check('.tab', '.tsv', '.csv', argparse.FileType('r')))

parser.add_argument('--union', '-u', default='N', choices=['Y','N'], help="Calculate union instead of intersection for SNP positions?")

parser.add_argument('--outputType', '-o', default='S', choices=['C', 'S'], help="'S' = output summary of intersection/union or 'C' = output single column")

args = parser.parse_args()

Raw = []
AdapterRem = []
BBDuk = []
FastP = []
Seqpurge = []
Skewer = []
Trimmomatic = []


with open(args.tableFile1.name, 'r') as table1:
    csvTable1 = csv.DictReader(table1, delimiter="\t")
    header1 = list(list(csvTable1)[0].keys())
    table1.seek(0)
    for row in csvTable1:
        for (key, val) in row.items():
            if(re.search(r'raw', key)):
                Raw.append(val)
            if(re.search(r'adapterRem', key)):
                AdapterRem.append(val)
            if(re.search(r'BBDuk', key)):
                BBDuk.append(val)
            if(re.search(r'FastP', key)):
                FastP.append(val)
            if(re.search(r'Seqpurge', key)):
                Seqpurge.append(val)
            if(re.search(r'Skewer', key)):
                Skewer.append(val)
            if(re.search(r'Trimmomatic', key)):
                Trimmomatic.append(val)

Intersec = []
Unionsec = []

if(args.union == 'Y'):
    Unionsec = Union(Raw, AdapterRem)
else:
    Intersec = Intersection(Raw, AdapterRem)

if(args.union == 'Y'):
    Unionsec = Union(Unionsec, BBDuk)
else:
    Intersec = Intersection(Intersec, BBDuk)

if(args.union == 'Y'):
    Unionsec = Union(Unionsec, FastP)
else:
    Intersec = Intersection(Intersec, FastP)

if(args.union == 'Y'):
    Unionsec = Union(Unionsec, Seqpurge)
else:
    Intersec = Intersection(Intersec, Seqpurge)

if(args.union == 'Y'):
    Unionsec = Union(Unionsec, Skewer)
else:
    Intersec = Intersection(Intersec, Skewer)



if(args.outputType == 'C'):
    if(args.union == 'Y'):
        for item in Unionsec:
            print(item)
    else:
        for item in Intersec:
            print(item)
elif(args.outputType == 'S'):
    if(args.union == 'Y'):
        print("SNPs union contains {} positions.".format(str(len(Unionsec))))
    else:
        print("SNPs intersection contains {} positions.".format(str(len(Intersec))))
    
