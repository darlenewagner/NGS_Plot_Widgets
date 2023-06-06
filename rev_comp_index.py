#!/usr/bin/python

import sys
import os.path
import argparse
import re
import logging
import warnings
import csv
import subprocess

## Reads 'SampleSheet.csv' in an Illumina NextSeq or MiSeq folder.  Takes reverse complement of index2 DNA (i7 or column 7)
## Note: This script uses regexes on lines 51 and 54 to restrict reverse complementing to 
##       sample ID rows beginning in a digit, "WOR", "RSV", or "Sample"'

## Function: A closure for .tsv or .csv extension checking

def tsv_check(expected_ext1, expected_ext2, openner):
    def extension(filename):
        if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2)):
            raise ValueError()
        return openner(filename)
    return extension

## Function: Complement DNA without Biopython

def complement(seq):
    complem = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    #for element in bases:
        #if element not in complement:
        #    print element  
    bases = [complem[base] for base in bases] 
    return(''.join(bases))


logger = logging.getLogger("rev_comp_index.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Read a comma-delimited file corresponding to a NextSeq run.  Sample name must begin in a digit, "WOR", "RSV", or "Sample"', usage="rev_comp_index.py SampleSheet.csv" )

parser.add_argument("tsvFile", type=tsv_check('.tsv', '.csv', argparse.FileType('r')))

args = parser.parse_args()


with open(args.tsvFile.name, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        row_list = list(row)
        if(len(row_list) == 0):
            continue
        if(re.search(r'(^\d+_|^M20|^75|^80|^WOR|^RSV|^Sample\d+|^AS_Untreated)', row_list[0])):
            newIdx = row[6][::-1]
            revIdx = complement(newIdx)
            if(re.search(r'(^WOR|^75|^80|Rubella|Measle)', row_list[1])):
                print(f'{row[0]},{row[1]},{row[2]},{row[3]},{row[4]},{row[5]},{row[6]},{revIdx},{row[8]},{row[9]}, ')
            elif(re.search(r'AS_Untreated', row_list[0])):
                print(f'{row[0]},{row[1]},{row[2]},{row[3]},{row[4]},{revIdx},{revIdx}, ')
            else:
                print(f'{row[0]},{row[1]},{row[2]},{row[3]},{row[4]},{row[5]},{row[6]},{revIdx},{row[8]}, ')
        else:
            i = 0
            while(i < len(row)):
                if(i < (len(row) - 1)):
                    print(f'{row[i]},', end="")
                else:
                    print(f'{row[i]} ')
                i += 1


