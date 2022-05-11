import pandas as pd

from pathlib import Path
import argparse, logging, warnings, json, os, re

## Function: A closure for file extension checking, for this script, .json

def ext_check(expected_ext, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext) or filename.lower().endswith(expected_ext3)):
                        raise ValueError()
                return openner(filename)
        return extension

## A script for displaying read1_adapter_counts and read2_adapter_counts from the .json output of fastP
parser = argparse.ArgumentParser(description='Retrieves items from large, nested .json by keyword (Expects single .json input)', usage="python piecemeal_parse_json.py filepath/filename.json")

parser.add_argument('filename', type=ext_check('.json', argparse.FileType('r')))

args = parser.parse_args()

print(args.filename.name)

basePath = os.getcwd()
inPath = os.path.split(args.filename.name)
myNameString = re.split(r'\/', args.filename.name)
newTitle = re.sub('_orig\.json', '', myNameString[len(myNameString) - 1])

# set path to file
p = Path(args.filename.name)

with p.open('r', encoding='utf-8') as f:
    data = json.loads(f.read())

type(data)

print(data.keys())

adapterR1 = data['adapter_cutting']['read1_adapter_counts']
adapterR2 = data['adapter_cutting']['read2_adapter_counts']

print(data['adapter_cutting']['read2_adapter_counts'])

halfTotalReads = int(data['read1_before_filtering']['total_reads'])

totalR1 = sum(int(sub) for sub in adapterR1.values())
totalR2 = sum(int(sub2) for sub2 in adapterR2.values())

print("%s\t%0.2f\t%0.2f" % (newTitle, totalR1*100/halfTotalReads, totalR2*100/halfTotalReads))

#df = pd.json_normalize(data)
#outPath = inPath[0] + "/" + newTitle + ".csv"
#df.to_csv(outPath, index=False, encoding='utf-8')
