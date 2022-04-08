# NGS_Plot_Widgets
### Usage Statement
```python fullPlotShuffledFastq.py filepath/filename.fastq(.gz) --outputType [F/J/N]```
### Installation Instructions
Use of a Python virtual environment is recommended where the Python version should be 3.9 or higher.
After cd into the subfolder where you want the Python virtual environment, write the following command:

```virtualenv -p /apps/x86_64/python/3.9.1/bin/python ./```

The specific path for virtualenv may differ according to where your python binary is installed on your system.
Next, install the two prerequisite Python modules, Biopython and Matplotlib:

```bin/pip install biopython```
```bin/pip install matplotlib```

Then, install NGS_Plot_Widgets by git clone:

```git clone https://github.com/darlenewagner/NGS_Plot_Widgets.git```

Finally, test fullPlotShuffledFastq.py using the included test fastq.gz:

```bin/python NGS_Plot_Widgets/fullPlotShuffledFastq.py NGS_Plot_Widgets/EnterovirusD70_SRR13402413_Pairs.fastq.gz```

### Long Description
fullPlotShuffledFastq.py computes sequence lengths and average PHRED for shuffled paired reads in fastq.  
It expects a single fastq(.gz) input and outputs a Readstatistics.README.txt, a Readstatistics.json, and 
a .png image showing PHRED quality histograms for forward (R1) and reverse (R2) reads, all in a folder
named after the input filename.fastq(.gz).  Number and location of output files can be varied by --outputType.
--outputType F for full output,... J for .json only, and N for no image.
A test input file, EnterovirusD70_SRR13402413_Pairs.fastq.gz, is included for testing fullPlotShuffledFastq.py
