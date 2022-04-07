# NGS_Plot_Widgets
usage: python fullPlotShuffledFastq.py filepath/filename.fastq(.gz) --outputType [F/J/N]
fullPlotShuffledFastq.py computes sequence lengths and average PHRED for shuffled paired reads in fastq.  
It expects a single fastq(.gz) input and outputs a Readstatistics.README.txt, a Readstatistics.json, and 
a .png image showing PHRED quality histograms for forward (R1) and reverse (R2) reads, all in a folder
named after the filename.fastq(.gz).  Amount and location of output can be varied by --outputType.
--outputType F for full output,... J for .json only, and N for no image.
