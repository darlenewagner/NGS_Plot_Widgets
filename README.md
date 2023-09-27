# NGS_Plot_Widgets

### 1. Installation Instructions
A <i>Python virtual environment</i> is recommended for installation of modules. The Python version should be 3.9 or higher.
After changing directory (cd) into the subfolder where you want the Python virtual environment, write the following command:

```virtualenv -p /apps/x86_64/python/3.9.1/bin/python ./```

The specific path for virtualenv may differ according to where your python binary is installed on your system.
Next, install the two prerequisite Python modules, Biopython and Matplotlib:

```bin/pip install biopython```<br/>
```bin/pip install matplotlib```

Then, install NGS_Plot_Widgets by git clone:

```git clone https://github.com/darlenewagner/NGS_Plot_Widgets.git```

Finally, test fullPlotShuffledFastq.py using the included test fastq.gz:

```bin/python NGS_Plot_Widgets/fullPlotShuffledFastq.py NGS_Plot_Widgets/EnterovirusD70_SRR13402413_Pairs.fastq.gz```

### 2. Description and usage of fullPlotShuffledFastq.py
fullPlotShuffledFastq.py computes sequence lengths and average PHRED for shuffled paired reads in fastq.  
It expects a single fastq(.gz) input and outputs a Readstatistics.README.txt, a Readstatistics.json, and 
a .png image showing PHRED quality histograms for forward (R1) and reverse (R2) reads, all in a folder
named after the input filename.fastq(.gz).  Number and location of output files can be varied by --outputType.
--outputType F for full output,... J for .json only, and N for no image.

```python fullPlotShuffledFastq.py filepath/filename.fastq(.gz) --outputType [F/J/N]```

### 3. Description of plotBedCoverage.py
plotBedCoverage.py creates a line plot .png image from a 3-column .bedGraph file created by bedtools genomecov.  
Plotting window can be varied by entering an integer after the optional '--window' parameter.

### 4. Venn diagram plotting utility for single nucleotide polymorphisms (SNPs) positions

vennDiagramPlotColumn.py creates a 2-set Venn diagram from two input files containing unique SNPs positions.  The script 
relies upon matplotlib-venn, which is separate from matplotlib.  In the usage example below, --outputType P determines that 
a matplotlib_venn plot will be created as output, --title "my title" is a user-supplied string for annotating both the plot and its 
filename, while --plotScale [W/U] give either a weighted or unweighted Venn diagram, respectively.

```bin/python vennDiagramPlotColumn.py SC2_MiSeq_SNPs.tsv SC2_iSeq_SNPs.tsv --outputType P --title "Coronavirus Venn" --plotScale U```

The files, SC2_MiSeq_SNPs.tsv and SC2_iSeq_SNPs.tsv are based upon output from the following command line processing of .vcf files:

```bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/AO\t%INFO/DP\n' sample_1.vcf.gz | perl -ne '@F=split(/\s+/, $_); printf "%\s\t%\d\t%\s\t%\s\t%\d\t%\d\t%0.4f\n", $F[0], $F[1], $F[2], $F[3], $F[4], $F[6], $F[5]/$F[6]' >> input1.table.tsv```
