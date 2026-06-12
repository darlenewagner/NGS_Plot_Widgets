# NGS_Plot_Widgets

### 1. Installation Instructions
Until the Singularity-containerized version is availalbe, a <i>Python virtual environment</i> is recommended for installation of modules. The Python version should be 3.13 or higher.

Install NGS_Plot_Widgets by git clone:

```git clone https://github.com/darlenewagner/NGS_Plot_Widgets.git```

Install python3.13+ locally or point to your systems python resource as shown:

```module load python/3.13.7```

Change working directory to NGS_Plot_Widgets:

```cd NGS_Plot_Widgets```

Set up your python3.13 virtual environment inside NGS_Plot_Widgets. <i>The exact path will vary according to python installation location</i>:

```virtualenv -p /apps/x86_64/Python/3.13.7-GCCcore-14.2.0/bin/python ./```

```python3.13 -m venv venv```

Activate your newly-installed virtual environment to add python modules essential for NGS_Plot_Widgets:

```source venv/bin/activate```

Upgrade the pip installer and then install the modules, biopython, matplotlib, and pandas:

```pip install --upgrade pip```

```pip install biopython```

```pip install matplotlib```

```pip install pandas```

Finally, test fullPlotShuffledFastq.py using the included test fastq.gz:

```bin/python NGS_Plot_Widgets/fullPlotShuffledFastq.py NGS_Plot_Widgets/EnterovirusD70_SRR13402413_Pairs.fastq.gz```

When finished using scripts in NGS_Plot_Widgets, deactivate the virtual environment:

```deactivate```

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
