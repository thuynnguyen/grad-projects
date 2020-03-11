# Analyzing Next-Generation Sequencing Data

### The PDF in this directory visually describes the computational workflow for analyzing the data I collected from experiments for my thesis project. 

### How does **countingAlleles.py** work? 


**countingAlleles.py** is a python script that identifies and counts DHFR mutants in a given FASTA file (a plain text document with millions of DHFR DNA sequences). With this script, I was able to count 6,360 DHFR mutants at each time point (0, 4, 8, 12, 20, and 24 hours) sampled in an experiment that measuring effects of bacterial growth rate in a saturation mutagenesis library of DHFR in two TYMS backgrounds. 

	>USAGE: python countingAlleles.py [FASTA-FILE] -o [OUTPUT] -q [QUALITY-SCORE-THRESHOLD: default = 20] -sl [DHFR Sub-library ID] 



1. Each entry in the fasta file has four lines. One is an identifier, one is the DNA sequence, the other is an accompanying ASCII sequence that encodes the Quality-Score of each base corresponding base in the sequence. The last line is a place holder "+". 

2. This script is specific to how the DHFR library is constructed, as four sub-libraries containing mutations for 40 positions in DHFR at a time. Each base call in the sequence read has an accompanying probaility of being incorrectly identified as an A,T,G, or C. 

3. If a single base call in the read has a Quality Score lower than the threshold, the entire read is considered a low-quality data point and removed from further analysis.

		This probability is encoded as a Quality-Score (see https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm). The default Quality-Score is 20, equivalent to -log(1/100). In other words, this means that a base call with a Quality Score of 20 has a 1 in 100 chance of being incorrect.  

4. The output is a tab-delimited file containing a list of mutants and the size of the population in the sample.

### How can **countingAlleles.py** improve? 

1. I would make this modular for any output of a next-generation sequencing run where the objective is to count all variants in a population. Currently, this script is highly specific to the model system of my thesis project, Dihydrofolate Reductase. The reference sequence is hard-coded into the program as are the indicies for the coding-regions of each sub-library. 

2. To improve the efficiency for processing a single FASTQ file on a local drive can take up to two hours for a 10-15 GB file. With a high-performance computing cluster, it takes approximately 20 minutes.


### The output of countingalleles.py is analyzed further in Jupyter Notebook. 

In this directory is a HTML file that is the exported version of the notebook. In this notebook, my data analysis process is driven by two major questions: 

1. Are the data high-fidelity? 
	- What proportion of reads were filtered out as low-quality reads in each sample? 
	- What is the coverage of each mutant? In the library, we aim for > 1000 reads for each mutant. The Lower the counts for a particular mutant, the noiser the fit for growth rate will be.  
	- Are all mutants present? 

2. How well are the fits of growth rate from linear regression of relative frequency over time? 

3. How comparable are the growth rates effects of DHFR mutants with known biochemical activity? 

4. How does the mutating the background thymidylate synthase to an inactive form affect the growth rates for each mutant in the library? 

