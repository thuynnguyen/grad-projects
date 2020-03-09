# Analyzing Next-Generation Sequencing Data

**countingAlleles.py** is a python script that identifies and counts DHFR mutants in a given FASTA file (a plain text document with millions of DHFR DNA sequences). With this script, I was able to count 6,360 DHFR mutants at each time point (0, 4, 8, 12, 20, and 24 hours) sampled in an experiment that measuring effects of bacterial growth rate in a saturation mutagenesis library of DHFR in two TYMS backgrounds. 

	>USAGE: python countingAlleles.py [FASTA-FILE] -o [OUTPUT] -q [QUALITY-SCORE-THRESHOLD: default = 20] -sl [DHFR Sub-library ID] 

How does **countingAlleles.py** work? 


1. Each entry in the fasta file has four lines. One is an identifier, one is the DNA sequence, the other is an accompanying ASCII sequence that encodes the Quality-Score of each base corresponding base in the sequence. The last line is a place holder "+". 

2. If 
