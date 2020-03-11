# A tutorial for implementing the Statistical Coupling Analysis with python (pySCA)

**See https://ranganathanlab.gitlab.io/pySCA/ for instructions on how to install pySCA.**

pySCA relies on a number of dependencies in order to run. For the purpose of this portfolio, this directory contains only HTML versions of the Jupyter Notebook.


### In part I:

    https://htmlpreview.github.io/?https://github.com/thuynnguyen/grad-projects/blob/master/pySCA_tutorial/general_pySCA_usage_pt1.html. 


The Jupyter Notebook helps biologists that are novices in computer science run the core scripts that execute the calculations of SCA on a multiple sequence alignment (MSA) for a protein of interest, all without using the command line. An MSA is a plain text file (e.g. FASTA) containing the amino acid sequences of a single protein from different species. In this tutorial, I am using an MSA containing sequences of the protein, Dihydrofolate Reductase, from thousands of species.

After executing the core python scripts to analyze the MSA, the notebook then helps guide the user to evaluating whether or not their MSA is appropriate for analysis with SCA.

### In part II

    http://htmlpreview.github.io/?https://github.com/thuynnguyen/grad-projects/blob/master/pySCA_tutorial/general_pySCA_usage_pt2.html 

In the second notebook (see rendering in link above) I guide the user through how to define the positions in the protein that are co-evolving from the spectral decomposition of the matrix of pair-wise positional correlations. I show how phylogenetic and functional annotations can be used to highlight groups co-evolving pairs of sequences.  
    
I used Markdown commentary in this tutorial to guide the user through the purpose, functions, and process of each pre-liminary pySCA script and step of the analysis.


### How can this tutorial be improved?
Because phylogenetic and functional annotations are specific to the model system and the alignment, the user must be familiar with using python in order to parse these annotations in Part II of the tutorial. Otherwise, it will be diffifcult for the user to evaluate how correlations between sequences could contribute to the signal of positional co-evolution. I can improve this tutorial by developing a parser that takes a user input of the format of the annotation and then uses this input to automatically parse the sequences according to phylogeny or function. 
