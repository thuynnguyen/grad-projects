#!/usr/bin/env python

"""
Adapted from a Jupter Notebook to a .py script to be executed in the command line. 

Inputs: a fastq file, Q-score filters each read, identifies and counts mutants in each sample. 
Outputs: a txt file with a list of amino acid mutants and one with nucleic acid mutants 

**USAGE**
	python countingAlleles.py fastqFilePathName -sl sublibraryID -o outputPathName [-q] qualityScoreThreshold 


**Keyword Arguments**
	--p	\*.fasta File of DNA sequences from Next-Generation Sequencing Run
	--sl	DHFR sublibrary ('SL1', 'SL2', 'SL3','SL4') 
	--o	Path for output file containing a list of amino acid mutants and one with nucleic acid mutants
	--q	Minimum Quality Score threshold. Default:20 (1/100 chance of incorrect base call)

**Example**::

  python countingAlleles.py path2fasta.fasta --sl sublibraryID --o path2output --q 30 

By Thuy N. Nguyen on 20191031
"""

#import modules 	
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from itertools import islice

parser = argparse.ArgumentParser()
parser.add_argument('p', help = '\*.fasta File of DNA sequences from Next-Generation Sequencing Run')
parser.add_argument('--sl',help = 'DHFR sublibrary: SL1, SL2, SL3,SL4')
parser.add_argument('--o',help = 'Path for output file containing a list of amino acid mutants and one with nucleic acid mutants')
parser.add_argument('--q',default = '20',help = 'Minimum Quality Score threshold. Default:20 (1/100 chance of incorrect base call)')
args = parser.parse_args()


#define core functions
def fileNames(line):
    spLine = line.split('\t')
    fwdFileName = spLine[0]
    sl_id = spLine[1]
    outName = spLine[2].strip('\n')
    return fwdFileName, sl_id, outName

def qscore_filter(qscores,qthreshold = 20):
    low_quality = 'False'
    for qscore in qscores: 
        if (ord(qscore) -33) <= qthreshold: #ord() converts ASCII encoding to numerical qscore
            low_quality = 'True'
    return low_quality

def return_stats():
    print( 'read_total: ' + str(read_total)) 
    print( 'read_fail: ' + str(read_fail))
    print( 'read_pass: '+ str(read_pass))
    print( 'read_unassigned: '+ str(read_unassigned))
    print( 'qscore_used was: '+ str(qthreshold))
    print( str(percent_pass)+'% reads pass qscore filter threshold - '+str(qthreshold))

def trim_read(read,coding_region_index):
    trim = read[coding_region_index[0]:coding_region_index[1]]         
    return trim

def translate_seq(seq):
    seq = Seq(seq,generic_dna)
    seq_translate = str(seq.translate().strip())
    return seq_translate


def id_sublibrary(nuc_seq):
    sl_match = 'na'
    for sl in ['sl1','sl2','sl3','sl4']:
        barcode = nuc_seq[barcode_ix[sl][0]:barcode_ix[sl][1]]
        barcode_ref = barcode_region[sl]
        if barcode == barcode_ref:
            sl_match = sl
    return sl_match

def identifyMutant(nuc_seq,nuc_ref,seq_len):
    mutants = []
    nuc_list = []
    
    aa_seq = translate_seq(nuc_seq)
    aa_ref = translate_seq(nuc_ref)
        
    def fwdMut():
        for ix,aa in enumerate(aa_seq):        
            codonIX = np.multiply(ix,3)
            if nuc_seq[codonIX:(codonIX+3)] != nuc_ref[codonIX:(codonIX+3)]:
                mutants.append(aa_ref[ix]+str(aa_pos[sl_id][ix])+aa)
                nuc_list.append(str(aa_pos[sl_id][ix])+nuc_seq[codonIX:(codonIX+3)])
                
    #compare reference and sequence read
    if aa_seq == aa_ref:
        nuc_list.append('WT') 
        mutants.append('WT')
    else:
        fwdMut()    
    return mutants, nuc_list

def record_mut_counts(mutant_dict):
    
    #for each fastq file, make a fresh mutant_counts dictionary
    #counting mutants from parsing the mutant dictionary 

    mutant_counts =  {'sl1':{},'sl2':{},'sl3':{},'sl4':{},'cl':{}}

    def count_mutant(mutant,sl):
        if mutant in mutant_counts[sl].keys():
            mutant_counts[sl][mutant] +=1
        else:
            mutant_counts[sl][mutant] = 1

    for sl in mutant_dict.keys(): #iterate through each sub-library
        for mut_list in mutant_dict[sl]:
            if len(mut_list) == 1:
                count_mutant(mut_list[0],sl)
            elif len(mut_list) >=2: 
                count_mutant('fail_multimutant',sl)
                
    return mutant_counts



def record_nuc_counts(nuc_list_dict):
    
    #for each fastq file, make a fresh mutant_counts dictionary
    #counting mutants from parsing the mutant dictionary 

    nuc_counts =  {'sl1':{},'sl2':{},'sl3':{},'sl4':{},'cl':{}}

    def count_nuc(nuc,sl):
        if nuc in nuc_counts[sl].keys():
            nuc_counts[sl][nuc] +=1
        else:
            nuc_counts[sl][nuc] = 1


    for sl in nuc_dict.keys(): #iterate through each sub-library
        for nuc_list in nuc_dict[sl]:
            if len(nuc_list) == 1:
                count_nuc(nuc_list[0],sl)
            elif len(nuc_list) >=2: 
                count_nuc('fail_multimutant',sl)
                
    return nuc_counts



def writeOutputFile(outName):

    #writes read statistics & mutant counts into one convenient .txt file  
    
    output_file = open(outName+'.txt','w')
    #write out statistics
    output_file.write('read_total:\t'+str(read_total)+'\t'+str(sl_id)+'\n')
    output_file.write('read_fail:\t'+str(read_fail)+'\n')
    output_file.write('read_pass:\t'+str(read_pass)+'\n')
    output_file.write('read_unassigned:\t'+str(read_unassigned)+'\n')
    output_file.write('qscore_used was:\t'+str(qthreshold)+'\n')
    output_file.write(str(percent_pass)+'\t% reads pass qscore filter threshold -\t'+str(qthreshold)+'\n')

    #write out mutant counts 
    for sl in mutant_counts.keys():
        for key in mutant_counts[sl].keys():
            output_file.write(sl+'\t'+key+'\t'+str(mutant_counts[sl][key])+'\n')

    output_file.close()
    
def writeOutputFile_nuc(outName):

    #writes read statistics & mutant counts into one convenient .txt file  
    
    output_file = open(outName+'nuc'+'.txt','w')
    #write out statistics
    output_file.write('read_total:\t'+str(read_total)+'\t'+str(sl_id)+'\n')
    output_file.write('read_fail:\t'+str(read_fail)+'\n')
    output_file.write('read_pass:\t'+str(read_pass)+'\n')
    output_file.write('read_unassigned:\t'+str(read_unassigned)+'\n')
    output_file.write('qscore_used was:\t'+str(qthreshold)+'\n')
    output_file.write(str(percent_pass)+'\t% reads pass qscore filter threshold -\t'+str(qthreshold)+'\n')

    #write out mutant counts 
    for sl in nuc_counts.keys():
        for key in nuc_counts[sl].keys():
            output_file.write(sl+'\t'+key+'\t'+str(nuc_counts[sl][key])+'\n')

    output_file.close()
    
###### INITIALIZE ############
# WITH WT AMPLICON SEQUENCE FROM 4N TO 4N AND INDICIES OF REGIONS IN READ YOU CARE ABOUT 
# sequence of the wt amplicon from 4N to 4N 
# see 190218_ampliconCodingIX and 190808_ampliconDesign

wt_ref = {}
wt_ref['sl1'] = 'NNNNACTTTAATAATGAGATATACCATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCNNNN'
wt_ref['sl2'] = 'NNNNCACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACNNNN'
wt_ref['sl3'] = 'NNNNGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCNNNN'
wt_ref['sl4'] = 'NNNNGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAACAGGCGTCGACAAGNNNN'

amplicon_length = {}  #expecation, for filtering reads that are shorter/longer than you expect
for sl in ['sl1','sl2','sl3','sl4']:
    amplicon_length[sl] = len(wt_ref[sl])

#indices of amplicon where the coding region of interest begins & ends in FWD READ
coding_ix = {}
coding_ix['sl1'] = [25,25+120] 
coding_ix['sl2'] = [23,23+120]
coding_ix['sl3'] = [22,22+120]
coding_ix['sl4'] = [26,26+120]

#assigning barcodes for each SL 
barcode_ix = {}
barcode_ix['sl1'] = [4,25]
barcode_ix['sl2'] = [4,23]
barcode_ix['sl3'] = [4,22]
barcode_ix['sl4'] = [4,26]

#pulling out the wt sequence of the coding regions of interest
coding_region = {}
barcode_region = {}

for sl in ['sl1','sl2','sl3','sl4']:
    coding_region[sl] = wt_ref[sl][coding_ix[sl][0]:coding_ix[sl][1]]
    barcode_region[sl] = wt_ref[sl][barcode_ix[sl][0]:barcode_ix[sl][1]]
    #OUT PUT AMINO ACID SEQUENCES 
    print( sl+' sequence: ' + translate_seq(coding_region[sl])), (sl+ ' barcode: '+barcode_region[sl])

#these cover positions 1-40 in DHFR
#sl1_aa_pos = np.arange(1,40+1) #make a list from 1:40 the indicies of this will map to indicies in coding region

aa_pos = {}
aa_pos['sl1'] = np.arange(1,40+1) 
aa_pos['sl2'] = np.arange(41,80+1)
aa_pos['sl3'] = np.arange(81,120+1)
aa_pos['sl4'] = np.arange(121,160+1)

if args.p:
    fastqFilePathName = args.p
    if fastqFilePathName.split('.')[-1] == 'fastq':
        qthreshold = int(args.q)
        sample_sl_id = args.sl
        outName = args.o
        
        # the next step is to open a fastq file and go through each read
        print('Start reading file: %s' % fastqFilePathName)
        openfastq = open(fastqFilePathName,'rU')
        fastqpath = openfastq.readlines()
        totalFileReads = len(fastqpath)/4
        openfastq.close()
        mutant_dict = {'sl1':[],'sl2':[],'sl3':[],'sl4':[]}
        nuc_dict =  {'sl1':[],'sl2':[],'sl3':[],'sl4':[]}
        
        read_unassigned = 0
        read_total = 0
        read_pass = 0
        read_fail = 0

        with open(fastqFilePathName) as fwdFile:
        #this bit is to read huge files 4 lines at a time. The second line should be coding \
        #Can implement line recognition. There is none here. Illumina + flash + ucombine outputs use the same order \
        #So no problem right now, just be aware. 
            while True:
                next_n_lines= list(islice(fwdFile, 4))
                if not next_n_lines:
                    break
                read_total += 1
                fwd_seq = next_n_lines[1].strip('\n')
                qscore_line_fwd = next_n_lines[3]

                if sample_sl_id == 'mix':
                    sl_id = id_sublibrary(fwd_seq)
                    if sl_id == 'na':
                        read_unassigned+=1
                        read_fail +=1
                        continue
                else:
                    sl_id = sample_sl_id


            #filter for amplicon sequence length
                seq_length = len(fwd_seq)
                if seq_length != amplicon_length[sl_id]:
                    read_unassigned += 1
                    read_fail +=1 
                    continue

                #trim to coding regions of interest
                fwd_seq_trim = trim_read(fwd_seq,coding_ix[sl_id])
                qscore_line_fwd_trim = trim_read(qscore_line_fwd,coding_ix[sl_id])
                #call function for qscore filter 
                fail = qscore_filter(qscore_line_fwd_trim,qthreshold=20)    

                #give a progress report...    
                if read_total%(200000) == 0:
                    print('On read %i of %i (%1.3f)' % (read_total,totalFileReads,float(read_total)/float(totalFileReads)))

                if fail == 'True':
                    read_fail+=1
                    continue

                else: 
                #function for determining mutant
                    reference_fwd = coding_region[sl_id]
                    mutants, nuc_list = identifyMutant(fwd_seq_trim,reference_fwd,seq_length)

                    mutant_dict[sl_id].append(mutants)
                    nuc_dict[sl_id].append(nuc_list)
        
        openfastq.close()
        mutant_counts = record_mut_counts(mutant_dict)
        nuc_counts = record_nuc_counts(nuc_dict)
        read_pass = float(read_total)-float(read_fail)
        percent_pass = float(read_pass)*100/float(read_total)

        return_stats()

        #writing statistics and mutant counts to a txt file 
        writeOutputFile(outName)
        writeOutputFile_nuc(outName)
        
    else: 
        print (' you must input a valid fastq file ' )