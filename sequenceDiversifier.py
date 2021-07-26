#### Run the degenerator program

# import modules
import numpy as np
from numpy import random
import scipy.stats
import pandas as pd
import collections
from random import randint
from collections import Counter
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC
import string
import RNA
import os
import re
import subprocess
from scipy import random
import json
import string
import csv
from Bio.Alphabet import IUPAC

from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Alphabet import generic_nucleotide

# ------------------------ Variable definition (codon table) --------------------------- #

# Assign full codon table on which to operate codon selection based on user's input
codon_table_full = { 'A': ['GCT','GCC','GCA','GCG'],
                'C': ['TGT','TGC'],
                'D': ['GAT','GAC'],
                'E': ['GAA','GAG'],
                'F': ['TTT','TTC'],
                'G': ['GGT','GGC','GGA','GGG'],
                'H': ['CAT','CAC'],
                'I': ['ATT','ATC','ATA'],
                'K': ['AAA','AAG'],
                'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
                'M': ['ATG'],
                'N': ['AAT','AAC'],
                'P': ['CCT','CCC','CCA','CCG'],
                'Q': ['CAA','CAG'],
                'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
                'S': ['AGT','AGC','TCT','TCC','TCA','TCG'],
                'T': ['ACT','ACC','ACA','ACG'],
                'V': ['GTT','GTC','GTA','GTG'],
                'W': ['TGG'],
                'Y': ['TAT','TAC'],
                '*': ['TAG','TAA','TGA'] }  


# Define functions

# ---------------------------------- Input functions: START ---------------------------- #

# functions to read in user input - codon usage tables from file, protein sequence, w/t cds and RSCU thresholds

def input_expression_host():
    """ Uses raw_input to read in expression host"""
    
    host = raw_input('Enter expression host (human, mouse OR cho): ')
    
    if host == 'human':
        path = 'human_codon_usage.csv'
    if host == 'mouse':
        path = 'mouse_codon_usage.csv'
    if host == 'cho':
        path = 'cho_codon_usage.csv'
    
    return path
    
def codon_table_cutoffs():
    """ Uses raw_input to let the user define RSCU cutoff values"""
    
    rar_thr = raw_input('Enter RSCU value below which codons are discared: ')
    
    gc_thr = raw_input('Enter RSCU value below which codons not ending in G/C are discarded: ')
    
    return rar_thr, gc_thr
    
def generate_codon_table(path, codon_table_full):
    """ 
    Uses functions input_expression_host and codon_table_cutoffs to generate a codon table.
    Args: codon_table_full, dictionary with amino acids as keys, each taking a list of codons as value;
    path: the path to a codon usage table. Variable path can be generated with input_expression_host function
    """
    
    # load csv file for expression host 
    rscu = pd.read_csv(path)
    
    # Define cutoff values and create codon_table
    cutoff_val = codon_table_cutoffs()
    rar_thr = float(cutoff_val[0])
    gc_thr = float(cutoff_val[1])

    # Extract amino acid symbols (one-letter code) from .csv file
    symb = rscu['AmOneLet']

    # Initialise a dictionary with unique amino acid symbols as keys
    codon_table = {}
    for s in symb:
        if s not in codon_table:
            codon_table[s] = []
            
    # Extract parameters from .csv file
    cod = rscu['Codon']
    val = rscu['RSCU']
    gc3 = rscu['GC3']
    
    # Fill codon table with codons above thresholds
    for index in range(len(symb)):
        # Only consider codons if above 'rare threshold'
        if val[index] > rar_thr:
            # If above threshold, take all cods with rscu above gc_thr or ones ending in g/c (regardless of val) 
            if val[index] > gc_thr or gc3[index] == 'Y':
                codon_table[symb[index]].append(cod[index])
                
                          
    # If RSCU thresholds set too high, most highly expressed codons are added to codon_table
    
    # Create a tuple with (symb, cod, val), used below to take most highly expressed codon
    bund = [ (symb[i], cod[i], val[i]) for i in range(len(cod)) ]
    
    # Cycle through keys (amino acids) in codon_table
    for key in codon_table:
        # If amino acids without corresponding codons are found..
        if codon_table[key] == []:
            if key != 'W' and key != 'M':
                print('Using only most highly expressed codon for amino acid %s' % key)
            # Add most highly expressed codon from codon usage table
            pool = [ bund[x] for x in range(len(bund)) if bund[x][0] == key ]
            c = max(pool, key=lambda item:item[2])[1]
            codon_table[key].append(c)
            
    return codon_table
    

    
def input_wt_codingSeq():
    """Uses raw_input to read in w/t coding sequence - optional."""
    wt = raw_input('Optional: enter wild type coding sequence (require for CAI and GC comparisons: ')
    
    # Make sure sequence is upper case, there are no blanks or newline characters
    wt = wt.upper().replace(' ','').replace('\n','')
    
    return wt
    

# -------------------------- Input functions' definitions: END ------------------------- #


# ---------------------------------- Main functions: START ----------------------------- #

### Main functions for degenerator program


def find_identities(l):
    """
    Takes in a list and returns a dictionary with seqs as keys and positions of identical elements in list as values.
    argvs: l =  list, e.g. mat[:,x]
    """

    # the number of items in the list will be the number of unique types
    uniq = [item for item, count in collections.Counter(l).items()]
    
    # Initialise a dictionary that will hold the results
    identDict = {}
    for item in uniq:
        identDict[item] = [ x for x in range(len(l)) if l[x] == item ]
        
    return identDict                
                
                

# A function that produces a generator
# Usage:        
# alternator = alternate()      
def alternate():
    """A function that produces a generator
    Usage: alternator = alternate() """ 

    while True:
        yield 0
        yield 1


def lib_generator(my_prot, seqVars, codon_table):
    """
    Function to generate set of diversified coding sequences. 
    argvs: my_prot = amino acid sequence of the protein in question (string),
    seqVars = number of coding sequence to generate (integer),
    codon_table: codon usage table for the specified expression host (dict)
    """

   

    ### --- Preliminary operations, e.g. variables' initiation
    # N is the number of codons, same as the number of aa's in our protein sequence
    N = len(my_prot)
    
    # Initialise data structure to store CDS variants
    # use np array with dtype = string of 3 chars, as each codon is stored separately
    mat = np.zeros((seqVars,N), dtype='S3') 
    
    # create a list (codVars) with all codon variants at each position
    codVars = []
    # Fill codVars with a loop
    for aa in my_prot:
        codVars.append(codon_table[aa])
        
    # Check that dimensions of codVars and recipient data structure mat are compatible
    if len(codVars) != mat.shape[1]:
        raise ValueError('list of codon vars and matrix incompatible')   
        
    ### --- Now filling first two position of recipient mat
    
    # Pos 1, where we ensure that ATG is the starting codon
    if codVars[0] != ['ATG']:
        raise ValueError('first codon is not ATG')
    mat[:,0] = codVars[0]

    # Pos 2, where we ensure available codon variants are spread across growing sequences
    cur = 0
    while cur < seqVars:
        # Trick to cycle through codon variants and assign them to mat
        codIndex = cur - (cur/len(codVars[1]))*len(codVars[1])
        mat[cur,1] = codVars[1][codIndex]
        cur += 1  
    
    ### --- Method for comparing sequences and sort them into different sets
    
    # cycle through codVars (= codon variants, list of lists) and fill mat 
    for k in range( 2, len(codVars) ):
        
        # Initialise final and initial positions for joining codons
        # iniPos will change within the while loop
        finPos = k 
        iniPos = k - 1
        
        # Initialise data frame to hold identity dictionaries, could be a list of dicts
        idents = []
        
        
        #### From here we have an iterative logic block, whereby
        #### 1) identities of stretches with increasing length are evaluated;
        #### 2) for each position to fill (k), a list of dictionaries is returned
        
        # Create a switch to exit while loop when appropriate
        switch = 1
        
        while switch == 1:
        
            # Initialise joint Codons list
            jointCods = []
            
            # Create the set that will be submitted to find_identities() function
            for x in range(seqVars):
                jointCods.append("".join([mat[x,p] for p in range(iniPos,finPos)]))   
                
            # Deploy find_identities() on jointCods
            identDict = find_identities(jointCods)         
            
            ## Add hypothetical case where there are no identities at all
            
            # Clause to break out of while loop: if there are no identities or if first
            # codon has been reached 
            
            ### Note: we do not want the last entry in idents to be singlets
            ### UNLESS, iniPos = k - 1, i.e. there no terminal identities
            
            if iniPos == k - 1 and max([len(x) for x in identDict.values()]) == 1:
                ### We are allowed to have singlets and return
                idents.append(identDict)
                ### Exit while loop by turning switch off
                switch = 0
                
            ### We do not want to append the last comparison with singlets, so we do
            ### not append with the following condition, just turn off the switch    
            if max([len(x) for x in identDict.values()]) == 1 or iniPos == 0:
                
                ### Exit while loop by turning switch off
                switch = 0
                
            elif max([len(x) for x in identDict.values()]) > 1:
                ### Add ident dict to idents
                idents.append(identDict)
                ### Do not turn switch off but move iniPos back instead
                iniPos = iniPos - 1
                
            else:
                raise ValueError('problem with identDict evaluation or logic')
        
        
        ########### ----------------- Block 1 ends here ----------------------############
        
        # Case where only one codon is available: just fill all seqs with that
        if len(codVars[k]) == 1:
            cod = codVars[k][0] 
            for pos in range(seqVars):
                mat[pos,k] = cod
                
        else:
            ########### ----------------- Block 2 ----------------- ######################
            
            # Start with special condition, note list around list
            if len(idents) == 1:
                outerDict_temp = {}
                for key in idents[0]:
                    outerDict_temp[key] = [ idents[0][key] ]
                   
            if len(idents) > 1:
                # Extract outer idents layer in dictionary
                outerDict = idents[0]
                outerDict_temp = {}
                # Create a second dictionary with same keys as above, but empty value-lists
                for item in outerDict:
                    outerDict_temp[item] = [] 
                    
                # Sequence clustering
                
                # Cycle through all layers in idents (including 3-bp homology layer, i.e. outermost layer)
                for cur in range(1, len(idents) + 1):
                    # print 'cur is:', cur
                    index = -cur
                    # Take sequences for which there are > 1 ID's
                    
                    longestHom = [ x for x in idents[index] if len(idents[index][x]) > 1 ]
                    # Cycle through items in longestHom   
                    for item in longestHom:
				
                        # Define last triplet as key for outerDict_temp (storage Dict)
                        currentKey = item[-3:]
					
                        # Define elements corresponding to item in longestHom
                        l = [ x for x in idents[index][item] ]
					
                        # if outerDict_temp[currentKey] is empty, just append
                        if outerDict_temp[currentKey] == []:
                            outerDict_temp[currentKey].append(l)
                            
                        else:
                            
					
                            # Two cases - 1) no item in l is present in outerDict_temp[currentKey]
						
                            # all items in outerDict_temp[currentKey] as a set
                            all_items = set( x for sublist in outerDict_temp[currentKey] for x in sublist )
						
                            # if l and all_items do not share elements
                            if set(l).isdisjoint(all_items):
						
                                # Append list l to outerDict_temp
                                outerDict_temp[currentKey].append(l)
						
                            # 2) they are not disjoint
                            else:
								
                                # Cycle through sublists in dict
                                for sublist in outerDict_temp[currentKey]:
									
                                    # if a sublist shares items with l
                                    if not set(l).isdisjoint(set(sublist)):
								
                                        # Find items NOT in common
                                        dif = set(l).difference(set(sublist))
									
                                        # Cycle through them
                                        for e in dif:
									
                                            # If e is not already in all_items
                                            if e not in all_items:
											 
                                                # Append element to sublist and update all_items
                                                sublist.append(e)
                                                all_items.add(e)
                                                
            ## Special condition - if singlets are present throughout the depth of idents,
            ## these are not considered. Fix by collecting them at the end
            for key in outerDict_temp:
                if outerDict_temp[key] == []:
                    # print 'FIXING BUG'
                    outerDict_temp[key] = [idents[0][key]]
                                    
            ####################### -------------------------- ###########################                
                
            ########### Module 3 -  filling module ###################
            ### Now we work on outerDict_temp
            
            # Initialise generator/alternator
            #print 'alternator is being initialised'
            alternator = alternate() 
            #print 'alternator is on: ', alternator
            
        
            codList = [ x for x in codVars[k] ]    
            
            # Shuffle codList
            random.shuffle(codList)
	
            # Split the shuffled list into two
            x1 = codList[0:len(codList)/2]
            x2 = codList[len(codList)/2:len(codList)]
            # Then we put them together
            pool = [x1, x2]
		
            # Initialise state, only done on first iteration
            if k == 2:
                state = alternator.next()
                
            # Before starting work on this 'k', ensure state is on 0
            if state == 1:
                state = alternator.next()
            
            # Cycle through items in dictionary, each set defined as workingList
            for item in outerDict_temp:
                workingList = outerDict_temp[item]     ###### SET OF LISTS
            
                # Cycle through sublists in workingList
                for sub in workingList:
                    # On making this transition, we re-define x1 and x2
                
                    for id in sub:                     ###### SINGLE ID      
                
                        ### Within the same sub list, we can flick between x1 and x2 when they are used up
                        
                        if pool[state] == []:
                            ### Switch to codon sublist not in use
                            state = alternator.next()
                        
                            ### Reload previous state
                            pool[state - 1] = [ x for x in codVars[k] if x not in pool[state] ]
                            
                        ## Standard lines for choosing codon and appending to matrix
                        if pool[state] == []:
                            state = alternator.next()
                        cod = random.choice(pool[state])
                        mat[id,k] = cod
                    
                        pool[state].remove(cod)
                           
                    ## Here we re-define x1 and x2, based on state
                    ######### Note: now adopting a new strategy, previous cod in unitTest if nec.
                    if state == 0:
                        # x2 becomes = leftovers of x1 + x2, and state switches
                        x2 = pool[state] + pool[state - 1]
                        x1 = [ x for x in codVars[k] if x not in x2]
                        pool = [x1,x2]
                        # Now we switch the state
                        state = alternator.next()
                        # 
                    if state == 1:
                        # x2 is leftovers of x2, state remains the same
                        x2 = pool[state]
                        x1 = [ x for x in codVars[k] if x not in x2]
                        pool = [x1,x2]
                        #### note: remaining in state 1
         
                        
        # Alternator must be initialised "manually" in the special case where there is only one
        # codon option at k = 2.              
        if k == 2 and len(codVars[k]) == 1:
            alternator = alternate()
            state = alternator.next()
            
                
    
    # Convert matrix mat to list of sequences
    mySeqs = []
    for ind in range(seqVars):
        s = ''
        for cod in mat[ind,:]:
            s = s + cod
        mySeqs.append(s)   
             
    return mySeqs

# --------------------------- Main functions' definitions: END ------------------------- #


# ------------------------ Analysis functions' definitions: START ---------------------- #

# Function to evaluate sequence identity - calculate the hamming distance btw two seqs
def hamming_distance(s1,s2):
    """ Takes two strings as input and
    returns the number of mismatches as type integer"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    
    
    
# Calculate pairwise hamming distances from a set of coding sequences   
def hamming_matrix(pool):
    """Takes in a list of sequences (pool) and returns 
    a matrix with pairwise hamming distances"""
    
    # Make sure pool variable is a list
    if type(pool) != list:
        raise ValueError('Pool is not a list')
    
    # Array dimension
    dimLen = len(pool)
    
    # Generate array of zeros
    hammingMatrix = np.zeros((dimLen,dimLen))
    
    # Fill hammingMatrix with hamming distances - symmetrical matrix
    for i,ele_1 in enumerate(pool):
    
        for j,ele_2 in enumerate(pool):
        
            # Only fill in upper triangle
            if j >= i:
                break
            # Calculate distance
            misMatches = hamming_distance(ele_1,ele_2)
            # Fill in array
            hammingMatrix[i,j] = misMatches
            # Same for symmetrical element
            hammingMatrix[j,i] = misMatches
            
    return hammingMatrix
    
def hamming_stats(pool):
    """Generates a hammingMatrix with hamming_matrix()
    and calculates mean distance, minimum and maximum 
    divergence (in %)"""
    
    # Generate array with hamming distances
    divArray = hamming_matrix(pool)
    # Store len of divArray
    dimA = divArray.shape[0]
    
    # Extract upper triangular elements
    iu = np.triu_indices(dimA)
    upTri = [ ele for ele in divArray[iu] if ele != 0 ] # excludes diagonal values
    
    # Calculate mean of upTri and express it in percentage
    meanHam = np.mean(upTri)
    
    # Express mean in percentage - meanHam/len(sequence)
    meanHamPer = float( meanHam/(len(pool[0]))) * 100.0
    
    # Extract the smallest hamming distance
    minHam = np.min(upTri)
    
    # Express min in percentage
    minHamPer = float( minHam/(len(pool[0]))) * 100.0
    
    maxHam = np.max(upTri)
    
    # Express max in percentage - meanHam/len(sequence)
    maxHamPer = float( maxHam/(len(pool[0]))) * 100.0
    
    
    return (meanHamPer,minHamPer,maxHamPer)
    
## Functions to evaluate the longest and average stretches of contiguous homology within a set of sequences

def longest_cont(s1,s2):
    """Works out the longest stretch of identical bases between two
    degenerate coding sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    res = ''
    for ch1, ch2 in zip(s1,s2):
        if ch1 == ch2:
            res += '1'
        else:
            res += '0'        
    #print res   ### Need to subdivide into chunks of contiguous 1's (and 0's)
            
    return np.max([len(x) for x in re.compile("(1+1)*").findall(res)])
    
def average_cont(s1,s2):
    """Works out the average stretch of identical bases between two
    degenerate coding sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    res = ''
    for ch1, ch2 in zip(s1,s2):
        if ch1 == ch2:
            res += '1'
        else:
            res += '0'      
              
    # Extract all homologous segments' lengths
    all_lens =  [len(x) for x in re.compile("(1+1)*").findall(res)]
    # Strip all 0's
    all_homs = [ x for x in all_lens if x!=0 ]   
    return np.mean(all_homs)
    
    
def longest_cont_matrix(pool):

    """Takes in a list of sequences (pool) and returns 
    a matrix with longest stretch of identity"""
    
    # Make sure pool variable is a list
    if type(pool) != list:
        raise ValueError('Pool is not a list')
    
    # Array dimension
    dimLen = len(pool)
    
    # Generate array of zeros
    longestMatrix = np.zeros((dimLen,dimLen))
    
    # Fill longestMatrix with identity values - symmetrical matrix
    for i,ele_1 in enumerate(pool):
    
        for j,ele_2 in enumerate(pool):
        
            # Only fill in upper triangle
            if j >= i:
                break
            # Calculate distance
            idents = longest_cont(ele_1,ele_2)
            # Fill in array
            longestMatrix[i,j] = idents
            # Same for symmetrical element
            longestMatrix[j,i] = idents
            
    return longestMatrix
    
def abs_longest(pool):
    """takes in a list of sequences (pool) and returns the longest stretch of homology
    that can be found in the set"""
    
    # Call longest_cont_matrix() and take the upper triangular region
    mat = np.triu(longest_cont_matrix(pool))
    store_seq = []
    # dump non-zero elements into store_seq
    for row in mat:
        for ele in row:
            if ele != 0:
                store_seq.append(ele)
    # return max
    return np.max(store_seq)    
    
def average_cont_matrix(pool):

    """Takes in a list of sequences (pool) and returns 
    a matrix with average stretches of identity"""
    
    # Make sure pool variable is a list
    if type(pool) != list:
        raise ValueError('Pool is not a list')
    
    # Array dimension
    dimLen = len(pool)
    
    # Generate array of zeros
    average_matrix = np.zeros((dimLen,dimLen))
    
    # Fill longestMatrix with identity values - symmetrical matrix
    for i,ele_1 in enumerate(pool):
    
        for j,ele_2 in enumerate(pool):
        
            # Only fill in upper triangle
            if j >= i:
                break
            # Calculate distance
            idents = average_cont(ele_1,ele_2)
            # Fill in array
            average_matrix[i,j] = idents
            # Same for symmetrical element
            average_matrix[j,i] = idents
            
    return average_matrix

## Pipeline for CAI calculation

# First, define functions that will generate RSCU and relative adaptiveness dictionaries

def RSCU(path):
    """Given a path (argv) to codon usage table of expression host, functions generates a dictionary
    with all RSCU values (return)."""
    
    # Read codon_usage table from data folder
    codon_usage = pd.read_csv(path)
    
    # Initialise dictionary where RSCU values will be stored {triplet:RSCU val}
    RSCU_dict = {}
    
    # Cycle through codon_usage, extract info and load onto dictionary
    for i in range(len(codon_usage['Codon'])):
        RSCU_dict[codon_usage['Codon'][i]] = codon_usage['RSCU'][i]
        
    return  RSCU_dict
    
def relative_adaptiveness(codon_table_full, path):
    """Calculates the relative adaptiveness of each codon given a path to a codon usage table."""
    
    # Define RSCU_dict using RSCU function
    RSCU_dict = RSCU(path)
    
    # For each amino acid, compute the max(RSCU) of synonymous codons
    # Initialise a dictionary, keys are cods, vals max RSCU
    maxRSCU_dict = {}
    for aa in codon_table_full:
        RSCUvals = [ RSCU_dict[cod] for cod in codon_table_full[aa] ]
        for cod in codon_table_full[aa]:
            maxRSCU_dict[cod] = max(RSCUvals)
    
    # Now define a relative adaptiveness dict, keys are cods, vals are ratios RSCU_cod/RSCU_max
    relative_adaptiveness_dict = {}
    for cod in RSCU_dict:
        relative_adaptiveness_dict[cod] = RSCU_dict[cod]/maxRSCU_dict[cod]
        
    return relative_adaptiveness_dict
    

# Calculate CAI for a coding sequence
# seq must be a dna string of len multiple of 3
def CAI_calculator(codon_table_full, path, seq):

    # Generate relative_adaptiveness dictionary
    relative_adaptiveness_dict = relative_adaptiveness(codon_table_full, path)

    if len(seq) % 3 != 0:
        raise ValueError('length of sequence is not a multiple of three')
    
    # Initialise a list where codon adaptiveness values will be stored     
    seq_adaptiveness = []

    # Traverse coding sequence and extract codons
    for i in range(0,len(seq),3):
        cod = seq[i:i+3]
        cod_adaptiveness = relative_adaptiveness_dict[cod]
        seq_adaptiveness.append(cod_adaptiveness)
    
    # Calculate the geometric mean for the list seq_adaptiveness (= CAI)
    CAI = scipy.stats.mstats.gmean(seq_adaptiveness)
    
    return CAI

# ------------------------- Analysis functions' definitions: END ----------------------- #

################################### RUN THE PROGRAM ######################################

# Define input

# Protein and number of degenerate sequences to generate
my_prot = raw_input('Paste protein sequence (press Return) then number of degenerate sequences to generate (press Return): ')  
# Make sure sequence is upper case, there are no blanks or newline characters
my_prot = my_prot.upper().replace(' ','').replace('\n','')
# Number of degenerate sequences to generate
seqVars = int(raw_input(''))
# Expression host (codon_table is generated accordingly)
path = input_expression_host()  # Note: path will also be used with RSCU
codon_table = generate_codon_table(path = path, codon_table_full = codon_table_full)

# Launch lib_generator() to generate sequence library
mySeqs = lib_generator(my_prot = my_prot, seqVars = seqVars, codon_table = codon_table)
# 
# print mySeqs
# 
# Process data 
# 1) Generate ID list
IDlist = [ 'seq' + str(i) for i in range(len(mySeqs)) ]
# 
# Generate CAI values
CAIlist = []
for seq in mySeqs:
    CAIlist.append(CAI_calculator(codon_table_full = codon_table_full, path = path, seq = seq))
#     
## Generate GC values
GClist = []
for seq in mySeqs:
    GClist.append(GC(seq))
#     
## Create dataframe with output sequences and statistics
output_data = {'seq ID': IDlist, 'Sequence': mySeqs, 'CAI': CAIlist, 'GC %': GClist}
df = pd.DataFrame(output_data)
print df 
# 
## Communicate hamming distance stats and longest stretch of homology
print 'Stats: Mean, minimum and maximum hamming distances in the sequence set are (per cent):', hamming_stats(mySeqs)
print 'Stats: Longest stretch of homology between any two sequences (in bp):', abs_longest(mySeqs)
# 
df.to_csv('codingSequenceVariants.csv', sep = ',')




































