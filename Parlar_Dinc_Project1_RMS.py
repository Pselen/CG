#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 22:17:58 2018

@author: selen parlar 150113049 & bilal dinc 150113008
"""
from random import randint,sample
import random
import time
import sys
random.seed(5)

# Declare variables
global k,t,implant,scores, n
nucleotides = ['A','C','G','T']
t = 20 # number of sequences in dna.
n = 100 # number of repeats
dna = []
implanted_dna = []
#file = sys.argv[1]
#k = int(sys.argv[2])

file = '/home/selen/Desktop/CG/test_data.txt'
k=10

#%%
# Read file line by line into dna.
with open(file, 'r') as f:
    for line in range(t):
        dna.append(f.readline().strip())
#%%
# Generate 10-mer randomly to implement. 
def generate_implant():
    implant = []
    for i in range(k):
        s = randint(0,3)
        implant.append(nucleotides[s])
    return implant
#%%
# Mutate implant randomly.
def generate_mutated_implant(implant):
    s = sample(range(0,9),4)
    mutation = {'A': 'T', 'C': 'G', 'G':'C', 'T': 'A'}
    for i in s:
        implant[i] = mutation.get(implant[i])
    return implant
#%%
# Implant the given kmer into dna.   
def implant(dna):     
    #implant = generate_mutated_implant()
    implant = generate_implant()
    for sequence in dna:
           mutated_implant = generate_mutated_implant(implant)
           print('Mutated 10-mer: '+ str(implant))
           start = randint(0,len(sequence)-k)
           sequence = sequence[:start] + ''.join(mutated_implant) + sequence[(start+len(mutated_implant)):]
           implanted_dna.append(sequence)
    
    return implanted_dna
#%%
# Returns a list of random kmers of each sequence in dna which represent motifs.
def random_motifs(dna):
    randoms = []
    for sequence in dna:
        start = randint(0, len(sequence)-k)
        randoms.append(sequence[start:start+k])
    return randoms    
#%%
# Calculate the occurrances of each nucleotide in a given motif. Returns a dict obj.
def Count(motifs):
    count = {}
    # Set values to 0.
    for n in "ACGT":
        count[n] = []
        for i in range (k):
            count[n].append(0) 
    for i in range(t):
        for j in range(k):
            n = motifs[i][j]
            count[n][j] += 1
    return count
#%%
# Calculate score of motifs.
def Score(motifs):
    count = Count(motifs)
    max_count=[]
    for j in range(k):
        temp = []
        for i in range(4):
            temp.append(count[nucleotides[i]][j])
        max_count.append(max(temp))
    score=0
    for i in range(k):
        score += t - max_count[i]
    return score
#%%
# Return profile matrix as a dictionary.
def Profile(motifs):
    count = Count(motifs)
    profile = count.copy()
    sum_count = 0
    for n in nucleotides:
        sum_count += count[n][0]
    for j in range(k):
        for i in range(4):
            profile[nucleotides[i]][j] = count[nucleotides[i]][j]/sum_count
    return profile
#%%
# Calculate probability of a kmer.
def probability_of_kmer(kmer, profile):
    probability = 1.0
    for i in range(len(kmer)):
        probability *= profile[kmer[i]][i]
    return probability

#%%
# Calculate all kmer's probabilities and find best kmer with the highest probability.
def most_probable_kmer(dna, profile):
    probable_kmers = []
    for sequence in dna:
        maxvalue = 0
        probable_kmer = ''
        for i in range(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            value = probability_of_kmer(kmer,profile)
            if  value > maxvalue:
                maxvalue = value
                probable_kmer = kmer
        probable_kmers.append(probable_kmer)
    return probable_kmers
        
#%% 
# Ramdomized motif search algorithm.
def randomized_motif_search(dna):
    best_motifs = random_motifs(dna)
    best_score = Score(best_motifs)
    while True:
        profile = Profile(best_motifs)
        motifs = most_probable_kmer(dna, profile)
        score = Score(motifs)
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            return best_motifs        
            
#%% 
# Repeat RMS for n times.
def repeated_randomized_motif_search(dna,n):
    best_score = float('inf')
    scores = []
    best_motifs = []
    i = 0
    print ('Scores in '+str(n)+' episodes:')
    while True:
        motifs = randomized_motif_search(dna)
        score = Score(motifs)
        if score < best_score:
            best_motifs = motifs
            best_score = score
            scores.append(best_score)
            print(best_score, end = ' ')
            i=0
        else:
            i +=1
        if i>n:
            break
    #print(scores)
    return best_motifs
#%%
# Generate consensus string depending on counts in motifs.
def Consensus(motifs):
    count = Count(motifs)
    consensus = ""
    for j in range(k):
        m=0
        frequent_symbol = ""
        for n in nucleotides:
            	if count[n][j] > m:
                    m = count[n][j]
                    frequent_symbol = n
        consensus += frequent_symbol
    return consensus
#%%
# Main 

# Implant mutated kmer into dna read from file.
dna = implant(dna)
# Implanted kmer ['C', 'G', 'A', 'A', 'G', 'A', 'G', 'A', 'G', 'T']
#dna = ['TTCTCGCCTACCTGCATGCCGAACTGCGGGGACCTCTCAAATCTTGAATACCTATTCAGCAATCGTCGGTTATCGAAGAGAGTGTCTGCCCCATTTCTAA', 'ATGCCGAACTGCGGGGACCTCTCAAATCTTGCGAAGAGAGTCAGCAATCGTCGGTTATCAGCCGCTATGTCTGCCCCATTTCTAGCAATCCTTTTCATGA', 'TCGAAGAGAGTGGTGTCGCGCTGACGCTGGTAGATGAGTAAATAGATATGCACCCTGAGGCGCGTATCACCTCTCTCAGTACACAAAAGCGTAATGAGTA', 'AAGGGAACCGTGGCAATTGAGGCACTACGAAGAGAGTGCAACATAGTGTCAATGCTCGGCGGTAGCATTCTGCAGGGGCAAGCGGTAACGGAACAAGAAA', 'TCCTCATCTTGGGTCCGCCGACATGCTATTGATGGCGGAGCGGTGTGGAACACGAAGAGAGTTAACCTGCCCTCCCAGTGCGAGCTACCTGTCTAGTTAA', 'TAGTTGCTGTAAATCAGAACCATTTTTATTGGGGGCGAAGAGAGTGCCGAGCAACTGCTGAGGGCAAGAACATAGCTACCATGACCAACATCAATAGCAA', 'AACGTTGTCCCGAAAAAGGCTGACGAAGAGAGTATGCAAGCCCCTACCTTCATACCGACATTCGACAACTTCTATAGCAAGTCCGAGGTAAACTACTAGA', 'TGCACCGGTGCACTGCCTCGTCAAAACCACTGACAGTCACGCTAAGCATCGAAGAGAGTTCGCTGTTGACGGACCCATTCGGAAGACAGGTTCACACTTA', 'TTCTACAATAAAAGTTGACCCGAAGAGAGTGGCATGCCTCGGACGCAGCGCATGCGCTAAGGAGATGAGTTGGTCTGGCGTTGCAAAGAGTCCCGCACAA', 'GTACATTTTCGAAGAGAGTTTCTTTCCAAGTTACGGGCTGTTGGAATAGCAACACACTCGACGCGCGTTGTCAAAACTGCGGAGTCACGTCAGGATGGAA', 'GCCTTGAACGAGCCACGCGAAGAGAGTCGAATTGAGCGGTCTGCCCAACTGGAATTCTAAGCACGTTTGCCTGAGAACGCGACCTAGATCGGCTAATCTA', 'CGATGCTGGTGGTCGTCCCACACCCGATTATACCAACTGGATAGGCAAACCATTCCGAGCAACTGAGAGCACGCGTGTTCGAAGAGAGTACATGTCAAAA', 'TGAAAGTTACGAGTGTAGCCTCACACGGTTAGTATAGATCGCACCGAGCTACGAAACAAGACCTTAGGTGTCCGGACTGCGAAGAGAGTTCTCCCAACAA', 'AGAGCGTAAATACTCTGGTGAATGGGGTGGCTGAGGGCGATTATTCGCAATCCGCGCGAAGAGAGTAAGTCATAGTCAACTGCTTTAGCAAGAGTATGGA', 'GCGACAAGCTGAGGCGCGAAGAGAGTATTTTCTCGGCGCTATTCTTTTTTTTCCGGTGCTTCGATCACTGGAATAGCATGAAAACCCAGGCTCCGTGAAA', 'ATCTGAGTGTCAACTGCGAAGAGAGTCAGTCCAGTTTCAAACCGCAAGGTACCCGTTGCCACAATATCGAGATCAGTGCTAGTTACATGACTCAACTGTA', 'CGAAGAGAGTCCCGCGTGAAACTAACTGGAATAGCGGTCTCTTATCAGGCCCCTCCGATTGCTGATCCAACCGTTTCTACTCGTCGAACTGTCACATACA', 'CGAAGAGAGTCAAACAACTAGACAGTATTATGTCTTGATATGACATCGAGCCATCAGACCCGCCATGGGGGCAGCCCTTTCTGGATGTATAGCTGCGCAA', 'AGTGGTTCTTGCGCTATGGTACGGGTCGAAGAGAGTTTTGAGACTCACAAATACCGCCCCAGGGTCTAAGATGCACCCGTCGGCAAGCAACTGGAATATA', 'TCTCCTTTTCAACACCCATCCTTATTACGAAGAGAGTAGTGGCGAGGCGAGTCTGTGTCAAGTCCACAGCATAAGGGCAAAACGAAGGGAGACCAATAAA']
print('Randomized Motif Search')
s = time.time()
best_motifs = repeated_randomized_motif_search(dna,n)
e = time.time()
print('\nExecution time: '+ str(e-s))
# Generate concensus string within best motifs.
consensus = Consensus(best_motifs)
print('\nMotifs:')
for motif in best_motifs:
    print(motif)
print('Consensus string: ' + consensus)