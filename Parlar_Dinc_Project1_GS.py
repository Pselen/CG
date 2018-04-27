#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 18:04:01 2018

@author: selen parlar 150113049 & bilal dinc 150113008
"""
# count ve pseudo count hatalari, sumcountlar farkli dogru degil.
from random import randint,sample
import random
import time
import sys

random.seed(5)
# Declare variables
global k,t,implant,scores,n
implanted_dna = []
dna = []

#file = sys.argv[1]
#k = int(sys.argv[2])
nucleotides = ['A','C','G','T']

t = 20 # number of sequences in dna.
n = 100
k = 10
file = '/home/selen/Desktop/CG/test_data.txt'
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
#%%
# Returns a list of random kmers of each sequence in dna which represent motifs.
def random_motifs(dna):
    randoms = []
    for sequence in dna:
        start = randint(0, len(sequence)-k)
        randoms.append(sequence[start:start+k])
    return randoms

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
def Profile_with_pseudocounts(motifs):
    count = Count_with_Pseudocounts(motifs)
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

#%% Gibbs sampler algorithm.
def gibbs_sampler(dna,n):
    t = len(dna)
    motifs = random_motifs(dna)
    best_motifs = motifs
    best_score = Score(best_motifs)
    for step in range(n):
        i = randint(0,t-1)
        profile = (Profile_with_pseudocounts(motifs[:i] + motifs[(i+1):]))
        probability_profile = [probability_of_kmer(dna[i][s:s+k],profile) for s in range(len(dna[i])-k+1)]
        p = random_probability(probability_profile)
        motifs[i]=dna[i][p:p+k]
        score=Score(motifs)
        if score < best_score:
            best_motifs = motifs[:]
            best_score = score
            print(best_score, end = ' ')
    return (best_motifs,best_score)
#%% Repeat GS for n times.
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
print('Gibbs Sampler')
s = time.time()
best_motifs =repeated_gibbs_sampler(dna,n)
e = time.time()
print('\nExecution time: '+ str(e-s))
# Generate concensus string within best motifs.
consensus = Consensus(best_motifs)
print('\nMotifs:')
for motif in best_motifs:
    print(motif)
print('Consensus string: ' + consensus)
