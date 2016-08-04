# This document has some very basic code require to have genome assembly logic

__author__ = "Shubham Gupta"

def overlap(a, b, minL = 3): # Suffix of 'a' and prefix of 'b'
    start = 0
    occurence = []
    while True:
        start = a.find(b[:minL], start)# Gives index where prefix of b is in a
        if start == -1:
            return 0
        if b.startswith(a[start:]):# This is to check that letter after minL are also present in b
            # if start+ minL is not the last letter
            return len(a)-start
        start += 1

overlap('TTACGT', 'ACGTGTAC')

from itertools import permutations
list(permutations([1,2,3,4],3))

def Basic_Overlap_match (reads, minL):
    olaps = {}
    for (a,b) in permutations(reads,2):
        olen = overlap(a,b,minL)
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps
Reads = ['TTACGT', 'ACGTGTAC', 'GTACTACG']
print(Basic_Overlap_match(Reads,3))


from itertools import permutations
def readFastq(filename):
    Qval = []
    genome = []
    with open(filename, 'r') as fileHandle:
        while True:
            fileHandle.readline()
            seq = fileHandle.readline().rstrip()
            fileHandle.readline()
            q = fileHandle.readline().rstrip()
            if seq == '':
                break
            genome.append(seq)
            Qval.append(q)
    return genome, Qval


def createKmerSet(genomeReads, minL):
    dictOfKmer = {}
    
    # Creating dictionary of all kmers and reads they are present in 
    for read in genomeReads:
        for i in range(len(read) - minL +1):
            kmer = read[i:i+ minL]
            if kmer not in dictOfKmer.keys():
                dictOfKmer[kmer] = set([read])
            dictOfKmer[kmer].add(read)
    
    olaps = {}
    for keys in dictOfKmer:
        keyval = dictOfKmer[keys]
        for (a,b) in permutations(keyval,2):
            olen = overlap(a,b,minL)
            if olen > 0 and a!=b:
                olaps[(a,b)] = olen
    return olaps.keys()
    
import datetime as d
t0 = d.datetime.now()
genomeReads,_ = readFastq('ERR266411_1.first1000.fastq')
reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
overlap_all_pairs = createKmerSet(genomeReads, 30)

t1 = d.datetime.now()
print (t1 - t0).total_seconds()
print len(overlap_all_pairs)