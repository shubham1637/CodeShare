
__author__ = "Shubham Gupta"

# This document contains some basic online algorithms for pattern matching
import random
seq = ''
for _ in range(10):
	seq += random.choice('ACGT')
print(seq)


def reverseComplement(DNA):
    Complement = {'A':'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    RS = ''
    for base in DNA:
        RS = Complement[base] + RS
    return RS
reverseComplement('ATCGTAGCTGA')

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as fileHandle:
        for line in fileHandle:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
gene = readGenome('lambda_virus.fa')
gene[:10]

import collections
collections.Counter(gene)

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
genome, Qstrng = readFastq('SRR835775_1.first1000.fastq')

genome[:4]
Qstrng[:4]
print(genome[:4])

def phred33ToQval(char):
    return ord(char)-33
phred33ToQval('#')

def createHist(List):
    Hist = [0]*50
    for STRNG in List:
        for char in STRNG:
            q = phred33ToQval(char)
            Hist[q] +=1
    return Hist
Hist = createHist(Qstrng)
print Hist

import matplotlib.pyplot as plot
plot.bar(range(len(Hist)),Hist)

import random
def generateReads(genome, numReads, readLen):
    ''' Generate reads from random positions in the given genome. '''
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome)-readLen) - 1
        reads.append(genome[start : start+readLen])
    return reads

genome = readGenome('phix.fa')
def naiveSearch(Pattern, genome):
    outIndices = []
    Match = False
    for i in range(len(genome)- len(Pattern) +1):
        for j in range(len(Pattern)):
            if not Pattern[j] == genome[i+j]:
                Match = False
                break
            Match = True
        if Match == True:
            outIndices.append(i)
    return outIndices

genomeReads,_ = readFastq('ERR266411_1.first1000.fastq')

numMatched = 0
occurence = []
Total_reads = 0
for r in genomeReads:
    occurence = naiveSearch(r,genome)
    numMatched += len(occurence)
print "%d / %d matched" %(numMatched, len(genomeReads)) 

# This is quite low. One reason could be sequencing error and other could be the differences between individuals means Test is virus and reference is bacterial genome.

numMatched = 0
occurence = []
Total_reads = 0
for r in genomeReads:
    r = r[:30]
    occurence = naiveSearch(r,genome)
    numMatched += len(occurence)
print "%d / %d matched" %(numMatched, len(genomeReads)) 

# This is good but not close to 100. Keep in mind that gene is double stranded and we could be reading only one strand.

numMatched = 0
occurence = []
Total_reads = 0
for r in genomeReads:
    r = r[:30]
    occurence = naiveSearch(r,genome)
    occurence.extend(naiveSearch(reverseComplement(r),genome))
    numMatched += len(occurence)
print "%d / %d matched" %(numMatched, len(genomeReads))   


def naiveStrandAware(Pattern, genome):
    outIndices = []
    reversePattern = reverseComplement(Pattern)
    OnlyOnce = (Pattern == reversePattern)
    outIndices1 = []
    Match = False
    for i in range(len(genome)- len(Pattern) +1):
        for j in range(len(Pattern)):
            if not Pattern[j] == genome[i+j]:
                Match = False
                break
            Match = True
        if Match == True:
            outIndices.append(i)
        if not OnlyOnce:
            for j in range(len(reversePattern)):
                if not reversePattern[j] == genome[i+j]:
                    Match = False
                    break
                Match = True
            if Match == True:
                outIndices1.append(i)        
    return outIndices, outIndices1


def naive_2mm(Pattern, genome):
    outIndices = []
    Match = False
    count = 0
    for i in range(len(genome)- len(Pattern) +1):
        count = 0
        for j in range(len(Pattern)):
            if not Pattern[j] == genome[i+j]:
                count += 1
                if count > 2:
                    Match = False
                    break
            Match = True
        if Match == True:
            outIndices.append(i)
    return outIndices

gene = readGenome('lambda_virus.fa')
occurence, revOccurence = naiveStrandAware('AGT',gene)
print "1st oc %d" %occurence[0]
print "2nf de %d" %revOccurence[0]

gene = readGenome('lambda_virus.fa')
occurence = naive_2mm('AGGAGGTT',gene)
print "1st oc %d" %len(occurence)


reads, Qval = readFastq('ERR037900_1.first1000.fastq')
PsnSeqQ = [0]*len(Qval[0])
for i in range(len(Qval)):
    j = 0
    for char in Qval[i]:
        q = phred33ToQval(char)
        PsnSeqQ[j] += q
        j += 1
print PsnSeqQ
PsnSeqQ.index(min(PsnSeqQ))
