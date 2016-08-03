#!/usr/bin/env python

"""kmer_index.py: A k-mer index for indexing a text."""

__author__ = "Shubham Gupta"

import bisect
class index(object):
    def __init__(self, Text, k):
        # Define class variable
        self.k = k
        self.index = []
        for i in range(len(Text)-k+1):
            self.index.append((Text[i:i+k],i)) #tuple in python
        self.index.sort()
        
    def query(self, Pattern):
        kmer = Pattern[:self.k]
        i = bisect.bisect_left(self.index,(kmer,-1))
        # bisect_left is usedd to find insertion position
        # but we are using it here to find out kmer instances
        # -1 is used so that we get 1st instance
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1 #Increment till we hit break or
            #in other words finishes our matching kmer
        return hits

def queryIndex(Text, Pattern, OBJECT):
    k = OBJECT.k
    offsets = []
    for i in OBJECT.query(Pattern):
        if Pattern[k:] == Text[i+k: i+len(Pattern)]:
            offsets.append(i) #Verification step
    return offsets

Text = 'ATGCTGATGCTAGCCAGATGCTGA'
Pattern = 'AGC'
Obj = index(Text,2)
queryIndex(Text,Pattern,Obj)

def ApproxQueryIndex(Text, Pattern, OBJECT, n):
    k = OBJECT.k
    Seg_len = len(Pattern)/(n+1)
    hitsNum = 0
    All_Match = set()
    if Seg_len < k:
        print "Segment length is less than kmer length"
        return
    offsets = []
    for i in range(n+1):
        offsets = []
        start = i*Seg_len
        end = min((i+1)*Seg_len, len(Pattern))
        SubPattern = Pattern[start: end]
        hits = OBJECT.query(SubPattern)
        hitsNum += len(hits)
        for h in hits:
            if h-start < 0 or h-start+ len(Pattern) > len(Text):
                continue
            
            mismatches = 0
            for j in range(start):
                if not Pattern[j] == Text[h-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            
            for j in range(end, len(Pattern)):
                if not Pattern[j] == Text[h-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            
            if mismatches <= n:
                All_Match.add(h-start)
    return list(All_Match), hitsNum

Text = readGenome('chr1.GRCh38.excerpt.fasta')
Pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
kmerLength = 8
Obj = index(Text,kmerLength)
# n should always be greater than 1
occurence, hitsNum = ApproxQueryIndex(Text,'GGCGCGGTGGCTCACGCCTGTAAT',Obj,2)
print len(occurence)
print hitsNum

#occurence = naive_2mm('GGCGCGGTGGCTCACGCCTGTAAT',Text)
#print len(occurence)

import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def query_subseq(Pattern, Text, OBJECT):
    k = OBJECT.k
    offsets = []
    numHits = 0
    numHits = len(OBJECT.query(Pattern))
    numHits += len(OBJECT.query(Pattern[1:]))
    numHits += len(OBJECT.query(Pattern[2:]))
    # This function is not complete just to get output
    return offsets, numHits

t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = query_subseq(p, t, subseq_ind)
print num_index_hits
