# This file contains code to calculate shortest common superstring (SCS)
# Brute force, Greedy Assembly and de Bruijn graph

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

import itertools
def shortest_ss(SetOfString):
    scs = None
    for ss in itertools.permutations(SetOfString):
        sup = ss[0]
        for i in range(1, len(ss)):
            olen = overlap(sup,ss[i], 1)
            sup += ss[i][olen:]
        if scs is None or len(sup)<len(scs):
            scs = sup
    return scs

import itertools
SetOfString = ['ACGTCAGTC', 'ACGTAGCGTAC', 'CAGTGCATGC']
shortest_ss(SetOfString)

def shortest_ss_list(SetOfString):
    scs = []
    for ss in itertools.permutations(SetOfString):
        sup = ss[0]
        for i in range(1, len(ss)):
            olen = overlap(sup,ss[i], 1)
            sup += ss[i][olen:]
        
        if not scs:
            scs.append(sup)
           
        for element in scs:
            if len(sup) < len(element):
                scs = []
                scs.append(sup)
            elif len(sup) == len(element) and sup != element:
                scs.append(sup)
    return scs
SetOfString = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
shortest_ss_list(SetOfString)


def Maximum_Overlap(SetOfString, minL):
    readA = ""
    readB = ""
    MaxOlen = 0
    for (a,b) in itertools.permutations(SetOfString, 2):
        olen = overlap(a, b, minL)
        if olen > MaxOlen:
            readA, readB = a,b
            MaxOlen = olen
    return readA, readB, MaxOlen


def greedyAssembly(SetOfString, minL):
    a, b, olen = Maximum_Overlap(SetOfString, minL)
    while olen>0:
        SetOfString.remove(a)
        SetOfString.remove(b)
        SetOfString.append(a+b[olen:])
        a, b, olen = Maximum_Overlap(SetOfString, minL)
    return "".join(SetOfString)

SetOfString = ['ACGTCAGTC', 'ACGTAGCGTAC', 'CAGTGCATGC']
print greedyAssembly(SetOfString, 3)
print shortest_ss(SetOfString)

def de_bruijn_ize(st, k):
    """ Return a list holding, for each k-mer, its left
        k-1-mer and its right k-1-mer in a pair """
    edges = [] # We want multiple edges between same nodes thats why it is taken as list
    nodes = set() # We dont want multiple nodes of same name
    for i in range(len(st) - k + 1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges

def visualize_de_bruijn(st, k):
    """ Visualize a directed multigraph using graphviz """
    nodes, edges = de_bruijn_ize(st, k)
    dot_str = 'digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str += '  %s [label="%s"] ;\n' % (node, node)
    for src, dst in edges:
        dot_str += '  %s -> %s ;\n' % (src, dst)
    return dot_str + '}\n'

nodes, edges = de_bruijn_ize("ACGCGTCG", 3)

%install_ext https://raw.github.com/cjdrake/ipython-magic/master/gvmagic.py

%load_ext gvmagic
%dotstr visualize_de_bruijn("ACGCGTCG", 3)
