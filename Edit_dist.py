# This file contains code to calculate edit distance between two strings

__author__ = "Shubham Gupta"

def edistRecursive(a,b):
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)
    
    if a[-1] == b[-1]:
        delt = 0
    else:
        delt = 1
    return min(edistRecursive(a[:-1],b[:-1]) + delt,
               edistRecursive(a,b[:-1]) + 1,
               edistRecursive(a[:-1],b) + 1)


def edistMatRecursive(a,b):
    mat = []
    for i in range(len(a)+1):
        mat.append([0]*(len(b)+1))
    
    for i in range(len(b)+1):
        mat[0][i] = i
        
    for i in range(len(a)+1):
        mat[i][0] = i
        
    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            distHor = mat[i][j-1] + 1
            distVer = mat[i-1][j] + 1
            if a[i-1] == b[j-1]:
                distComm = mat[i-1][j-1]
            else:
                distComm = mat[i-1][j-1] + 1
            mat[i][j] = min (distHor,distVer,distComm)
    
    return mat[-1][-1]

Strng1 = 'Shakespeare'
Strng2 = 'shakes peare'

import datetime as d
t0 = d.datetime.now()
dist = edistRecursive(Strng1,Strng2)
t1 = d.datetime.now()
print (t1 - t0).total_seconds()
print dist

#76.733 sec
#2


t0 = d.datetime.now()
dist = edistMatRecursive(Strng1,Strng2)
t1 = d.datetime.now()
print (t1 - t0).total_seconds()
print dist

#0.001 sec
#2
alphabet = ['A', 'C', 'G', 'T']
penalty = [[0, 4, 2, 4, 8], \
           [4, 0, 4, 2, 8], \
           [2, 4, 0, 4, 8], \
           [4, 2, 4, 0, 8], \
           [8, 8, 8, 8, 8]]

def globalAlignment(Pattern, Text):
    mat = []
    for i in range(len(Pattern)+1):
        mat.append([0]*(len(Text)+1))
    
    # Removed the code that initializes first row with Position index
    
    # Filling the first column with the penalty values of skipping that character
    for i in range(1, len(Pattern)+1):
        mat[i][0] = mat[i-1][0] + penalty[alphabet.index(Pattern[i-1])][-1]
        
    for i in range(1, len(Pattern)+1):
        for j in range(1, len(Text)+1):
            distHor = mat[i][j-1] + penalty[-1][alphabet.index(Pattern[j-1])]
            distVer = mat[i-1][j] + penalty[alphabet.index(Pattern[i-1])][-1]
            if Pattern[i-1] == Text[j-1]:
                distComm = mat[i-1][j-1]
            else:
                distComm = mat[i-1][j-1] + penalty[alphabet.index(Pattern[i-1])][alphabet.index(Text[j-1])]
            mat[i][j] = min (distHor,distVer,distComm)
    
    return mat[-1][-1]

a = 'TACCAGATTTCGA'
b = 'TACCAGATTTCAA'
globalAlignment(a, b)
