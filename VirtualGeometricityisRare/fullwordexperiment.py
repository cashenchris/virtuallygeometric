#!/usr/bin/env python

## usage:  ./fullwordexperiment.py rank (minlength maxlength samplesize)

import sys
import matplotlib.pyplot as plt
import freegroup
import whiteheadgraph.split.split as split
from math import exp, log

######
rank = int(sys.argv[1])
try:
    minlength = int(sys.argv[2])
except:
    minlength = 1
try:
    maxlength = int(sys.argv[3])
except:
    maxlength = 40
try:
    sample_size = int(sys.argv[4])
except:
    sample_size = 100


outputname = 'rank'+str(rank)+'fullwords.txt'

####

outfile = open(outputname,'a')
F=freegroup.FGFreeGroup(numgens=rank)
length=minlength
print 'length','full','notfull'

while True:
    if length>maxlength:
        break
    full=0
    notfull=0
    for i in range(sample_size):
        w=F.random_word(length)
        if split.contains_all_3_letter_subwords(w): # w, as a cyclic word, contains all 3 letter subwords
            full+=1
        else:
            notfull+=1
        print length, full, notfull,'\r',
        sys.stdout.flush()
    print
    outfile.write(str(length)+ '   '+ str(full)+ '   '+ str(notfull)+ '   ' + '\n')
    outfile.close()
    outfile = open(outputname,'a')
    length+=1



