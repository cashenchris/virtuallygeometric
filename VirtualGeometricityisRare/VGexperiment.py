#!/usr/bin/env python

## usage:  ./VGexperiment.py rank (minlength maxlength samplesize)

import sys
import freegroup
import virtuallygeometric
from time import localtime,strftime
import whiteheadgraph.split.split as split




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



outputname = 'rank'+str(rank)+'_threeway.txt'
badwordname = 'rank'+str(rank)+'_threeway_bad.txt'

######
outfile = open(outputname,'a')
length=minlength
F = freegroup.FGFreeGroup(numgens=rank)
print 'length','geometric','vg','not_vg'
outfile.write(strftime("%Y-%m-%d,%H:%M",localtime())+'\n')
while True:
    if length>maxlength:
        break
    geometric = 0
    vg = 0
    not_vg = 0
    for i in range(sample_size):
        w = F.random_word(length)
        if not split.missing_3_letter_subwords(w):
            not_vg+=1
            print length,geometric,vg,not_vg,'\r',
            sys.stdout.flush()
            continue
        try:
            heegaard_yes = virtuallygeometric.is_orientably_geometric([w.alpha()])
        except RuntimeError:
            heegaard_yes = False
        if heegaard_yes:
            geometric+=1
            print length,geometric,vg,not_vg,'\r',
            sys.stdout.flush()
            continue
        try:
            if virtuallygeoemtric.is_virtually_geometric(F,[w]):
                vg+=1
            else:
                not_vg+=1
        except StandardError as foo:
            # try to keep a record of the bad words
            print foo, w  
            wordsfile = open(badwordname,'a') 
            wordsfile.write(str(type(foo))+foo.message+'  '+w.alpha()+'\n')
            wordsfile.close()
        print length,geometric,vg,not_vg,'\r',
        sys.stdout.flush()
    outfile.write(str(length)+ '   '+ str(geometric)+ '   '+ str(vg)+ '   ' + str(not_vg)+ '\n')
    outfile.close()
    outfile = open(outputname,'a')
    print length,geometric,vg,not_vg,'\r',
    sys.stdout.flush()
    print
    length+=1
outfile.close()
