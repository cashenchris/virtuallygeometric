#!/usr/bin/env python2

## usage:  ./vg.py word1, word2, ...
## decide if wordlist is virtually geometric
## wordi is a string where a-z represent generators of free group and A-Z their inverses

## examples:
## $ ./vg.py a aabAAAB
## True
## $ ./vg.py aabbccacb
## False

import freegroup
import virtuallygeometric
import sys

####

F=freegroup.FGFreeGroup(numgens=26)
if len(sys.argv)>1:
    wordlist=[F.word(str(x)) for x in sys.argv[1:]]
    print virtuallygeometric.is_virtually_geometric(F,wordlist)
else:
    print "Usage: ./vg.py word1, word2, ..."
