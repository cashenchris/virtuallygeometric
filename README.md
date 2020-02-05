This code is for testing whether a multiword in a free group is virtually geometric. Towards this end it implements various algorithms related to elements and subgroups of finitely generated free groups, including computing the JSJ decomposition of a free group relative to a multiword. 

Requires python 2.7 and the packages 'pexpect' and 'networkx<=1.9.1'. API changes in networkx1.10 break things.  

All of the geometricity and virtual geometricity code requires the program 'heegaard', which is available at [http://www.math.uic.edu/t3m/](http://www.math.uic.edu/t3m/). You must also edit the file ./geometric/heegaard.py and fill in the path to the heegaard binary.


Examples of testing virtual geometricity from the command line:

```
#!bash
$ ./vg.py a aabAAAB
True
$ ./vg.py aabbccacb
False
```

The folder VirtualGeometricityisRare contains scripts for recreating the raw data and graphs from the paper:
Virtual geometricity is rare, Cashen and Manning,
LMS Journal of Computation and Mathematics, 18 (2015), no. 1, 444–455.
https://doi.org/10.1112/S1461157015000108

This repository was migrated from the one referenced in the paper (https://bitbucket.org/christopher_cashen/virtuallygeometric).

Questions/Comments to Chris Cashen or Jason Manning
