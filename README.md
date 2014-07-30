This code is for testing whether a multiword in a free group is virtually geometric. Towards this end it implements various algorithms related to elements and subgroups of finitely generated free groups, including computing the JSJ decomposition of a free group relative to a multiword. 

All of the geometricity and virtual geometricity code requires the program 'heegaard', which is available at [http://www.math.uic.edu/t3m/](http://www.math.uic.edu/t3m/). You must also edit the file ./geometric/heegaard.py and fill in the path to the heegaard binary.

The 'pexpect' and 'networkx' packages are also required.

Examples of testing virtual geometricity from the command line:

```
#!bash
$ ./vg.py a aabAAAB
True
$ ./vg.py aabbccacb
False
```

Questions/Comments to Chris Cashen or Jason Manning