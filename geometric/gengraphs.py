import itertools
from stallings import xgraph

if 'combinations' not in dir(itertools):
    def combinations(iterable, r):
        # combinations('ABCD', 2) --> AB AC AD BC BD CD
        # combinations(range(4), 3) --> 012 013 023 123
        # Copied this code from python 2.6 documentation.
        pool = tuple(iterable)
        n = len(pool)
        indices = range(r)
        yield tuple(pool[i] for i in indices)
        while 1:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i+1, r):
                indices[j] = indices[j-1] + 1
            yield tuple(pool[i] for i in indices)


    def permutations(iterable, r=None):
        # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
        # permutations(range(3)) --> 012 021 102 120 201 210
        # Copied this code from python 2.6 documentation.
        pool = tuple(iterable)
        n = len(pool)
        r = n if r is None else r
        indices = range(n)
        cycles = range(n, n-r, -1)
        yield tuple(pool[i] for i in indices[:r])
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i+1:] + indices[i:i+1]
                    cycles[i] = n - i
                else:
                    j = cycles[i]
                    indices[i], indices[-j] = indices[-j], indices[i]
                    yield tuple(pool[i] for i in indices[:r])
                    break
            else:
                return

    def product(*args, **kwds):
        # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)
else:
    combinations = itertools.combinations
    permutations = itertools.permutations
    product = itertools.product

def negcyc(elts):
    pool = tuple(elts)
    r = len(pool)
    ind = [r-1]+range(r-1)
    return tuple(pool[i] for i in ind)

    
def deccparts(elements,max_size=None):
    """
    Partitions of a list, with largest parts first.
    """
    n = len(elements)
    if n == 0:
	return
    max_size = n if max_size==None else min(n,max_size)
    for part1 in range(1,max_size+1):
	if part1==n:
	    yield [elements]
	else:
	    for smallerpart in deccparts(elements[part1:],part1):
		yield [elements[:part1]] + smallerpart

def perm_cc(elements):
    """
    Iterator giving permutations of elements, but only one per conjugacy class.
    """
    pool = list(elements)
    n = len(pool)
    for partlist in deccparts(elements):
	nextperm = ()
	for subs in partlist:
	    nextperm += negcyc(subs)
	yield nextperm
    

def graphs_old(r,k):
    """
    This should be an iterator, which spits out xgraphs which k-fold
    cover the wedge of r circles.
    """
    vertices = range(k)
    # To make the graph, we need r (not necessarily distinct)
    # permutations of k elements, permuted somehow.  This could be
    # more efficient; we only really need r-tuples of permutations *up
    # to conjugacy in S_k*, not arbitrary r-tuples.
    for perms in combinations(product(permutations(vertices),repeat=r),1):
        for match in permutations(perms):
            # match[0][i][j] tells us the target of the edge labeled
            # (i+1) coming from the vertex marked j
            gr = xgraph()
            for v in vertices:
                gr.add_node(v)
            for i in range(r):
                for j in vertices:
                    gr.add_edge(j,match[0][i][j],i+1)
            if gr.is_connected():
                yield gr

def graphs(r,k):
    """
    This is an iterator, which spits out xgraphs which k-fold
    cover the wedge of r circles.
    """
    vertices = range(k)
    # To make the graph, we need r (not necessarily distinct)
    # permutations of k elements, permuted somehow.  We won't miss any
    # graphs if we just choose one per conjugacy class on the first
    # permutation.
    for perm1 in perm_cc(vertices):
	for otherperms in product(permutations(vertices),repeat=r-1):
	    perms = (perm1,)+otherperms
	   # for match in permutations(perms):
		# match[0][i][j] tells us the target of the edge labeled
		# (i+1) coming from the vertex marked j
	    gr = xgraph()
	    for v in vertices:
		gr.add_node(v)
	    for i in range(r):
		for j in vertices:
		    gr.add_edge(j,perms[i][j],i+1)
	    if gr.is_connected():
		yield gr




