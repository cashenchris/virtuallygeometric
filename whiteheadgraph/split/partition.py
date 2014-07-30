import networkx as nx

class Partition(object):
    """ 
    A tuple of disjoint nonempty frozensets.
    """

    def __init__(self,listofsetsorlists=[set([])],disjoint=False):
        tupleoffrozensets=tuple([frozenset(x) for x in listofsetsorlists if set(x)!=set([])])
        assert(len(frozenset.union(*tupleoffrozensets))==sum(len(x) for x in tupleoffrozensets)) # Assert the parts are disjoint
        self.parts=tupleoffrozensets          
            
    def __repr__(self):
        return str(self.parts)

    def __hash__(self):
        return hash(self.parts)
    
    def __eq__(self,other):
        if len(self.parts)!=len(other.parts):
            return False
        else:
            return all([self.parts[i]==other.parts[i] for i in range(len(self.parts))])
        
    def elements(self):
        return set.union(*[set(x) for x in self.parts])

    def which_part(self, element):
        return filter(lambda x: element in self.parts[x], range(0,len(self.parts)))[0]


def partition_union(*args):
    """
    Finest partition that coarsens a sequence of partitions.
    """
    def partition_union2(P1,P2):
        """
        Finest partition that coarsens both P1 and P2.
        """
        connectiongraph=nx.Graph()
        for p in P1.parts:
            connectiongraph.add_star([(1,e) for e in p])
        for p in P2.parts:
            connectiongraph.add_star([(2,e) for e in p])
        connectiongraph.add_edges_from([((1,e),(2,e)) for e in P1.elements()])
        components=[comp for comp in nx.connected_components(connectiongraph)]
        return Partition([set([comp[i][1] for i in range(len(comp))]) for comp in components])
    return reduce(partition_union2,args)
    
def make_partition(*args):
    """
    Finest partition that contains each set of args in a part.
    """
    totalspace=set.union(*[set(arg) for arg in args])
    partitionlist=[]
    for arg in args:
        partslist=[set(arg)]
        partslist.extend([ set([i]) for i in totalspace-set(arg) ])
        partitionlist.append(Partition(partslist))
    return partition_union(*partitionlist)
        

def compatible(P1,P2,splicemap):
    """
    Check if splicemap maps parts of P1 into parts of P2
    """
    
    return all([len(set([P2.which_part(splicemap[x]) for x in p]))==1 for p in P1.parts])

def is_part_bijection(P1,P2,splicemap):
    """
    Check if splicemap maps parts of P1 bijectively to parts of P2.
    """
    # First clause says parts map into parts.
    # Second clause says map is 1 to 1 on parts.
    # Third clause says map is onto on parts.
    return all([len(set([P2.which_part(splicemap[x]) for x in p]))==1 for p in P1.parts]) and all([len(set([P1.which_part(x) for x in P1.elements() if P2.which_part(splicemap[x])==i]))==1  for i in range(len(P2.parts))]) and all([len(set([x for x in P1.elements() if P2.which_part(splicemap[x])==i]))>0 for i in range(len( P2.parts))])

def is_reordered_partition(P1,P2):
    """
    Decide if P1 and P2 are same up to reordering parts.
    """
    fP1=set([frozenset(x) for x in P1.parts])
    fP2=set([frozenset(x) for x in P2.parts])
    return fP1==fP2

def compatible_coarsenings(P1,P2,partmap,splicemap):
    """
    Finest coarsenings of P1 and P2 compatible with splicemap and partitionmap.
    """
    newP1=Partition([p for p in P1.parts])
    newP2=Partition([p for p in P2.parts])
    while not compatible(newP1,newP2,splicemap):
        connectiongraph=nx.Graph()
        connectiongraph.add_edges_from([((1,newP1.which_part(i)),(2,newP2.which_part(splicemap[i]))) for i in range(0,len(splicemap))])
        components=[comp for comp in nx.connected_components(connectiongraph)]
        newP2=Partition([frozenset.union(*[newP2.parts[comp[i][1]] for i in filter(lambda x: comp[x][0]==2, range(0,len(comp)))]) for comp in components])
        newP1=Partition([frozenset.union(*[newP1.parts[partmap[comp[i][1]]] for i in filter(lambda x: comp[x][0]==2, range(0,len(comp)))]) for comp in components])
    return (newP1,newP2)

def partcd(P1,partsmap,P2,coarsemap2):
    """
    Compute bottom left entry of commutative partition diagram given top and right sides.
    """
    imageparts=len(set(coarsemap2))
    newpartslist=[[] for i in range(imageparts)]
    for j in range(len(P1.elements())):
        newpartslist[coarsemap2[partsmap[P1.which_part(j)]]]+=[j]
    return Partition(newpartslist)
