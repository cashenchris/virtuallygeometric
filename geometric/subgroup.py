from __future__ import division
import copy
from word import *
from stallings import *

def cyclicreduce(letters):
    if letters[0]+letters[-1]==0:
        return cyclicreduce(letters[1:-1])
    else:
        return letters


class conjclass(word):
    def __init__(self,w):
        if type(w)==word:
            shortest = word(cyclicreduce(w.letters))
        else:
            shortest = word(cyclicreduce(word(w).letters))
        testword = shortest
        for i in range(len(shortest.letters)):
            firstletter = word(testword.letters[:1])
            testword = firstletter**-1 * testword * firstletter
            if testword<shortest:
                shortest = testword
            self.letters = shortest.letters

                
#####################################################################
# New section on subgroups.
#####################################################################


class subgroup(object):
    """
    Subgroup of a free group, initialized by list of words.  
    Main data is a graph, obtained by folding.  Goal is to be
    able to define method subgroup.lift(word), which is supposed to
    give a list of words (or conjugacy classes?)  describing the lift
    of a loop to a cover.  Attached to a subgroup will be a rooted
    graph and a collection of generators.
    """
    def __init__(self,wordlist_or_graph):
        if str(type(wordlist_or_graph)).find('xgraph')>=0:
            graphcopy = copy.copy(wordlist_or_graph)
            subgroup.graph = graphcopy
        else:
            cover = xgraph(wordlist_or_graph)
            cover.fold()
            subgroup.graph = cover

    def __copy__(self):
        newgraph = copy.copy(self.graph)
        return subgroup(newgraph)
    
    def __deepcopy__(self):
        return self.__copy__()

    def is_finite_index(self,rank):
        """
        Must specify rank of ambient free group.  If rank is less than
        the number of generators appearing in the subgroup, the answer
        may be wrong.
        """
        answer = True
        for vertex in self.graph.nodes():
            if len(self.graph.neighbors(vertex))!=2*rank:
                answer = False
                break
        return answer

    def rank(self):
        v = len(self.graph.nodes())
        e = len(self.graph.edges())//2
        return e-v+1

    def index(self,rank):
        """
        Must specify rank of ambient free group.  If rank is less than
        the number of generators appearing in the subgroup, the answer
        may be wrong.
        """
        if self.is_finite_index(rank):
            return len(self.graph.nodes())
        else:
            return 'infinity'

    def lifts(self,w):
        """
        Return a list of words in a free generating set for the
        subgroup.  These words represent curves whose union is the
        preimage of the curve represented by w in the initial group.
        """
        roots = self.graph.nodes()[:] 
        tree = max_tree(self.graph)
        ## find the edges not in the spanning tree, and give them labels
        treeedges = set(tree.edges())
        alledges = set(self.graph.edges())
        newedges = list(alledges - treeedges)
        newgens={}; lab = 1
        for k in range(len(newedges)):   # This loop populates the dictionary (edge not in spanning tree) -> (nonzero integer)
            if newedges[k][2]>0:
                newgens[newedges[k]]=lab
                newgens[(newedges[k][1],newedges[k][0],-newedges[k][2])]=-(lab)
                lab=lab+1
        lifts_list = []
        while roots!=[]:  # We need to have a representative starting at an arbitrary vertex.
            root = roots.pop()
            current = root
            eat = w**1  # raise to power 1 to make sure we are copying and not disturbing w
            lift = word([])
            while len(eat)!=0 or current!=root:
                if len(eat)==0:  # in this case, we went through the word, but didn't end up back at root.
                    eat = w**1   # start eating another copy of w.
                    del roots[roots.index(current)]  # If we were to start at current,
                    # we'd just end up with a conjugate word, so we don't do it.
                else:
                    letter = eat.pop()
                    # now look for edge in the edges from current
                    edges = self.graph.dict[current]
                    for edge in edges:
                        if edge[1]==letter:
                            # we've found the right edge, so we traverse it, adding a letter to lift if it's not in the tree.
                            tup = (current,edge[0],letter)
                            if tup in newedges:
                                lift = lift*word([newgens[tup]])
                            current = edge[0]
                            break
            lifts_list.append(lift)
        return lifts_list

                            
                    
