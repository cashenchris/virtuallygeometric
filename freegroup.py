from numpy import sign
import random
import copy
from fractions import gcd
from numpy import sqrt
from group import *
import networkx as nx
import itertools
import pylab
from numpy import prod
#from fish import ProgressFish


def ishomologicallytrivial(w):
    reorder = copy.copy(w.letters)
    reorder.sort()
    reorder = freereduce(reorder)
    if len(reorder)>0:
        return False
    else:
        return True

def shortlexcompare(l1,l2):
        " shortlex order on sequences of integers "
        if l1 == l2:
            return 0
        elif len(l1)!=len(l2):
            return len(l1)-len(l2)
        else:
            # This is probably too clever.  Also it might be better to order differently.
            return -2*int(l1<l2)+1

class FGFreeGroup(FPGroup):
    """
    A finitely generated free group.
    """
    def __init__(self,gens=[], **kwargs):
        if gens:
             self.rank=len(gens)
        else:
             self.rank=kwargs['numgens']
        FPGroup.__init__(self,gens, [], **kwargs)

    #def __repr__(self):
    #    return "< "+", ".join([str(g) for g in self.gens])+" | "+", ".join([self.word(r)() for r in self.rels])+" >"


    def is_free(self):
        return True
    
    def word__ne__(self,v, w):
        return v.letters!=w.letters
        
    def word__eq__(self,v, w):
        return v.letters==w.letters

    def word__cmp__(self,v, w):
        return shortlexcompare(v.letters,w.letters)

    def word__hash__(self,w):
        return tuple(w.letters).__hash__()

    def randomcommutator(F,length,stopafter=10000):
        rndwd = F.randomword(length)
        for n in range(stopafter):
            if ishomologicallytrivial(rndwd):
                #print('Found commutator on '+str(n)+'th try.')
                return rndwd
            else:
                rndwd = F.randomword(length)
        #print('Failed to find commutator.')
        return []


    def guessmean(F,length):
    # Here's some experimental verification that the expected progress of
    # a random walk of length n in a free group of rank k is about
    # (1-1/k)n.  (Not quite right for small values of n.)
        result = []
        for n in range(100):
            result.append(F.randomword(length).letters)
        result = map(len,result)
        return sum(result)/100



    def are_conjugate(F,v,w):
        v1=copy.copy(F.cyclic_reduce(v).letters)
        w1=copy.copy(F.cyclic_reduce(w).letters)
        equal=False
        if len(v1)!=len(w1):
            return equal
        else:
            for i in range(len(v1)):
                if v1==w1:
                    equal=True
                    break
                temp = v1.pop()
                v1.insert(0, temp)
        return equal

    def get_conjugator(F,vi,wi):
        """
        Find an element u of a free group F such that inverse(u)*vi*u=wi. Return None if not conjugate.
        """
        v=vi**(1)
        w=wi**(1)
        uprefix=F.word([])
        umiddle=F.word([])
        usuffix=F.word([])
        if v.letters==w.letters:
            return uprefix*umiddle*usuffix
        elif len(v)==0 or len(w)==0:
            return None
        else:    
            while v.letters[0]==-v.letters[-1]:
                x=F.word(v.letters[-1:])
                uprefix=x**(-1)*uprefix
                v=v.cycle(-1)
            while w.letters[0]==-w.letters[-1]:
                y=F.word(w.letters[-1:])
                usuffix=y*usuffix
                w=w.cycle(-1)
            if len(w)!=len(v):
                return None
            else:
                for i in range(len(v)):
                    if v.letters==w.letters:
                        return uprefix*umiddle*usuffix
                    else:
                        z=F.word(v.letters[-1:])
                        v=v.cycle(-1)
                        umiddle=umiddle*z**(-1)
                else:
                    if v.letters==w.letters:
                        return uprefix*umiddle*usuffix
                    else:
                        return None
                
    def is_power(F,v,w):
        """
        Decide if v is a power of w.
        """
        v1=copy.copy(v)
        w1=copy.copy(w)
        w1bar=w1**(-1)
        pos=F.word([])
        neg=F.word([])
        ispower=False
        while len(v1)>=len(pos):
            if v1.letters==pos.letters or v1.letters==neg.letters:
                ispower=True
                break
            pos=pos*w1
            neg=neg*w1bar

        return ispower



    def is_conjugate_into_one(F,v,w):
        """
        Decide if v is conjugate into <w>.
        """
        v0,v1=F.cyclic_reducer(v)
        w0,w1=F.cyclic_reducer(w)
        x=F.word([])
        power=1
        if len(v1)==0:
            conjinto=True
            sgn=1
        elif len(w1)==0:
            conjinto=False
        elif len(v1)%len(w1):
            conjinto=False
        else:
            power=len(v1)//len(w1)
            sgn=1
            x=F.get_conjugator(v1,w1.__pow__(sgn*power))
            if x is None:
                sgn=-1
                x=F.get_conjugator(v1,w1.__pow__(sgn*power))
            if x is None:
                conjinto=False
            else:
                conjinto=True
        if conjinto:
            return dict([('conjugator',v0**(-1)*x*w0),('power',sgn*power)])
        else:
            return dict([])


        
            
        

    def is_conjugate_into(F,v, *wordlist):
        """
        Decide if v is conjugate into <w> for some w in wordlist.
        """
        for i in range(len(wordlist)):
            CIO=F.is_conjugate_into_one(v,wordlist[i])
            if CIO:
                CIO['index']=i
                return CIO
        return dict([])



    def is_subword(F,v,w,orientable=True):
        """
        Decide if v is a subword of w.
        """
        lv = len(v.letters)
        lw = len(w.letters)
        if orientable:
            return any((v.letters == w.letters[i:i+lv]) for i in xrange(len(w.letters)-lv+1))
        else:
            vbackward=[i for i in reversed(v.letters)]
            return any((v.letters == w.letters[i:i+lv] or vbackward == w.letters[i:i+lv]) for i in xrange(len(w.letters)-lv+1))


    def abelianization(F,w):
        """
        Compute abelianzation of a word w in a free group F.
        """
        if len(w.letters)==0:
            return []
        else:
            powers=[0]*F.rank
            lets=copy.copy(w.letters)
            while lets:
                l=lets.pop()
                powers[abs(l)-1]=powers[(abs(l)-1)]+sign(l)
            return powers


    def max_root(F,w,uptoconjugacy=False, withpower=True):
        """
        Find an indivisible root of w in a free group F.
        """
        w1=F.cyclic_reduce(w)
        conjugator=F.word(w.letters[0:(len(w)-len(w1))//2])
        n=len(w1.letters)
        if n==0:
            theroot=w1
            thepower=0
        elif n==1:
            theroot=w1
            thepower=1
        else:
            abpowergcd=reduce(gcd,F.abelianization(w1))
            possiblerootlengths=[i for i in range(1,1+int(sqrt(n))) if (n%i==0 and abpowergcd%(n//i)==0)]+[n//i for i in range(int(sqrt(n)),0,-1) if (n%i==0 and abpowergcd%(i)==0)]
            for i in possiblerootlengths:
                if w1.letters[0:i]*(n//i)==w1.letters:
                    theroot=F.word(w1.letters[0:i])
                    thepower=n//i
                    break
        if withpower:
            if uptoconjugacy:
                return theroot, thepower
            else:
                return conjugator*theroot*conjugator**(-1), thepower
        else:
            if uptoconjugacy:
                return theroot
            else:
                return conjugator*theroot*conjugator**(-1)
                
    def conjugate_root(F,w, withpower=False):
        """
        Find the conjugate of the indivisible root of w in F that is lex minimal.
        """
        root, power=F.max_root(w,True, True)
        invroot=root**(-1)
        lexmin=root
        if shortlexcompare(invroot.letters,lexmin.letters)<0:
            lexmin=invroot
        if len(root)>1:
            for i in range(1,len(root)):
                root=root.cycle(1)
                invroot=invroot.cycle(1)
                if shortlexcompare(root.letters,lexmin.letters)<0:
                    lexmin=root
                if shortlexcompare(invroot.letters,lexmin.letters)<0:
                    lexmin=invroot
        if withpower:
            return lexmin, power
        else:
            return lexmin
                

    def degree(F,w):
        """
        Max n such that w is an nth power.
        """
        return F.max_root(w,uptoconjugacy=True)[1]




#---------------------------
# subgroup stuff

def find_basis_by_folding(F,wl):
    """
    Find a basis for the subgroup of F generated by the words in wl.
    """
    basis=[]
    graph = StallingsGraph(wl)
    graph.fold()
    G=FGSubgroupOfFree(F,inclusionlist=[],graph=graph)
    inclusion=G.get_inclusion(F)
    return [inclusion(G.word([i])) for i in range(1,1+G.rank)]
    
def get_edge_path(tree,origin,terminus):
    """
    Find the unique geodesic edgepath joining orign to terminus in tree.
    """
    if origin==terminus:
        return []
    paths=[[edge] for edge in tree.out_edges(origin,keys=True, data=True)]
    while True:
        for path in paths:
            if path[-1][1]==terminus:
                return path
        else:
            newpaths=[path+[edge] for path in paths for edge in tree.out_edges(path[-1][1],keys=True,data=True)  if edge[3]['superlabel']!=-path[-1][3]['superlabel']]
            paths=newpaths
        
def get_tree_edge_path(graph,origin,terminus):
    """
    Get unique edgepath joining origin to terminus consisting of edges of max subtree  marked by data 'treelabel'=0.
    Assume graph is connected and max subtree is marked.
    """
    if origin==terminus:
        return []
    paths=[[edge] for edge in graph.out_edges(origin,keys=True, data=True) if edge[3]['treelabel']==0]
    while True:
        for path in paths:
            if path[-1][1]==terminus:
                return path
        else:
            newpaths=[path+[edge] for path in paths for edge in graph.out_edges(path[-1][1],keys=True,data=True)  if edge[3]['superlabel']!=-path[-1][3]['superlabel'] and edge[3]['treelabel']==0 ]
            paths=newpaths
 
    
class FGSubgroupOfFree(FGSubgroup, FGFreeGroup):
    """
    A finitely generated free subgroup of a free group.
    Initialization:
    supergroup = the free supergroup

    inclusionlist = optional list of words in the supergroup that are a basis for the subgroup. It is an error if this is not a basis. Use FGSubgroupOfFreeFrom for a wordlist that is not necesarily a basis.

    Alternately, group can be initialized from a StallingsGraph given with keyword 'graph'.

    If no graph and no inclusionlist are given the subgroup is trivial.

    If both an inclusionlist and a graph are given the inclusionlist is deemed to be the preferred basis.
    
    Part of the data that the subgroup carries is a marked Stallings graph as self.graph. Edges have data 'treelabel' and 'superlabel'. The 'superlabel' is the corresponding basis element of the supergroup. The 'treelabel' is 0 if the edge belongs to a fixed maximal subtree. The remaining edges have treelabel in +-{1,...,rank}. If edge e has treelabel i this means that we have chosen a basis for the fundametal group of the graph so that the i-th element is the loop consisting of the unique segment in the max subtree from vertex 0 to the origin of e, then the edge e, then the unique segment in the maximal subtree from the terminus of e back to 0.

    marking[i]['initial'] is unique oriented edgepath in the marked maximal subtree from the root vertex to the origin vertex of the edge with treelabel=i+1
    marking[i]['final'] is unique oriented edgepath in the marked maximal subtree from the terminal vertex of the edge with treelabel=i+1 to the root vertex
    inversemarking[i] is the word of the subgroup (written in terms of the preferred basis for the subgroup) corresponding to the edge with treelabel=i+1 

    If an inclusionlist was given at initialization then the prefered basis is the inclusionlist, and the graphbasis may be different. There are conversion methods graphBasisToGivenBasis and givenBasisToGraphBasis to go back and forth.
    If no inclusionlist is given then the prefered basis is set to be equal to the graph basis.

    """
    def __init__(self, supergroup, inclusionlist=[], **kwargs):
        try:
            verbose=kwargs['verbose']
        except KeyError:
            verbose=False
        self.supergroup=supergroup
        if inclusionlist:
            if verbose:
                print "Received input inclusionlist"
            FGSubgroup.__init__(self,supergroup,inclusionlist,**kwargs)
            marking=[]
            inversemarking=[]
            if verbose:
                print "Constructing Stallings Graph"
            cover= StallingsGraph(inclusionlist, 0, marking,inversemarking, self, verbose)
            cover.fold(marking, inversemarking, verbose=verbose)
            self.graph = cover
            self.inversemarking=inversemarking
            self.rank=1+len(self.graph.edges())//2 - len(self.graph.nodes())
            if self.rank!=len(inclusionlist):
                    raise ValueError("Rank has dropped. Input was not a basis.")
        else:
            try:
                self.graph=kwargs['graph']          
            except KeyError: # no graph given and no inclusionlist means trivial group
                self.rank=0
                self.graph=StallingsGraph()
                self.inversemarking=[]
                FGSubgroup.__init__(self,supergroup,inclusionlist,**kwargs)
            else:
                if verbose:
                    print "Received input Stallings Graph"
                self.rank=1+len(self.graph.edges())//2 - len(self.graph.nodes())
                try: # see if the graph has a maximal subtree and generators labeled
                    for e in self.graph.out_edges(self.graph.basepoint,data=True):
                        e[2]['treelabel']
                        break
                except KeyError: # treelabels not set, set them
                    if verbose:
                        print "Marking max tree in Stallings Graph"
                    mark_max_tree(self.graph, verbose=verbose)
                self.graph.fold(verbose=verbose)
                # now want to write marked loops in the graph as words in the supergroup.
                if verbose:
                    print "Construcing inclusionlist"
                spathtov=super_path_to_vertex(self.graph) # dict spathtov[v]=superletters along edge path from basepoint to v through marked maxtree
                for i in range(1,1+self.rank): 
                    for edge in self.graph.edges(keys=True,data=True):# find the edge that's the i-th basis element of the subgroup
                        if self.graph[edge[0]][edge[1]][edge[2]]['treelabel']==i:
                            break
                    else:
                        raise KeyError(i)
                    thesuperlettersforthisloop=spathtov[edge[0]]+[edge[3]['superlabel']]+[-x for x in reversed(spathtov[edge[1]])]
                    inclusionlist.append(self.supergroup.word(thesuperlettersforthisloop))
                FGSubgroup.__init__(self,supergroup,inclusionlist,**kwargs)
                self.inversemarking=[self.word([i]) for i in range(1,1+self.rank)]      
        

    def super_letters_to_graph_basis(self, wordorletters, startingvert=None):
        """
        Take a word in the supergroup that is an element of the subgroup and write it in terms of the graph basis.
        """
        if startingvert is None:
            startingvert=self.graph.basepoint
        try: # see if the graph has a maximal subtree and generators labeled
            for e in self.graph.out_edges(self.graph.basepoint,data=True):
                e[2]['treelabel']
                break
        except KeyError: # treelabels not set, set them
            mark_max_tree(self.graph)
        letters=self.supergroup.word(wordorletters).letters
        graphletters=[]
        currentvert=startingvert
        while letters:
            nextletter=letters.pop(0)
            for edge in self.graph.out_edges(currentvert,keys=True,data=True):
                if edge[3]['superlabel']==nextletter:
                    if edge[3]['treelabel']: # if this is not an edge of the marked maximal subtree
                        graphletters.append(edge[3]['treelabel'])
                    currentvert=edge[1]
                    break
            else:
                return None
        if currentvert!=startingvert:
            return None
        return graphletters

    def given_basis_to_graph_basis(self, wordorletters):
        """
        Take a word in the given basis and write it in terms of the graph basis.
        """
        return self.super_letters_to_graph_basis(self.inclusion(self.word(wordorletters)))

    def graph_basis_to_given_basis(self, graphletters):
        """
        Take a list of integers from {-rank..-1}u{1,..rank} interpreted as a loop in the graph and return the corresponding word with repsect to the given basis.
        """
        theletters=[x for x in graphletters]
        theword=self.word([])
        while theletters: # translate graphletters to groupletters
            nextletter=theletters.pop(0)
            if nextletter>0:
                theword=theword*self.inversemarking[nextletter-1]
            elif nextletter<0:
                theword=theword*(self.inversemarking[-nextletter-1])**(-1)
            else:
                raise TypeError
        return theword      

    def super_letters_to_given_basis(self, wordsorletters):
        """
        Take a word in the supergroup that is an element of the subgroup and write it in terms of the given basis for the subgroup.
        """
        return self.graph_basis_to_given_basis(self.super_letters_to_graph_basis(wordsorletters))
    
    def is_finite_index(self):
        return self.index()<float('inf')
    
    def index(self):
        """
        Return index of self in supergroup.
        """
        iscover=True
        for v in self.graph.nodes_iter():
            if len(self.graph.out_edges(v))<2*self.supergroup.rank:
                iscover=False
                break
        if not iscover:
            return float('inf')
        else:
            return len(self.graph.nodes())
    

    def lifts(self,w, verbose=False):
        """
        Return a list of words in a free generating set for the
        subgroup.  These words represent curves whose union is the
        preimage of the curve represented by w in the initial group.
        """
        if verbose:
            print "Finding lifts of word "+w()
        theletters=[i for i in w.letters]
        theSG=self.graph
        verts=set(theSG.nodes())
        liftsInSubgroup=[]
        vertssofar=0
        if verbose:
            numverts=len(verts)
            #fish=ProgressFish(total=numverts)
        while verts:
            thisvert=verts.pop()
            thisword=[]
            nextvert, nextgraphletters=followPath(theSG,theletters,thisvert, withgraphletters=True)
            assert(nextvert is not None) # Should not fail if self is a finite index subgroup.
            thisword+=nextgraphletters
            vertssofar+=1
            while nextvert!=thisvert:
                verts.remove(nextvert)
                nextvert, nextgraphletters=follow_path(theSG,theletters,nextvert,withgraphletters=True)
                assert(nextvert is not None)
                vertssofar+=1
                thisword+=nextgraphletters
            liftsInSubgroup.append(self.graph_basis_to_given_basis(thisword))
            #if verbose:
                #fish.animate(amount=vertssofar)
        return liftsInSubgroup


    def contains_word(self,w):
        """
        Decide if a word w defined in a supergroup is contained in a subgroup self.
        """
        theword=self.restrictWord(w)
        if theword is None:
            return False
        else:
            return True
        
    def contains_subgroup(self,G):
        """
        Decide if a subgroup G of a supergroup is contained in a subgroup self.
        """
        F=common_ancestor(self,G)
        if F is None:
            return False
        return all([self.contains_word(w) for w in [G.get_inclusion(F)(G.word([i])) for i in range(1,1+G.rank)]])
        
    def equivalent_subgroup(self,G):
        """
        Decide if G is the same subgroup as self.
        """
        return self.contains_subgroup(G) and G.contains_subgroup(self)
    
    # Be careful about defining comparison operators here. Our subgroups come equipped with a preferred generating set.
    # equivalentSubgroup(G,H)=True if G and H define the same (in group theoretic sense) subgroup of some common supergroup.
    # G is H =True if G and H are the same object for python.
    # We could also call G and H "the same" if the inclusion into a common supergroup matched up the preferred generating sets.
    # What ought G==H to mean?


    def restrict_word(self,w):
        """
        If w and self are contained in a common free group try to write w as a word in self. Return None if w is not in subgroup
        """
        if w.group is self:
            return w
        elif w.group is self.supergroup: # simple one step restriction
            graphletters=self.super_letters_to_graph_basis(w) # translate w to graph basis
            if graphletters is None:
                return None
            return self.graph_basis_to_given_basis(graphletters)
        else:  
            G=common_ancestor(self,w.group)
            wgpinG=w.group.get_inclusion(G)
            winG=wgpinG(w)
            superchain=[self]
            while superchain[-1] is not G:
                superchain.append(superchain[-1].supergroup)
            return self.restrict_word(superchain[-2].restrict_word(winG)) # inner call is a simple one step restriction
        
        
    def find_conjugate_in_subgroup(self,w):
        """
        Return a conjugate of w in self, or None if none exists.
        """
        if w.group is self:
            return self.cyclic_reduce(w)
        else:
            G=common_ancestor(self,w.group)
            winG=G.cyclic_reduce(w.group.get_inclusion(G)(w))
            newself=FGSubgroupOfFreeFrom(G,self) # this is a copy of self as an immediate subgroup of G, possibly with a different generating set
            for v in newself.graph:
                graphletters=newself.super_letters_to_graph_basis(winG,v)
                if graphletters is not None: # winG makes a loop at v, found a conjugate of w in newself
                    conj=newself.word(newself.graph_basis_to_given_basis(graphletters))
                    return self.restrict_word(conj)
            else:
                return None
            
        
        

        


    
def FGSubgroupOfFreeFrom(F,wlorsg,**kwargs):
    """
    Take a list of words or a subgroup of F and return subgroup with basis.
    """
    if hasattr(wlorsg,'gens'):
        inclusion=wlorsg.get_inclusion(F)
        wl=[inclusion(wlorsg.word([i])) for i in range(1,1+len(wlorsg.gens))]
    else:
        wl=wlorsg
    thisgraph=StallingsGraph(wl)
    thisgraph.fold()
    mark_max_tree(thisgraph)
    return FGSubgroupOfFree(F,inclusionlist=[],graph=thisgraph,**kwargs)


class StallingsGraph(nx.MultiDiGraph):
    """
    A Stallings Graph for a finitely generated subgroup of a free group.
    Each edge has a unique key. Inverse edges have keys differing by sign.
    Each edge has data 'superlabel' that is +-integer representing a generator of the supergroup.
    Each edge optionally carries data 'treelabel'. A chosen maximal subtree can be set to have treelabels=0, and remaining edges are given treelabel +- 1..rank of subgroup.

    Initialization:
    wordlist is a list of words in a free group generating the subgroup.
    Optionally specify the name of the basepoint, default=0.

    marking and inversemarking and subgroup are for tracking the relationship between the preferred basis for the subgroup and the one that comes from the complement of a max tree in the Stallings Graph. These are used by FGSubgroupOfFree.__init__().
    """
    def __init__(self,wordlist=[], basepoint=0, marking=None, inversemarking=None, subgroup=None, verbose=False):
        nx.MultiDiGraph.__init__(self)
        self.basepoint=basepoint
        self.add_node(self.basepoint)
        counter=1
        if verbose:
            lengthwordlist=len(wordlist)
            #fish=ProgressFish(total=lengthwordlist)
        for i in range(1,1+len(wordlist)):
            if verbose:
                pass
                #fish.animate(amount=i)
            w=wordlist[i-1]
            try:
                inversemarking.append(subgroup.word([i]))
            except AttributeError:
                pass
            currentvert=self.basepoint
            newletters=[x for x in w.letters]
            try:
                marking.append(dict([('initial',[]),('final',[])]))
            except AttributeError:
                pass
            while len(newletters)>1:
                nextletter=newletters.pop(0)
                self.add_node(counter)
                try:
                    k=1+max([edge[2] for edge in self.edges(keys=True)])
                except ValueError: # there aren't any edges yet
                    k=1
                self.add_edge(currentvert,counter,superlabel=nextletter, treelabel=0, key=k)
                try:
                    marking[i-1]['initial'].append((currentvert,counter,k))
                except TypeError:
                    pass
                currentvert=counter
                counter+=1
            if newletters:# word was nonempty
                nextletter=newletters.pop(0)
                self.add_edge(currentvert,self.basepoint,superlabel=nextletter, treelabel=i)
            else: #word was empty
                pass
    

    def add_edge(self, origin, terminus, **kwargs):
        data=dict([(x,kwargs[x]) for x in kwargs if x!='key'])
        if 'key' in kwargs:
            key=kwargs['key']
        elif not self.edges():
            key=1
        else:
            key=1+max([edge[2] for edge in self.edges(keys=True)])
        nx.MultiDiGraph.add_edge(self,origin,terminus, key, **data)
        # now add the inverse edge, reversing key and treelabel, superlabel, if they exist.
        try:
            data['treelabel']=-data['treelabel']
        except KeyError:
            pass
        try:
            data['superlabel']=-data['superlabel']
        except KeyError:
            pass
        nx.MultiDiGraph.add_edge(self,terminus,origin, -key, **data)

    def del_edge(self,origin, terminus, superlabel): # delete an arbitrary edge with given superlabel and its inverse
        for k in self[origin][terminus]: # find a suitable key k
            if self[origin][terminus][k]['superlabel']==superlabel:
                break
        else:
            raise KeyError
        nx.MultiDiGraph.remove_edge(self,origin,terminus,k)
        nx.MultiDiGraph.remove_edge(self,terminus,origin,-k)

    def delete_edge(self,edge): # delete an edge and its inverse
        nx.MultiDiGraph.remove_edge(self,edge[0],edge[1],edge[2])
        nx.MultiDiGraph.remove_edge(self,edge[1],edge[0],-edge[2])

    def del_node(self,vert):
        nx.MultiDiGraph.remove_node(self,vert)

    def clone_node(self,vert ,new_node=None):
        if new_node in self.nodes():
            raise KeyError(str(new_node)+" is already in the graph")
        if new_node is None:
            new_node=1+max(self.nodes())
        self.add_node(new_node)
        for edge in self.out_edges(vert,keys=True,data=True):
            if edge[1]==vert: # this edge is a loop at v
                if edge[2]>0:
                    self.add_edge(new_node,new_node,key=edge[2],**edge[3])
            else:
                self.add_edge(new_node,edge[1],key=edge[2],**edge[3])
        return new_node
        

    def switch_nodes(self,one,two):
        """
        switch two nodes, carrying edges along.
        """
        new = self.clone_node(one)
        self.del_node(one)
        self.clone_node(two,one)
        self.del_node(two)
        self.clone_node(new,two)
        self.del_node(new)
        

    def neighbors(self,vertex):
        """
        return as list
        """
        return list(self.distinct_neighbors(vertex))

    def distinct_neighbors(self,vertex):
        """
        return as set
        """
        return set([edge[1] for edge in self.out_edges(vertex)]+[edge[0] for edge in self.in_edges(vertex)])

    def posedges(self):
        return [edge[:3] for edge in self.edges(keys=True, data=True) if edge[3]['superlabel']>0]

    def fold(self, marking=None, inversemarking=None, verbose=False):
        try: # see if the graph has a maximal subtree and generators labeled
            for e in self.out_edges(self.basepoint,data=True):
                e[2]['treelabel']
                break
        except KeyError: # treelabels not set, set them
            if verbose:
                print "Marking max tree"
            mark_max_tree(self, verbose=verbose)
        if verbose:
            #fish=ProgressFish() # progress counter
            print "Folding. Number of vertices remaining:"
        while foldonce(self, marking, inversemarking)!=None:
            if verbose:
                pass
                #fish.animate(self.number_of_nodes())

    def core(self, verbose=False):
        """
        Change self to be core.
        """
        # Pass to connected component containing basepoint and trim off any valence 1 vertices other than basepoint.
        # Need this, for instance, in finding Stallings Graph of intersection of two subgroups.
        simplegraph=nx.Graph(self) # An undirected copy of self without any multiedges.
        if verbose:
            print "Finding connected component of the basepoint"
        simplecoregraph=simplegraph.subgraph(nx.node_connected_component(simplegraph,self.basepoint)) # list vertices in the same component as basepoint.
        if verbose:
            print "Removing valence 1 vertices"
            totalverts=simplecoregraph.number_of_nodes()
            #fish=ProgressFish(total=totalverts)
            deletedverts=0
        vertstoremove=[v for (v,d) in simplecoregraph.degree_iter() if v!=self.basepoint and d==1]
        while vertstoremove:
            while vertstoremove:
                nextvert=vertstoremove.pop()
                simplecoregraph.remove_node(nextvert)
                if verbose:
                    deletedverts+=1
                    #fish.animate(amount=deletedverts)
            vertstoremove=[v for (v,d) in simplecoregraph.degree_iter() if v!=self.basepoint and d==1]
        vertsinthecore=set(simplecoregraph.nodes())
        self.remove_nodes_from([v for v in self if v not in vertsinthecore])

    def draw(self):
        # Note: Networkx draws all edges as straight lines. This doesn't work so well for multieges or loops.
        # To avoind this problem we subdivide edges and label the mid-edge vertex with the edgelabel.
        positivegraph=nx.MultiDiGraph()
        labels=dict()
        for v in self.nodes_iter():
            positivegraph.add_node(v)
            labels[v]=v
        for e in self.edges_iter(data=True):
            if e[2]['superlabel']>0:
                positivegraph.add_edge(e[0],(e[0],e[1],e[2]['superlabel']))
                positivegraph.add_edge((e[0],e[1],e[2]['superlabel']),e[1])
                labels[(e[0],e[1],e[2]['superlabel'])]=e[2]['superlabel']
        pylab.figure(1)
        pos=nx.spring_layout(positivegraph)
        nx.draw_networkx_nodes(positivegraph,pos,self.nodes())
        nx.draw_networkx_nodes(positivegraph,pos,set(positivegraph.nodes())-set(self.nodes()),node_color='w')
        nx.draw_networkx_labels(positivegraph,pos,labels)
        nx.draw_networkx_edges(positivegraph,pos)
        #edge_labels=dict([((u,v),d['superlabel']) for u,v,d in positivegraph.edges(data=True)])
        #nx.draw_networkx_edge_labels(positivegraph,pos,edge_labels=edge_labels)
        
def unfolded(theSG):
    "Return key for vertex which is not folded "
    for vertex in theSG.nodes():
        edges=theSG.out_edges(vertex,keys=True,data=True)
        for i in range(len(edges)-1):
            for j in range(i+1,len(edges)):
                if edges[i][3]['superlabel']==edges[j][3]['superlabel']:
        # We're going to delete one of the targets if they differ.  Make sure not to delete the source.
                    target1 = edges[i][1]; target2 = edges[j][1]
                    if target1!=target2 and target2==vertex:
                        return (vertex,edges[j],edges[i])
                    elif target1==target2 and edges[i][3]['treelabel']==0:
                        return (vertex,edges[j],edges[i])
                    else:
                        return (vertex,edges[i],edges[j])
    return None

def tighten(edgepath):
    """
    Remove backtracking from an edgepath
    """
    tightpath=[]
    while edgepath:
        nextedge=edgepath.pop(0)
        if not tightpath:
            tightpath.append(nextedge)
        else:
            if tightpath[-1][2]==-nextedge[2]: # backtracking
                tightpath.pop()
            else:
                tightpath.append(nextedge)
    return tightpath
    
    
    
def foldonce(theSG, marking=None, inversemarking=None):
    """
    Do a single fold (identify two edges) if possible.  Returns 1 if something happened. Returns 'None' otherwise.
    """
    # marking[i]['initial'] is unique oriented edgepath in the marked maximal subtree from the root vertex to the origin vertex of the edge with treelabel=i+1
    # marking[i]['final'] is unique oriented edgepath in the marked maximal subtree from the terminal vertex of the edge with treelabel=i+1 to the root vertex
    # inversemarking[i] is the word of the subgroup corresponding to the edge with treelabel=i+1 
    place_to_fold = unfolded(theSG)
    if place_to_fold is None:
        return None
    else:
        source,edge1,edge2 = place_to_fold
        target1 = edge1[1]
        target2 = edge2[1]
        foldletter = edge1[3]['superlabel']
        # First delete extra edge coming from source
        if marking: 
            # we've arrange in function unfolded that if edge2 is a loop then so is edge1, and if target1==target2 then if edge1 has treelabel==0  so does edge2
            if target1==target2:
                raise TypeError("Rank dropped. Input was not a basis.")
            elif edge2[0]==edge2[1]:
                raise RuntimeError("edge2 should not be a loop unless edge1 is. Something went wrong in unfolded.")
            elif edge2[3]['treelabel']!=0: # in this case edge2 is marked as a basis element. We don't want to lose it. Need to move this label somewhere else.

                # we've set things up so that marking[tl-1]['initial'])+edge2+len(marking[tl-1]['final'] consists of a leading segment followed by a simple loop followed by the leading segment backwards. edge2 is in the simple loop and is not the only edge in the simple loop. If initial segment has length <= final segment then the first edge of the final segment is part of the simple loop. Otherwise, the last edge of the initial segment is part of the simple loop.
                tl=abs(edge2[3]['treelabel'])
                if edge2[3]['treelabel']<0:
                    inversemarking[tl-1]=inversemarking[tl-1]**(-1)
                    tempmark=[x for x in marking[tl-1]['initial']]
                    marking[tl-1]['initial']=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['final'])]
                    marking[tl-1]['final']=[(e[1],e[0],-e[2]) for e in reversed(tempmark)]
                
                if len(marking[tl-1]['initial'])<=len(marking[tl-1]['final']):
                    newedge=marking[tl-1]['final'][0]
                else:
                    newedge=marking[tl-1]['initial'][-1]
                change_of_marking(theSG,marking,inversemarking,edge2,newedge)
                          
            # now edge2 is not a loop and has treelabel=0
            
            #make sure folding edge2 leaves us with a connected max subtree
            tl=abs(edge1[3]['treelabel'])
            if tl!=0 and edge1[0]!=edge1[1]:
                if sign(edge1[3]['treelabel'])<0: # for simplicity reverse this loop
                    inversemarking[tl-1]=(inversemarking[tl-1])**(-1)
                    dummymarking=[e for e in marking[tl-1]['initial']]
                    marking[tl-1]['initial']=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['final'])]
                    marking[tl-1]['final']=[(e[1],e[0],-e[2]) for e in reversed(dummymarking)]
                    theSG[edge1[0]][edge1[1]][edge1[2]]['treelabel']=tl
                    theSG[edge1[1]][edge1[0]][-edge1[2]]['treelabel']=-tl
                    
                newrepresentativeedge=edge1                  
                if marking[tl-1]['initial']:
                    if -edge2[2]==marking[tl-1]['initial'][-1][2]:
                        if marking[tl-1]['final']:
                            if marking[tl-1]['final'][0][2] in [-e[2] for e in marking[tl-1]['initial']]:
                                newrepresentativeedge=marking[tl-1]['initial'][-2]
                            else:
                                newrepresentativeedge=marking[tl-1]['final'][0]
                        else:
                            newrepresentativeedge=marking[tl-1]['initial'][-2]
                if -edge2[2] in [e[2] for e in marking[tl-1]['final']]:
                    newrepresentativeedge=marking[tl-1]['final'][0]
                if newrepresentativeedge!=edge1:
                    change_of_marking(theSG,marking, inversemarking,edge1,newrepresentativeedge)
                
                                       
                

            # if edge1 has treelabel==0 we can go ahead and replace edge2 by edge1 in the marking and tighten
            # if edge1 has treelabel!=0 we must reroute marking around edge2
            newmarking=[x for x in marking]
            newinversemarking=[x for x in inversemarking]
            if edge1[3]['treelabel']==0:
                for i in range(len(marking)):
                    for j in range(len(marking[i]['initial'])):
                        if marking[i]['initial'][j][2]==edge2[2]:
                            newmarking[i]['initial'][j]=(edge1[0],edge1[1],edge1[2])
                            newmarking[i]['initial']=tighten(marking[i]['initial'])
                            break
                        if marking[i]['initial'][j][2]==-edge2[2]:
                            newmarking[i]['initial'][j]=(edge1[1],edge1[0],-edge1[2])
                            newmarking[i]['initial']=tighten(marking[i]['initial'])
                            break
                    for j in range(len(marking[i]['final'])):
                        if marking[i]['final'][j][2]==edge2[2]:
                            newmarking[i]['final'][j]=(edge1[0],edge1[1],edge1[2])
                            newmarking[i]['final']=tighten(marking[i]['final'])
                            break
                        if marking[i]['final'][j][2]==-edge2[2]:
                            newmarking[i]['final'][j]=(edge1[1],edge1[0],-edge1[2])
                            newmarking[i]['final']=tighten(marking[i]['final'])
                            break

            else:
                tl=abs(edge1[3]['treelabel'])
                orientation=sign(edge1[3]['treelabel'])
                if orientation>0:
                    initialsegment=marking[tl-1]['initial']
                    finalsegment=marking[tl-1]['final']
                else:
                    initialsegment=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['final'])]
                    finalsegment=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['initial'])]
                # if edge1 is a loop edge2 could still lie on the initial or final segment
                if edge1[0]==edge1[1] and finalsegment:
                    if finalsegment[0][2]==edge2[2]:
                        assert(initialsegment[-1][2]==-edge2[2])
                        finalsegment.pop(0)
                        initialsegment.pop()
                
                for i in range(len(marking)):
                    for j in range(len(marking[i]['initial'])):
                        if marking[i]['initial'][j][2]==edge2[2]:
                           newmarking[i]['initial']=tighten([(e[1],e[0],-e[2]) for e in reversed(finalsegment)]+marking[i]['initial'][j+1:])
                           newinversemarking[i]=inversemarking[tl-1]**(-1*orientation)*inversemarking[i]
                           break
                        if marking[i]['initial'][j][2]==-edge2[2]:
                           newmarking[i]['initial']=tighten(initialsegment+marking[i]['initial'][j+1:])
                           newinversemarking[i]=inversemarking[tl-1]**(orientation)*inversemarking[i]
                           break
                    for j in range(len(marking[i]['final'])):
                        if marking[i]['final'][j][2]==edge2[2]:
                            newmarking[i]['final']=tighten(marking[i]['final'][:j]+[(e[1],e[0],-e[2]) for e in reversed(initialsegment)])
                            newinversemarking[i]=inversemarking[i]*inversemarking[tl-1]**(-1*orientation)
                            break
                        if marking[i]['final'][j][2]==-edge2[2]:
                            newmarking[i]['final']=tighten(marking[i]['final'][:j]+finalsegment)
                            newinversemarking[i]=inversemarking[i]*inversemarking[tl-1]**orientation
                            break
            for i in range(len(marking)):
                marking[i]=newmarking[i]
            for i in range(len(inversemarking)):
                inversemarking[i]=newinversemarking[i]
                 
                
        theSG.delete_edge(edge2)
        # If target1=target2, we are done.
        if target1==target2:
            return theSG
        # We are going to delete target2, but we don't want to delete the basepoint.
        if target2==theSG.basepoint:
            theSG.switch_nodes(target1,target2)
            dummy = target2
            target2 = target1
            target1 = dummy
        # Now redirect edges from target2
        for edge in theSG.out_edges(target2,keys=True,data=True):
            if edge[1]==target2: # for loops at target2
                theSG.add_edge(target1,target1,key=edge[2],**edge[3])
            else:
                theSG.add_edge(target1,edge[1],key=edge[2],**edge[3])
        theSG.del_node(target2)
        # finally, if we're tracking the marking we must replace target2 by target 1 in the marking
        if marking:
            for i in range(len(marking)):
                for j in marking[i].keys():
                    for k in range(len(marking[i][j])):
                        if marking[i][j][k][0]==target2:
                            marking[i][j][k]=(target1,marking[i][j][k][1],marking[i][j][k][2])
                        if marking[i][j][k][1]==target2:
                            marking[i][j][k]=(marking[i][j][k][0],target1,marking[i][j][k][2])
        return 1

def change_of_marking(theSG,marking,inversemarking,edge1,edge2):
    """
    Change the marking on the graph theSG so that edge2 represents the loop previously represented by edge1.
    It is necessary that edge2 does actually live on that loop!
    """
    tl=abs(theSG[edge1[0]][edge1[1]][edge1[2]]['treelabel'])
    assert(tl>0)
    orientation1=sign(theSG[edge1[0]][edge1[1]][edge1[2]]['treelabel'])
    # find edge2
    initial, final, orientation = False, False, None
    for i in range(len(marking[tl-1]['initial'])):
        if edge2[2]==marking[tl-1]['initial'][i][2]:
            initial=True
            orientation = 1
            initialsegment=marking[tl-1]['initial'][:i]
            finalsegment=marking[tl-1]['initial'][i+1:]+[edge1]+marking[tl-1]['final']
            break
        if -edge2[2]==marking[tl-1]['initial'][i][2]:
            initial=True
            orientation=-1
            initialsegment=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['final'])]+[(edge1[1],edge1[0],-edge1[2])]+[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['initial'][i+1:])]
            finalsegment=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['initial'][:i])]
            break
    for i in range(len(marking[tl-1]['final'])):
        if edge2[2]==marking[tl-1]['final'][i][2]:
            final=True
            orientation = 1
            initialsegment=marking[tl-1]['initial']+[(edge1[0],edge1[1],edge1[2])]+marking[tl-1]['final'][:i]
            finalsegment=marking[tl-1]['final'][i+1:]
            break
        if -edge2[2]==marking[tl-1]['final'][i][2]:
            final=True
            orientation=-1
            initialsegment=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['final'][i+1:])]
            finalsegment=[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['final'][:i])]+[(edge1[1],edge1[0],-edge1[2])]+[(e[1],e[0],-e[2]) for e in reversed(marking[tl-1]['initial'])]
            break
    assert(initial^final)# edge2 can't represent the loop if it is on both the initial and final segment, but must be on one of them
    theSG[edge2[0]][edge2[1]][edge2[2]]['treelabel']=tl
    theSG[edge2[1]][edge2[0]][-edge2[2]]['treelabel']=-tl
    theSG[edge1[0]][edge1[1]][edge1[2]]['treelabel']=0
    theSG[edge1[1]][edge1[0]][-edge1[2]]['treelabel']=0
    # now look through the markings. Whenever +- edge2 appears it must be replaced, since edge2 is no longer in the max subtree
    for i in range(len(marking)):
        if i==tl-1:
            marking[i]['initial']=initialsegment
            marking[i]['final']=finalsegment
        else:
            newmarking=[x for x in marking]
            newinversemarking=[x for x in inversemarking]
            for j in range(len(marking[i]['initial'])):
                if marking[i]['initial'][j][2]==edge2[2]:
                    newmarking[i]['initial']=tighten([(e[1],e[0],-e[2]) for e in reversed(finalsegment)]+marking[i]['initial'][j+1:])
                    newinversemarking[i]=(inversemarking[tl-1])**(-1*orientation1)*inversemarking[i]
                    break
                if marking[i]['initial'][j][2]==-edge2[2]:
                    newmarking[i]['initial']=tighten(initialsegment+marking[i]['initial'][j+1:])
                    newinversemarking[i]=(inversemarking[tl-1])**(orientation1)*inversemarking[i]
                    break
            for j in range(len(marking[i]['final'])):
                if marking[i]['final'][j][2]==edge2[2]:
                    newmarking[i]['final']=tighten(marking[i]['final'][:j]+[(e[1],e[0],-e[2]) for e in reversed(initialsegment)])
                    newinversemarking[i]=(inversemarking[i])*((inversemarking[tl-1])**(-1*orientation1))
                    break
                if marking[i]['final'][j][2]==-edge2[2]:
                    newmarking[i]['final']=tighten(marking[i]['final'][:j]+finalsegment)
                    newinversemarking[i]=(inversemarking[i])*(inversemarking[tl-1])**(orientation1)
                    break
            for i in range(len(marking)):
                marking[i]=newmarking[i]
            for i in range(len(inversemarking)):
                inversemarking[i]=newinversemarking[i]
        
        

def super_path_to_vertex(theSG):
    """
    Returns a dictionary  v:lists of superletters that describe the path through marked max subtree to vertex v.
    (If theSG represents a finite index subgroup these are the coset representatives.)
    """
    reps=dict()
    tree=theSG.copy()
    sphere=dict()
    for e in tree.edges_iter(keys=True,data=True):
        if e[3]['treelabel']!=0:
            tree.remove_edge(e[0],e[1],e[2])
        elif e[3]['superlabel']<0:
            tree.remove_edge(e[0],e[1],e[2])
            # What's left is a maximal subtree with all superlabels positive
    predecessor=dict()
    currentlevel=set([theSG.basepoint])
    sphere[0]=[theSG.basepoint]
    nextradius=1
    while currentlevel:
        nextlevel=set([])
        sphere[nextradius]=[]
        while currentlevel:
            thisvert=currentlevel.pop()
            newneighbors=tree.distinct_neighbors(thisvert)
            for neighbor in list(newneighbors):
                try:
                    edgekey=tree[thisvert][neighbor].keys()[0]
                    predecessor[neighbor]=(thisvert,tree[thisvert][neighbor][edgekey]['superlabel'])
                except KeyError:
                    edgekey=tree[neighbor][thisvert].keys()[0]
                    predecessor[neighbor]=(thisvert,-tree[neighbor][thisvert][edgekey]['superlabel'])
                sphere[nextradius].append(neighbor)
            nextlevel|=newneighbors
            tree.remove_node(thisvert)
        currentlevel|=nextlevel
        nextradius+=1
    reps[theSG.basepoint]=[]
    for radius in range(1,nextradius):
        for v in sphere[radius]:
            reps[v]=reps[predecessor[v][0]]+[predecessor[v][1]]
    return reps

def coset_reps(theSG):
    """
    Assuming theSG represents a finite index subgroup, returns a list of lists of superletters representing coset representatives for G\F
    """
    return super_path_to_vertex(theSG).values()
    
def follow_path(theSG,superletters,origin, withgraphletters=False):
    """
    Start at vertex origin. Follow edges labeled by superletters. Return terminal vertex or None given superletters do not form an edgepath.
    """
    currentvertex=origin
    theletters=[i for i in superletters]
    thegraphletters=[]
    while theletters:
        nextletter=theletters.pop(0)
        for e in theSG.out_edges(currentvertex, keys=True, data=True):
            if e[3]['superlabel']==nextletter:
                currentvertex=e[1]
                if e[3]['treelabel']:
                    thegraphletters.append(e[3]['treelabel'])
                break
        else:
            currentvertex=None
            break
    if withgraphletters:
        return currentvertex, thegraphletters
    else:
        return currentvertex
   


    
def mark_max_tree(theSG,rootvertex=None, verbose=False):
    """
    In the data of each edge add 'treelabel'. Edges of a maximal subtree get treelabel=0. Other edges get unique labels 1..rank(theSG).
    The specified rootvertex is the root of the maxtree.
    Only marks component containing the rootvertex.
    Any existing treelabels are deleted.
    """
    if rootvertex is None:
        rootvertex=theSG.basepoint
    for e in theSG.edges(keys=True):
        try:
            del theSG[e[0]][e[1]][e[2]]['treelabel']
        except KeyError:
            pass
    unvisited=set(theSG.nodes())
    totalnumverts=len(unvisited)
    if verbose:
        pass
        #fish=ProgressFish(total=totalnumverts) # progress bar
    newverts=[rootvertex]
    unvisited.remove(rootvertex)
    rank=0
    while newverts:
        if verbose:
            pass
            #fish.animate(amount=totalnumverts-len(unvisited))
        edges=theSG.out_edges(newverts,keys=True,data=True)
        newverts=[]
        for edge in edges:
            if 'treelabel' in edge[3]:
                pass # we've already visited this edge or its inverse
            elif edge[1] in unvisited: # this edge is part of the maxtree. label it and its inverse 0
                theSG[edge[0]][edge[1]][edge[2]]['treelabel']=0
                theSG[edge[1]][edge[0]][-edge[2]]['treelabel']=0
                unvisited.remove(edge[1])
                newverts.append(edge[1])
            else: # this edge is not in the maxtree and does not have a label. label it.
                rank+=1
                if edge[2]>0:
                    theSG[edge[0]][edge[1]][edge[2]]['treelabel']=rank
                    theSG[edge[1]][edge[0]][-edge[2]]['treelabel']=-rank
                else:
                    theSG[edge[0]][edge[1]][edge[2]]['treelabel']=-rank
                    theSG[edge[1]][edge[0]][-edge[2]]['treelabel']=rank
    for edge in theSG.edges(keys=True,data=True):# the subtree is complete. any remaining edges are not in it
        if e[0] in unvisited:
            pass # the graph is disconnected and this edge is not in the component containing the basepoint
        elif 'treelabel' in edge[3]:
            pass # we've already visited this edge or its inverse
        else: # this edge is not in the maxtree and does not have a label. label it.
            rank+=1
            if edge[2]>0:
                theSG[edge[0]][edge[1]][edge[2]]['treelabel']=rank
                theSG[edge[1]][edge[0]][-edge[2]]['treelabel']=-rank
            else:
                theSG[edge[0]][edge[1]][edge[2]]['treelabel']=-rank
                theSG[edge[1]][edge[0]][-edge[2]]['treelabel']=rank


def subgroup_intersection(*subgroups, **kwargs):
    """
    Given a finitely generated free group F with finitely generated subgroups, return the intersection as a subgroup of F.
    """
    # Build a Stallings Graph for the intersection from (folded) Stallings Graphs of the subgroups.
    # Vertices are tuples of vertices of the subgroup graphs.
    # basepoint is tuple of basepoints.
    # Edge with superlabel=i from (u_1, u_2,...) to (v_1, v_2,...) if for all j there is an edge in subgroup graph j with superlabel=1 from u_j to v_j.
    # Then take core graph of the basepoint.
    # This is a (folded) Stallings Graph for the intersection.
    try:
        verbose=kwargs['verbose']
    except KeyError:
        verbose=False
    if len(subgroups)<2:
        return subgroups
    F=subgroups[0].supergroup
    assert(all([F is G.supergroup for G in subgroups]))
    basepoint=tuple([G.graph.basepoint for G in subgroups])
    sg=StallingsGraph([],basepoint)
    labellededges=dict()
    for i in range(1,1+F.rank):
        labellededges[i]=[[edge for edge in subgroup.graph.edges(data=True) if edge[2]['superlabel']==i] for subgroup in subgroups]
    if verbose:
        totalnumberedges=sum([prod([len(j) for j in labellededges[i]]) for i in labellededges])
        print str(totalnumberedges)+" edges in the product graph"
        #fish=ProgressFish(total=totalnumberedges)
        numberedges=0
    for i in range(1,1+F.rank):
        for hyperedge in itertools.product(*labellededges[i]):
            origin=tuple([e[0] for e in hyperedge])
            terminus=tuple([e[1] for e in hyperedge])
            sg.add_edge(origin,terminus, superlabel=i)
            if verbose:
                numberedges+=1
                #fish.animate(amount=numberedges)
    if verbose:
        print "Coring the product graph"
    sg.core(verbose=verbose)
    return FGSubgroupOfFree(F,inclusionlist=[],graph=sg)

    
    

def normal_core(G, verbose=False):
    """
    Normal core of G in its supergroup.
    """
    #Take every word in G's max tree and use it as a conjugator. Intersect all the conjugates.
    F=G.supergroup
    assert(G.index()<float('inf'))
    gens=[G.inclusion(G.word([i])) for i in range(1,1+G.rank)] # Basis for G as words in F.
    reps=[F.word(ll) for ll in coset_reps(G.graph)] # Coset representatives for F/G.
    conjugates=[]
    if verbose:
        print "Constructing conjugates"
    for rep in reps:
        conjugategenerators=[rep**(-1)*gen*rep for gen in gens]
        conjugate=FGSubgroupOfFree(F,conjugategenerators)
        conjugates.append(conjugate)
    if verbose:
        print "Intersecting conjugates:"
    core=subgroup_intersection(*conjugates, verbose=verbose)
    return core


def primitive_lift_cover(wordlist, verbose=False):
    """
    A finite index subgroup of w.group in which w lifts to a primitive element.
    """
    # take a labeled loop representing the word and complete it randomly to a covering graph of the F rose.
    F=wordlist[0].group
    sg=StallingsGraph(wordlist)
    sg.fold()
    deficiencies=dict()
    for i in range(1,1+F.rank):
        deficiencies[i]=set(sg.nodes())
        deficiencies[-i]=set(sg.nodes())
    for v in sg.nodes_iter():
        for e in sg.out_edges(v,keys=True,data=True):
            deficiencies[e[3]['superlabel']].remove(v)
    thiskey=1+max([0]+[e[2] for e in sg.edges(keys=True)])
    thistreelabel=1+max([0]+[e[3]['treelabel'] for e in sg.edges(keys=True,data=True)])
    for i in range(1,1+F.rank):
        while(deficiencies[i]):
            u=deficiencies[i].pop()
            v=deficiencies[-i].pop()
            sg.add_edge(u,v,key=thiskey, treelabel=thistreelabel, superlabel=i)
            thiskey+=1
            thistreelabel+=1
    plc=FGSubgroupOfFree(F,inclusionlist=[],graph=sg)
    return plc

def normal_primitive_lift_cover(wordlist, verbose=False):
    return normal_core(primitive_lift_cover(wordlist), verbose=verbose)


def clean_cover(wordlist, verbose=False):
    """
    Given a finitely generated free group F with a list of non-trivial words, return a finitely generated subgroup of F that is a clean cover.
    """
    # intersect normal primitive lift covers of each word in the wordlist.
    # this seems to be faster than building a primitive lift cover for the whole wordlist and then taking normal core
    assert(wordlist)
    F=wordlist[0].group
    assert(all([w.group is F for w in wordlist]))
    wl=[w for w in wordlist]
    nplcs=[]
    if verbose:
        print "Constructing normal primitive lift covers for each word"
    while wl:
        w=wl.pop()
        if w.letters:
            if verbose:
                print "Constructing nplc for word "+w()
            nplc=normal_primitive_lift_cover([w], verbose=verbose)
            nplcs.append(nplc)
    if verbose:
        print "Intersecting the covers"
    cc=subgroup_intersection(*nplcs, verbose=verbose)
    return cc


