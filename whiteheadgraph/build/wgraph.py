import copy
import freegroup
import group
import networkx as nx
import orderedmultigraph as omg
import whiteheadreduce as wreduce
#from fish import ProgressFish



def simplify_wordlist(F, wl1,wl2=[],withmap=False, verbose=False): 
    """
    Find cyclically reduced roots of words in wl1 that are needed to generate distinct conjugacy classes of maximal cyclic subgroups not already generated by roots of elements in wl2.
    """
    progressbar=verbose and (len(wl1)>30)
    if progressbar:
        pass
        #fish=ProgressFish(total=2*len(wl1)+len(wl2))
    powers=[1]*len(wl1)       
    alreadyhave=[F.word([])]*len(wl1)
    for i in range(len(wl2)):
        root, power=F.conjugate_root(wl2[i], withpower=True)
        alreadyhave.append(root)
        powers.append(power)
        if progressbar:
            pass
            #fish.animate(amount=i)
    candidates=[]
    for i in range(len(wl1)):
        root, power=F.conjugate_root(wl1[i], withpower=True)
        candidates.append(root)
        powers[i]=power
        if progressbar:
            pass
            #fish.animate(amount=i+len(wl2))
    need=[]
    selfwordmap=[None]*len(wl1)
    simpwordmap=[None]*len(wl1)
    for wordindex in range(len(wl1)):
        nextword=candidates[wordindex]
        if len(nextword)==0:
            selfwordmap[wordindex]=0
            simpwordmap[wordindex]=0
        else:
            for i in range(len(alreadyhave)):
                if nextword.letters==alreadyhave[i].letters:
                    selfwordmap[wordindex]=i
                    if i<len(wl1): # we found a conjugate power in wl1, note we already know simpwordmap[i]
                        simpwordmap[wordindex]=simpwordmap[i]
                    break
            else: # didn't find nextword in alreadyhave
                need.append(nextword)
                alreadyhave[wordindex]=nextword
                selfwordmap[wordindex]=wordindex
                simpwordmap[wordindex]=len(need)-1
        if progressbar:
            pass
            #fish.animate(amount=wordindex+len(wl1)+len(wl2))
    if withmap:
        return need, {'maptoself':selfwordmap, 'maptosimp':simpwordmap, 'power':powers[:len(wl1)]}
    # for word w=w1[i], 'maptoself'[i]= index of first word to which w is conjugate power in the list w1+w2
    # while 'maptosimp'[i]=index of the image of w in returned list need, which may be None if w is conjugate power of something in wl2
    # if wl2=[] then these are the same
    else:
        return need
freegroup.FGFreeGroup.simplify_wordlist=simplify_wordlist

def blind_simplify_wordlist(F,wordlist):
    """
    Return a minimal list of elements generating maximal cyclic subgroup such that every word in wordlist is conjugate into one of the subgroups.
    The difference between this and simplifyWordlist is that we do not track the correspondence between words in the input and the output.
    """
    Abelianizations=dict()
    for w in wordlist:
        root=F.conjugate_root(w)
        abel=tuple(F.abelianization(root))
        if abel in Abelianizations:
            Abelianizations[abel].add(tuple(root.letters))
        else:
            Abelianizations[abel]=set([tuple(root.letters)])
    needed=[]
    for abel in Abelianizations:
        while Abelianizations[abel]:
            lets=Abelianizations[abel].pop()
            needed.append(F.word(lets))
    return needed


def are_equivalent_wordlists(F, wl1,wl2):
    """
    Decide if given wordlists generate same conjugacy classes of maximal cyclic subgroups in a free group F.
    """
    return not bool(F.simplify_wordlist(wl1,wl2)) and not bool(F.simplify_wordlist(wl2,wl1))
freegroup.FGFreeGroup.are_equivalent_wordlists=are_equivalent_wordlists

class WGraph(omg.OrderedMultiGraph):
    """ 
    Thw Whitehead Graph is an omg.OrderedMultiGraph constructed from a list of group.word in a free group F together with a dictionary of splicing maps splicemaps={vert:splicemap for splicing vert to -vert,...}
    """
    # set simplified=True to skip simplification step
    # set autominimize=True to automatically reduce to minimal complexity
    def __init__(self, originalwordlist, autominimize=False, simplified=False, verbose=False):
        omg.OrderedMultiGraph.__init__(self)
        self.group=originalwordlist[0].group
        self.rank=self.group.rank
        self.splicemaps={}
        if simplified:
            simplifiedwordlist=originalwordlist
        else:
            if verbose:
                print "Simplifying wordlist"
            simplifiedwordlist=simplify_wordlist(self.group,originalwordlist)
        if autominimize:
            if verbose:
                print "Performing peak reduction"
            self.wordlist=wreduce.whitehead_minimal(self.group,simplifiedwordlist,simplified=True,verbose=verbose,blind=True)['wordlist']
        else:
            self.wordlist=simplifiedwordlist
                    
        for i in set(range(-self.rank,self.rank+1))-set([0]):
            self.add_vertex(i)
            self.splicemaps[i]=[]
            
        edgecounter=0
        if verbose:
            print "Constructing Whitehead Graph"
            manywords=bool(len(self.wordlist)>999)
            if manywords:
                pass
                #wgfish=ProgressFish(total=len(self.wordlist))
        for i in range(len(self.wordlist)):
            w=self.wordlist[i]
            if verbose and manywords:
                pass
                #wgfish.animate(amount=i)
            if len(w)>0:
                firstplace=len(self.splicemaps[-w.letters[0]])
                self.splicemaps[-w.letters[0]]+=[None]
                if len(w)>1:
                    for i in range(1,len(w)):
                        self.add_edge(-w.letters[i-1],w.letters[i],'e'+str(edgecounter))
                        edgecounter+=1
                        self.splicemaps[w.letters[i]]+=[len(self.splicemaps[-w.letters[i]])]
                        self.splicemaps[-w.letters[i]]+=[len(self.splicemaps[w.letters[i]])-1]
                self.add_edge(-w.letters[-1],w.letters[0],'e'+str(edgecounter))
                edgecounter+=1
                self.splicemaps[w.letters[0]]+=[firstplace]
                self.splicemaps[-w.letters[0]][firstplace]=len(self.splicemaps[w.letters[0]])-1

    def inv(self,vert):
        return -vert

    def next_edge(self,edge,vert):
        return self.edge_order(self.inv(vert))[self.splicemaps[vert][self.edge_order(vert).index(edge)]]
            
    def get_wordlist(self):
        """
        return a wordlist generating this Whitehead graph
        """
        if len(self.edges())==0:
            return [self.group.word([])]
        traversed=set()
        wordlist=[]
        while traversed!=set(self.edgekeys):
            firstedge=(set(self.edgekeys)-traversed).pop()
            firstvert=self.origin(firstedge)
            nextvert=self.opposite_end(firstedge,firstvert)
            currentword=[nextvert]
            traversed.add(firstedge)
            nextedge=self.next_edge(firstedge,nextvert) 
            while nextedge!=firstedge:
                nextvert=self.opposite_end(nextedge,-nextvert)
                currentword.append(nextvert)
                traversed.add(nextedge)
                nextedge=self.next_edge(nextedge, nextvert)
            wordlist.append(self.group.word(currentword))
        return wordlist

    def __copy__(self):
        return WGraph(self.wordlist)

    def permute_edge_order_at_vertex(self,thisvertex,inputpermutation,lastfixed=False):
        """
        Change the ordering of incident edges at thisvertex by applying thispermutation.
        Update splicemap accordingly.
        If zerofixed=True then thispermution is a permutation of the first valence-1 edges, with that last being fixed.
        """
        valence=self.valence(thisvertex)
        if lastfixed:
            thispermutation=inputpermutation+(valence-1,)
        else:
            thispermutation=inputpermutation
        newedgeorder=dict([(thispermutation[i],self.node[thisvertex]['edgeorder'][i]) for i in range(valence)])
        newsplicemap=dict([(thispermutation[i],self.splicemaps[thisvertex][i]) for i in range(valence)])
        newinversesplicemap=dict([(j,thispermutation[self.splicemaps[-thisvertex][j]]) for j in range(valence)])
        self.node[thisvertex]['edgeorder']=[newedgeorder[i] for i in range(valence)]
        self.splicemaps[thisvertex]=[newsplicemap[i] for i in range(valence)]
        self.splicemaps[-thisvertex]=[newinversesplicemap[i] for i in range(valence)]

    def make_edge_orders_consistent(self):
        """
        For each basis element of F change the edge order at the inverse vertex to make the splicing maps reverse edge orders.
        """
        for v in range(1,self.rank+1):
            thisvertexpermutation=[(-1-(self.splicemaps[v].index(i)))%self.valence(v) for i in range(self.valence(v))]
            self.permute_edge_order_at_vertex(-v,thisvertexpermutation)

        
        

def wgrow_one(G,W,vert):
    """
    Extend generalized Whitehead graph G by splicing on Whitehead graph W at vertex vert.
    """
    assert(vert in G)
    assert(type(vert)==tuple)
    prefix = vert
    Gsplicevert= vert
    Wsplicevert= -vert[-1]
    return omg.splice(G,W,Gsplicevert,Wsplicevert,W.splicemaps[vert[-1]],(),prefix)

def wgrow_word(W,w):
    """
    Returns generalized Whitehead graph over the segment w
    """
    Gedges=[((edge[0],), (edge[1],), (edge[2],), edge[3]) for edge in W.edges_iter(keys=True, data=True)]
    Gedgeorders=dict(((vert,),[(edge,) for edge in W.node[vert]['edgeorder']]) for vert in W.nodes_iter())  #{(vert,):[(edge,) for edge in W.node[vert]['edgeorder']] for vert in W.nodes_iter()}
    G=omg.OrderedMultiGraph(Gedges,Gedgeorders)
    for vert in W.nodes_iter(): # We have missed isolated vertices. Add them back.
        if W.valence(vert)==0:
            G.add_vertex((vert,))
    # G is a deep copy of W so that vertex and edge names are all length 1 tuples
    if len(w)>0:
        prefix=()
        for i in range(len(w)):
            prefix+=(w.letters[i],)
            G=wgrow_one(G,W,prefix)
    return G






def simplify_and_minimize_wordlist(F,wordlist, simplified=False, minimized=False, extrawordlist=None,blind=False,cutvertsonly=False, verbose=False):
    """
    Simplify and Whiteheadminimize the wordlist wl, and return newwordlist, wordmap, sequence of minimizing Whitehead automorphisms so that minimizingautomorphism(wl[i]) is conjugate to newwordlist[j]**p. where wordmap[i]=(j,p)
    """
    if not simplified:
        simplifiedwordlist, simpmap = simplify_wordlist(F,wordlist,withmap=True)
        wordmap=[]
        for i in range(len(wordlist)):
            wordmap.append((simpmap['maptosimp'][i], simpmap['power'][i]))
    else:
        simplifiedwordlist=wordlist
        wordmap=[(i,1) for i in range(len(simplifiedwordlist))]
    wm=dict([('wordmap',wordmap),('connected',None),('wordlist',simplifiedwordlist),('minimizingautomorphism',None),('inverseminimizer',None), ('whiteheadsequence',None),('extrawordlist',extrawordlist)])
    if minimized:
        if not blind:
            wm['minimizingautomomorphism']=group.Automorphism(F)
            wm['inverseminimizer']=group.Automorphism(F)
            wm['whiteheadsequence']=[]
    else:
        wm.update(wreduce.whitehead_minimal(F,simplifiedwordlist, extrawordlist=extrawordlist,simplified=True, blind=blind, cutvertsonly=cutvertsonly, verbose=verbose))    
    return wm




def wgparse(F, wgraphorwordlist, extrawordlist=None, simplifyandminimize=False, simplified=False, minimized=False, blind=False, cutvertsonly=False, verbose=False):
    """
    Take as input either a whiteheadgraph or wordlist and return a dict whose entries include 'WhiteheadGraph' and 'wordlist'
    Note that if a wgraphorwordlist is a wordlist then the wordlist and Whitehead graph returned will be in terms of the group F.

    simplified=True means the imput wordlist is already simplified, do not spend time resimplifying

    minimized=True means the input wordlist is already minimized, do not spend time looking for minimizations

    simplifyandminimize=True means simplify and minimize the wordlist (or just minimize, if simplified=True)

    cutvertsonly=True means only reduce until the resulting Whitehead graph has no cut vertices, not necessarily all the way to minimal complexity

    blind=True means do not track minimizing autmorphism

    extrawordlist is a list of words in F. If simplifyandminimize=True the reducing automorphisms will also be applied to them.
    """
    results=dict([('wordlist',None),('minimizingautomorphism',None),('inverseminimizer',None),('wordmap',None),('extrawordlist',extrawordlist)])
    try:
        results['originalwordlist']=wgraphorwordlist.wordlist
        results['wordlist']=results['originalwordlist']
        results['WhiteheadGraph']=wgraphorwordlist
    except AttributeError: # its not a whiteheadgraph. had better be a wordlist
        try:
            inc=wgraphorwordlist[0].group.get_inclusion(F)
        except KeyError:
            results['originalwordlist']=[F.word([])] # empty wordlist
        else:
            results['originalwordlist']=[inc(w) for w in wgraphorwordlist] # translates all the words to be in the common supergroup F
        if (not simplifyandminimize) or (simplified and minimized):
            results['wordlist']=results['originalwordlist']
            results['WhiteheadGraph']=WGraph(results['wordlist'], simplified=simplified)
    if simplifyandminimize:
        smwl=simplify_and_minimize_wordlist(F,results['originalwordlist'], simplified=simplified, minimized=minimized,extrawordlist=extrawordlist, blind=blind,cutvertsonly=cutvertsonly, verbose=verbose)
        results.update(smwl)
        if not (simplified and minimized):
            results['WhiteheadGraph']=WGraph(results['wordlist'], simplified=True,verbose=verbose)
    return results        