import copy
import group
import freegroup
import whiteheadgraph.build.wgraph as wg
import whiteheadgraph.build.orderedmultigraph as omg
import networkx as nx
import partition as part
import whiteheadgraph.build.whiteheadreduce as wreduce
import AutF
import graphofgroups as gog
#from fish import ProgressFish
import enumeratewords




class TooBigError(Exception):
     def __init__(self, value,howbig=None):
         self.value = value
         self.howbig=howbig
     def __str__(self):
         return repr(self.value)

def is_primitive(F,w):
    """
    Decide if w is a primitive word in F.
    """
    r,p=F.max_root(w,uptoconjugacy=True,withpower=True)
    if p != 1:
        return False
    W=wg.WGraph([r],autominimize=True)
    if len(W.edges())==1:
        return True
    else:
        return False
freegroup.FGFreeGroup.is_primitive=is_primitive

def is_subbasic(F,wordlist):
    """
    Decide if the wordlist determines the same list of conjugacy classes of maximal cyclic subgroups as a basis.
    """
    simplifiedwordlist=wg.blind_simplify_wordlist(F,wordlist)
    W=wg.WGraph(simplifiedwordlist,autominimize=True)
    return len(simplifiedwordlist)==len(W.edges())
freegroup.FGFreeGroup.is_subbasic=is_subbasic

def splits_freely_rel(F, wordlist, simplified=False, minimized=False, verbose=False):
    """
    Decide if F  splits freely relative to the multiword.
    """
    return not wreduce.whitehead_minimal(F,wordlist, blind=True,stopatdisconnected=True,simplified=simplified,verbose=verbose,cutvertsonly=True)['connected']
freegroup.FGFreeGroup.splits_freely_rel=splits_freely_rel

def is_circle(F, whiteheadgraphorwordlist, simplified=False, minimized=False, verbose=False):
    """
    Decide if Whitehead graph is a circle.
    """
    W=wg.wgparse(F,whiteheadgraphorwordlist, simplified=simplified, minimized=minimized,verbose=verbose,simplifyandminimize=True,blind=True)['WhiteheadGraph']
    return W.is_circle()
freegroup.FGFreeGroup.is_circle=is_circle
freegroup.FGFreeGroup.is_QH=is_circle


def get_free_splitting_rel(F, originalwordlist, simplified=False, minimized=False, verbose=False, withwordmap=False, blind=False, printresult=False):
    """
    Find canonical maximal free splitting of free group of given rank relative to words.

    Returns graphofgroups with one vertex v0 with trivial stabilizer, other vertex stabilizers non-trivial, and all edge stabilizers trivial. Non-trivial vertex stabilizers do not split further.

    withwordmap=True then also return a list whose i-th entry is a tuple (vertex, word w in vertex stabilizer, power) such that word i of the original wordlist is conjugate to w**(+-power) in stabilizer of vertex.

    blind=False then return splitting of F relative to the original wordlist, ie, in the original basis
    blind=True then return a splitting of F relative to the minimized wordlist

    Use simplified=True if wordlist is already simplified to avoid spending time resimplifying
    Use minimized=True if wordlist is already Whitehead reduced at least enough so that the Whitehead graph is cut vertex free and connected components are inverse closed 
    """
    wgp=wg.wgparse(F,originalwordlist, simplifyandminimize=True, simplified=simplified, minimized=minimized,verbose=verbose, cutvertsonly=True)
    W=wgp['WhiteheadGraph']
    wordlist=wgp['wordlist']
    wordmap=wgp['wordmap']
    wordgens=[set([abs(i) for i in w.letters]) for w in wordlist] # the generators of the new basis that are used in the wordlist
    usedpartition=part.make_partition(*wordgens) # partition of the used generators according to whether they are used in a common word
    partcontainingword=[]
    for i in range(len(wordlist)):
        for j in range(len(usedpartition.parts)): 
            if abs(wordlist[i].letters[0]) in usedpartition.parts[j]:
                partcontainingword.append(j) # says word wordlist[i] is in part j
                break
        else:
            raise IndexError('There is something wrong. We should always find the generator in one of the parts.')
            
    unused=set(range(1,1+F.rank))-usedpartition.elements()
    # thepartition is a list of disjoint subsets of generators of the free group, but not in the original basis
    # to get the answer in terms of the original basis need to apply inverse sequence of whitehead automorphisms
    if not blind:
        obused=[set([wgp['inverseminimizer'](F.word([i])) for i in p]) for p in usedpartition.parts]
        obunused=set([wgp['inverseminimizer'](F.word([i])) for i in unused])
    else:
        obused=[set([F.word([i]) for i in p]) for p in usedpartition.parts]
        obunused=set([F.word([i]) for i in unused])
    if verbose:
        print "Constructing graph of groups."
    Gamma=gog.FPGraphOfGroups()
    if len(usedpartition.parts)>1:
        Gamma.add_vertex('v0', freegroup.FGSubgroupOfFree(F,[]))
        for i in range(len(obused)):
            Gamma.add_vertex('v'+str(i+1), freegroup.FGSubgroupOfFree(F,[w for w in obused[i]]))
            Gamma.add_edge('v0','v'+str(i+1), freegroup.FGSubgroupOfFree(F,[]),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v0']['group']),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v'+str(i+1)]['group']))
        for w in obunused:
            Gamma.add_edge('v0','v0',freegroup.FGSubgroupOfFree(F,[]),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v0']['group']),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v0']['group']), label=w() )
    elif (len(usedpartition.parts)==1):
        Gamma.add_vertex('v1', freegroup.FGSubgroupOfFree(F,[w for w in obused[0]]))
        for w in obunused:
            Gamma.add_edge('v1','v1',freegroup.FGSubgroupOfFree(F,[]),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v1']['group']),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v1']['group']), label=w() )
    else:
         Gamma.add_vertex('v0', freegroup.FGSubgroupOfFree(F,[]))
         for w in obunused:
            Gamma.add_edge('v0','v0',freegroup.FGSubgroupOfFree(F,[]),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v0']['group']),group.Homomorphism(freegroup.FGSubgroupOfFree(F,[]), Gamma.node['v0']['group']), label=w() )
    if not withwordmap:
        if printresult:
            print Gamma
        return Gamma
    else:
        if verbose:
            print "Finding conjugates of words in the vertex groups."
            #fishy=ProgressFish(total=len(wordmap))
        for i in range(len(wordmap)):
            if verbose:
                pass
                #fishy.animate(amount=i)
            whichvert='v'+str(1+partcontainingword[wordmap[i][0]])
            whichgroup=Gamma.localgroup(whichvert)
            if not blind:
                whichword=whichgroup.find_conjugate_in_subgroup(wgp['inverseminimizer'](wordlist[wordmap[i][0]]))
            else:
                whichword=whichgroup.find_conjugate_in_subgroup(wordlist[wordmap[i][0]])
            wordmap[i]=(whichvert, whichword, wordmap[i][1])
        if printresult:
            print Gamma
            if blind:
                print "\n".join(["Image of '"+str(originalwordlist[i]())+"' is conjugate to ('"+str(wordmap[i][1]())+"')**('"+str(wordmap[i][2])+"') in vertex "+wordmap[i][0] for i in range(len(originalwordlist))])
            else:
                print "\n".join(["'"+str(originalwordlist[i]())+"' is conjugate to ('"+str(wordmap[i][1]())+"')**('"+str(wordmap[i][2])+"') in vertex "+wordmap[i][0] for i in range(len(originalwordlist))])
        return Gamma, wordmap
freegroup.FGFreeGroup.get_free_splitting_rel=get_free_splitting_rel

def gives_cut(F,whiteheadgraphorwordlist,inputw,returnnumbercomponents=False,simplified=False, minimized=False, verbose=False):
    """
    Check if endpoints of word w give cut point/pair in decomposition space for Whitehead graph W.
    """
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, extrawordlist=[inputw], simplified=simplified, minimized=minimized,verbose=verbose,simplifyandminimize=True, blind=True)
    W=wgp['WhiteheadGraph']
    wordlist=wgp['wordlist']
    w=wgp['extrawordlist'][0]
    if len(w)==0:
        return False
    else:
        prefixw=F.word(w.letters[0:len(w)-1])
        G=wg.wgrow_word(W,prefixw)
        v1=tuple(w.letters)
        v2=(-w.letters[-1],)
        # G is a generalized whitehead graph. The w action will identify vertices v1 and v2.
        # G-{v1,v2} is a findamental domain for the action of w on the infinte whitehead graph over <w>.
        # We can compute the number of complementary components of this infinte whitehead graph by understanding how the splicemap interacts with connected components of the finite whitehead graph.
        
        components=G.connected_components_minus_two_vertices(v2,v1)
        edgelist1=[[e for e in range(0, G.valence(v1)) if G.opposite_end(G.edge_order(v1)[e],v1) in components[i]] for i in range(0,len(components))] # edgelist1[k] is list of loose edges at v1 that connected to vertices in component k
        edgelist2=[[e for e in range(0, G.valence(v2)) if G.opposite_end(G.edge_order(v2)[e],v2) in components[i]] for i in range(0,len(components))]      # edgelist2[k] is list of loose edges at v2 that connected to vertices in component k             
        missededges=set(G.incident_edges(v1))-set.union(*[set([G.edge_order(v1)[e] for e in p]) for p in edgelist1]) # edges that go directly between v1 and v2
        while missededges!=set([]):
            nextedge=missededges.pop()
            edgelist1+=[[G.edge_order(v1).index(nextedge)]]
            edgelist2+=[[G.edge_order(v2).index(nextedge)]]
        P1=part.Partition(edgelist1)
        P2=part.Partition(edgelist2)
        partitionmap=range(0,len(P1.parts))
        splicemap=W.splicemaps[w.letters[-1]]
        (newP1,newP2)=part.compatible_coarsenings(P1,P2,partitionmap,splicemap)
        if F.is_conjugate_into(w,*W.wordlist): # If w is conjugate into the wordlist then one part of the partition does not correspond to a component, just a segregated edge.
            numberofcomponents=len(newP1.parts)-1
        else:
            numberofcomponents=len(newP1.parts)
        if returnnumbercomponents:
            return numberofcomponents
        else:
            return bool(numberofcomponents -1)
        
        

def number_complementary_components(F, whiteheadgraphorwordlist,w,simplified=False, minimized=False, verbose=False):
    """
    Find numberof complementary components of <w> in W
    """
    return gives_cut(F,whiteheadgraphorwordlist,w,returnnumbercomponents=True,simplified=simplified, minimized=minimized,verbose=verbose)

def find_cut_points(F, whiteheadgraphorwordlist, simplified=False, minimized=False, verbose=False):
    """
    Return list of reperesentatives of the cut points.
    """
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, simplified=simplified, minimized=minimized,verbose=verbose, blind=True, simplifyandminimize=True)
    W=wgp['WhiteheadGraph']
    wordlist=wgp['wordlist']
    wordmap=wgp['wordmap']
    # every cut point is stabilized by a conjugate of one of the generating words, so just check if they give cuts.
    cutpoints=[]
    manywords=bool(len(wordlist)>9)
    if verbose and manywords:
        #fish=ProgressFish(total=len(wordlist))
        print "Checking if generating words give cut points."
    for j in range(len(wordlist)):
        if verbose and manywords:
            pass
            #fish.animate(amount=j)
        w=wordlist[j]
        if gives_cut(F,W,w, simplified=True,minimized=True,verbose=verbose):
            for i in range(len(wgp['originalwordlist'])):
                if wordmap[i][0]==j:
                    if verbose:
                        print "The word at index "+str(i)+" gives a cut point."
                    break
            else:
                raise KeyError
            cutpoints.append(F.max_root(wgp['originalwordlist'][i], uptoconjugacy=True)[0])
    return cutpoints

def crossing_cut_pairs(F, whiteheadgraphorwordlist,w1,w2, simplified=False, minimized=False, verbose=False, theyareknowncutpairs=False):
    """
    Decide if w1 and w1 give crossing cut pairs for W.

    set theyareknowncutpairs=True if w1 and w2 are already known to give cup pairs.
    """
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, simplified=simplified, minimized=minimized,verbose=verbose, blind=True)
    W=wgp['WhiteheadGraph']
    w1=F.cyclic_reduce(w1)
    w2=F.cyclic_reduce(w2)
    if theyareknowncutpairs: # skip cut pair check
        return not gives_cut(F,W.wordlist+[w1],w2) # If they cross and we add one of them to the wordlist then the second word no longer gives a cut.
    else:    # First check if they are even cut pairs.
        if not number_complementary_components(F,W,w1)==2 and number_complementary_components(F,W,w2)==2: 
            return False
        else:
            return not gives_cut(F,W.wordlist+[w1],w2) # If they cross and we add one of them to the wordlist then the second word no longer gives a cut.

def gives_splitting(F, whiteheadgraphorwordlist,w,simplified=False, minimized=False, verbose=False):
    """
    Decide if W splits over <w>.
    """
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, extrawordlist=[w], simplified=simplified, minimized=minimized,verbose=verbose, blind=True)
    return gives_cut(F,wgp['WhiteheadGraph'],wgp['extrawordlist'][0]) and not crossing_cut_pairs(F,wgp['WhiteheadGraph'],wgp['extrawordlist'][0],wgp['extrawordlist'][0]) # <w> gives a splitting if it gives a cut and it doesn't cross itself.


    
def push_forward_partition(W,v0,P0,v1,precomputedcomponents=None):
    """
    Find partitions newP0 of of v0 edges and P1 of v1 edges compatible with P0 and connectivity in W-{v0,v1}
    """
    assert(v0!=v1)
    # This function gets called a lot. Let's allow it to remember the components computed in the next step.
    if precomputedcomponents:
        if (v0,v1) in precomputedcomponents:
            components=precomputedcomponents[(v0,v1)]
        else:
            components=W.connected_components_minus_two_vertices(v0,v1)
            precomputedcomponents[(v0,v1)]=components
    else:
        components=W.connected_components_minus_two_vertices(v0,v1)
    newgraph=nx.Graph()
    for c in components:
        newgraph.add_star([(n,'vert') for n in c])
    for p in P0.parts:
        newgraph.add_star([(e,'edge0') for e in p])
    for i in range(W.valence(v1)):
        if W.opposite_end(W.incident_edges(v1)[i],v1)==v0:
            newgraph.add_edge((i,'edge1'),(W.incident_edges(v0).index(W.incident_edges(v1)[i]),'edge0'))
        else:
            newgraph.add_edge((i,'edge1'),(W.opposite_end(W.incident_edges(v1)[i],v1),'vert'))
    for i in range(W.valence(v0)):
        if W.opposite_end(W.incident_edges(v0)[i],v0)!=v1:
            newgraph.add_edge((i,'edge0'),(W.opposite_end(W.incident_edges(v0)[i],v0),'vert'))
    newcomponents=nx.connected_components(newgraph)
    edgelist0=[]
    edgelist1=[]
    for newcomp in newcomponents:
        part0=[]
        part1=[]
        for n in newcomp:
            if n[1]=='edge0':
                part0+=[n[0]]
            if n[1]=='edge1':
                part1+=[n[0]]
        edgelist0+=[part0]
        edgelist1+=[part1]
    #for i in range(len(newcomponents)):
        #part0=[]
        #part1=[]
        #for n in newcomponents[i]:
        #    if n[1]=='edge0':
        #        part0+=[n[0]]
        #    if n[1]=='edge1':
        #        part1+=[n[0]]
        #edgelist0+=[part0]
        #edgelist1+=[part1]
    newP0=part.Partition(edgelist0)
    part0coarseningmap=[]
    for i in range(len(P0.parts)):
        for e in P0.parts[i]: break
        part0coarseningmap+=[newP0.which_part(e)]
    return (part0coarseningmap, newP0 , part.Partition(edgelist1))
            


    
def find_universal_splitting_words(F, W, wordlist, DoNotVerifyTwoComponentWords=False, StopAtFirstCut=False, MinNumComponents=2, simplified=False, minimized=False, verbose=False, check3LetterSubwords=True,cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None):
    #repalces findCutPairs
    """
    Find cut pairs for a whitehead graph.
    """
    if check3LetterSubwords:
        if contains_all_3_letter_subwords(*wordlist): # if the multiword contains all 3 letter subwords of F then it is rigid.
            return ({'cutpoints':set([]),'uncrossed':set([]),'othercuts':set([])},True)
    rank=F.rank
    precomputedcomponents=dict() # this will be a dict with key (v0,v1) containing the components of W-{v0,v1}, populated as needed
    whiteheadgraphiscomplete=False
    simplegraph=nx.Graph(W)
    if len(simplegraph.edges())==rank*(2*rank-1):
        whiteheadgraphiscomplete=True
    del simplegraph
    directions=range(-rank,rank+1)
    directions.remove(0)
    buds=set({})
    SM=nx.DiGraph()
    # SM is a finite state machine.
    # Nodes correspond to free group generator with partition of edges of Whitehead graph crossing that edge in cayley tree.
    # Edges correspond to a legal turn in the free group, together with a coarsening map of partitions coming from connectivity in the Whitehead graph.
    # If we find a loop in the state machine so that the partitions have more than one part then we have found a cut word.
    # Buds are newly added nodes that need to have outgoing edges computed.
    for urvert in range(-rank,0):
        urpart=part.Partition([[i] for i in range(W.valence(urvert))])
        SM.add_node((urvert, urpart))
        buds.add((urvert, urpart))
    def extendSM(W,SM,buds,directions, maxlength=None, MinNumComponents=2):
        """
        Recursively extend finite state machine holding edge partitions.
        """
        # buds are the newly added states that we still need to compute outgoing edges for
        # newbuds will be the buds in the next iteration
        # maxlength is bound on number of steps to extend the state machine from the original states
        newbuds=set({})
        while buds:
            thisbud=buds.pop()
            indirec=thisbud[0]
            inpart=thisbud[1]
            for outdirec in directions:
                if outdirec!=indirec: # don't backtrack
                    (coarseningmap,coarsenedinpart,outpart)=push_forward_partition(W,indirec,inpart,outdirec,precomputedcomponents)
                    keepgoing=True
                    if whiteheadgraphiscomplete and len(outpart.parts)==2:
                            if len(outpart.parts[0])==1 or len(outpart.parts[1])==1:
                                keepgoing=False
                    if len(outpart.parts)<MinNumComponents:
                        keepgoing=False
                    if keepgoing: # only keep going if we potentially have at least  minnumbercomponents components in this direction
                        newindirec=-outdirec
                        newinpartslist=[[] for i in range(len(outpart.parts))]
                        for i in range(W.valence(newindirec)):
                            newinpartslist[outpart.which_part(W.splicemaps[newindirec][i])]+=[i]
                        newinpart=part.Partition(newinpartslist)
                        newbud=(newindirec,newinpart)
                        # check if newbud is already in SM
                        for n in SM:
                            if newbud[0]==n[0]:
                                if part.is_reordered_partition(newbud[1],n[1]): # newbud is already in SM, so don't add a new vertex
                                    SM.add_edge(thisbud,n)#,{'label':outdirec,'coarseningmap':newcoarseningmap})
                                    break
                        else: # newbud is a new node.
                            SM.add_node(newbud)
                            SM.add_edge(thisbud,newbud)#,{'label':outdirec,'coarseningmap':coarseningmap})
                            newbuds.add(newbud)
        # we have now extended all the original buds, but if we created new vertices we need to recurse.
        buds.update(newbuds)
        if buds:
            if maxlength==None:
                extendSM(W,SM,buds,directions,maxlength, MinNumComponents)
            elif maxlength>0:
                extendSM(W,SM,buds,directions,maxlength-1, MinNumComponents)
  
    extendSM(W,SM,buds,directions,cutpairsearchrecursionlimit, MinNumComponents)
    if buds:
        surethatsall=False # If there are still buds that means we reached maxlength and didn't finish building the state machine. It may contain some cycles, so that we still get cut pairs. If not, it may be because we didn't let the state machine build far enough out.
    else:
        surethatsall=True
    cycles=nx.simple_cycles(SM) # get the simple cycles in the state machine.
    cutpoints=set([])
    uncrossed=set([])
    othercuts=set([])
    foundacut=False
    for cycle in cycles: 
        thewordletters=[]
        if len(cycle)==1 or cycle[0]!=cycle[-1]: # in networkx >=1.8 cycle is a list of nodes with no repetition
            for i in range(len(cycle)): # read off the word of the free group from the cycle in the state machine
                thewordletters+=[-cycle[i][0]]
        else: # networkx <=1.7 cycle is a list of vertices with first and last nodes equal
            for i in range(1,len(cycle)): # read off the word of the free group from the cycle in the state machine
                thewordletters+=[-cycle[i][0]]
        theword=F.conjugate_root(F.word(thewordletters))
        wordinlist=bool(F.is_conjugate_into(theword,*wordlist)) # see if theword is in the generating wordlist
        if wordinlist:
            complementarycomponents=len(cycle[0][1].parts)-1 # if theword is in the wordlist then one component is just the word itself and not a complementary component
        else:
            complementarycomponents=len(cycle[0][1].parts)
        if wordinlist and complementarycomponents>1:
            cutpoints.add(tuple(theword.letters))
            foundacut=True
        elif complementarycomponents>2:
            uncrossed.add(tuple(theword.letters))
            foundacut=True
        elif complementarycomponents>1:
            othercuts.add(tuple(theword.letters))
            foundacut=True
        if foundacut and StopAtFirstCut:
                return ({'cutpoints':set([F.word(w) for w in cutpoints]),'uncrossed':set([F.word(w) for w in uncrossed]),'othercuts':set([F.word(w) for w in othercuts])},surethatsall)
                        
    potentiallyuncrossed=list(othercuts-uncrossed)
    if verbose:
        print "Found "+str(len(cutpoints))+" cut points, "+str(len(uncrossed))+" uncrossed cut pairs, and "+str(len(othercuts))+" other potential cuts."
    if maxnumberof2componentcutstoconsider:
        if len(othercuts)>maxnumberof2componentcutstoconsider:
            raise TooBigError(str(len(othercuts))+" potential cut pairs is beyond limit set by 'maxnumberof2componentcutstoconsider'",len(othercuts))

    reducedcutpoints=set([F.word(t) for t in cutpoints])
    if DoNotVerifyTwoComponentWords:
        reduceduncrossed=set([F.word(t) for t in uncrossed])
        reducedothercuts=set([F.word(t) for t in potentiallyuncrossed])
    else:
        uncrossed|=verify_uncrossed_splitting_words(F,wordlist,potentiallyuncrossed,verbose=verbose)
        reduceduncrossed=set([F.word(t) for t in uncrossed])
        reducedothercuts=set([F.word(t) for t in othercuts-uncrossed])
    if verbose:
        print "Found "+str(len(reducedcutpoints)+len(reduceduncrossed))+" splitting elements."
    return ({'cutpoints':reducedcutpoints,'uncrossed':reduceduncrossed,'othercuts':reducedothercuts},surethatsall)

def verify_uncrossed_splitting_words(F,wordlist,potentiallyuncrossed, verbose=False):
    """
    Take a list of letter tuples that are known to give cup pairs and check if they cross each other. Return set containing those that are uncrossed.
    """
    uncrossed=set([])
    # check if any of the othercuts are really 2 component uncrossed cut pairs.
    if verbose and len(potentiallyuncrossed)>20:
        pass
        #fish= ProgressFish(total=len(potentiallyuncrossed))
    for i in range(len(potentiallyuncrossed)):
        thisword=potentiallyuncrossed[i]
        if not crossing_cut_pairs(F,wordlist, F.word(thisword),F.word(thisword), theyareknowncutpairs=True): # If it crosses itself we can not split over thisword.
            for otherword in potentiallyuncrossed: # If it doesn't cross itself it does give a splitting, but maybe in a surface component if thisword crosses something else.
                if crossing_cut_pairs(F,wordlist,F.word(thisword),F.word(otherword),theyareknowncutpairs=True):
                    break
            else:
                uncrossed.add(thisword)
        if verbose and len(othercuts)>20:
            pass
            #fish.animate(amount=i+1)
    return uncrossed



def is_rigid_rel(F, whiteheadgraphorwordlist, simplified=False, minimized=False, verbose=False,cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None):
    """
    Decide if W is rigid.
    """
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, blind=True, simplified=simplified, minimized=minimized,verbose=verbose, simplifyandminimize=True)
    W=wgp['WhiteheadGraph']
    wordlist=wgp['wordlist']
    if not wgp['connected']:
        return False
    elif W.is_circle():
        return False
    else:
        try:
            cuts,surethatsall=find_universal_splitting_words(F,W,wordlist, DoNotVerifyTwoComponentWords=True, StopAtFirstCut=True, simplified=True, minimized=True, verbose=verbose,cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider)
        except TooBigError:
            raise TooBigError("Could not determine with 'cutpairsearchrecursionlimit'="+str(cutpairsearchrecursionlimit))
        if set.union(cuts['cutpoints'],cuts['uncrossed'], cuts['othercuts']): # found some cuts
            return False
        else: # we didn't find any cuts
            if surethatsall: #  the state machine was complete, so we're sure there really are no cuts
                return True
            else: # otherwise we quit looking because the state machine got too big, so we can't be sure that there are really no cuts
                raise TooBigError("Could not determine with given value of 'maxnumberof2componentcutstoconsider'="+str(maxnumberof2componentcutstoconsider))
freegroup.FGFreeGroup.is_rigid_rel=is_rigid_rel
    
def smash(prefix,oldname):
    try:
        return prefix+oldname
    except TypeError:
        try:
            return prefix+(oldname,)
        except TypeError:
            try:
                return (prefix,)+oldname
            except TypeError:
                return (prefix,)+(oldname,)

def get_relative_cyclic_splitting_over(F, W, wordlist, splittingword, nameprefix='', extrawordlist=None, simplified=False, minimized=False, verbose=False,cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None):
    """
    Return graph of groups splitting of F relative to a multiword over the cyclic subgroup < splittingword >, and 
    a list whose ith element is (imagevertex,imageword in imagevertex group) of the ith element of the original wordlist.
    Extrawordlist is a list of words in F whose images we also want to track. These words must be known to be elliptic in the splitting. If extrawordlist is present then a second wordlist similar to above is also returned.
    """
    # names of vertices prefixed by nameprefix
    # splitting word should already be a nontrivial cyclically reduced and not a proper power, wordlist should be simplified and Whitehead minimized. Should be no free splittings.
    w=splittingword
    #--------- from givesCut
    prefixw=F.word(w.letters[0:len(w)-1])
    G=wg.wgrow_word(W,prefixw)
    v1=tuple(w.letters)
    v2=(-w.letters[-1],)
    # G is a generalized whitehead graph. The w action will identify vertices v1 and v2.
    # G-{v1,v2} is a fundamental domain for the action of w on the infinte whitehead graph over <w>.
    # We can compute the number of complementary components of this infinte whitehead graph by understanding how the splicemap interacts with connected components of the finite whitehead graph.
        
    components=G.connected_components_minus_two_vertices(v2,v1)
    edgelist1=[[e for e in range(0, G.valence(v1)) if G.opposite_end(G.edge_order(v1)[e],v1) in components[i]] for i in range(0,len(components))] # edgelist1[k] is list of loose edges at v1 that connected to vertices in component k
    edgelist2=[[e for e in range(0, G.valence(v2)) if G.opposite_end(G.edge_order(v2)[e],v2) in components[i]] for i in range(0,len(components))]      # edgelist2[k] is list of loose edges at v2 that connected to vertices in component k             
    missededges=set(G.incident_edges(v1))-set.union(*[set([G.edge_order(v1)[e] for e in p]) for p in edgelist1]) # edges that go directly between v1 and v2
    while missededges!=set([]):
        nextedge=missededges.pop()
        edgelist1+=[[G.edge_order(v1).index(nextedge)]]
        edgelist2+=[[G.edge_order(v2).index(nextedge)]]
    P1=part.Partition(edgelist1)
    P2=part.Partition(edgelist2)
    partitionmap=range(0,len(P1.parts))
    splicemap=W.splicemaps[w.letters[-1]]
    (newP1,newP2)=part.compatible_coarsenings(P1,P2,partitionmap,splicemap)
    if F.is_conjugate_into(w,*W.wordlist): # If w is conjugate into the wordlist then one part of the partition does not correspond to a component, just a segregated edge.
        numberofcomponents=len(newP1.parts)-1
    else:
        numberofcomponents=len(newP1.parts)
    #----------
    if numberofcomponents==1 or (numberofcomponents==2 and crossing_cut_pairs(F,W,w,w)): # if this is true the splitting is trivial
    # the function givesSplitting is a conjuntion of givesCut and not crossingCutPairs. The work of givesCut is to compute numberofcomponents, which we've already done here in extra detail.
        return None
    # otherwise we really do get a splitting over the splitting word, so carry on

        
    partsplicemap=[None]*len(newP1.parts)
    for i in range(len(partsplicemap)):
        for j in newP1.parts[i]: break # take a sample from newP1parts[i]
        partsplicemap[i]=newP2.which_part(splicemap[j]) # see which part of newP2 it splices to

    # components are components just over the fundamental domain for the w action on its axis, not over the infinite w axis, so components may be strictly finer than the partition
    componentstopartition=[None]*len(components)
    for i in range(len(components)):
        verts=set(components[i])
        while verts:
            vert=verts.pop()
            if v2 in G.neighbors(vert):
                e=G[v2][vert].keys()[0]
                componentstopartition[i]=newP2.which_part(G.edge_order(v2).index(e))
                break
        else:
            raise KeyError("something wrong, any vertex in component "+str(i)+" adjacent to v2")

    # big components are vertices that are eventually connected over the whole w axis, though they may not be connected in the fundamental domain of the w action on the w axis.
    # note there is an empty bigcomponent if some component does not have any vertices in the fundamental domain of the w action.
    # indices of bigcomponents are the parts of the newP1 partition
    bigcomponents=[set([])]*len(newP1.parts)
    for i in range(len(components)):
        bigcomponents[componentstopartition[i]]=bigcomponents[componentstopartition[i]]|set(components[i])

    complementsofw=dict()
    for i in range(len(bigcomponents)):
        verticesinthiscomponent=[F.word(c) for c in bigcomponents[i]]
        if not(len(verticesinthiscomponent)==0 and partitionmap[partsplicemap[i]]==i): # if this i is not the part coming from the splitting word closing up upon itself
            complementsofw[i]=[1,verticesinthiscomponent]
    # the keys of complementsofw correspond to indices in bigcomponents, except that we've dropped the non-component corresponding to the splitting word, if there is one.

    def which_complement(inputz):
        # find the big component containing z
        z=F.word(inputz)
        assert(F.degree(w)==1)
        if F.is_power(z,w):
            raise RuntimeError('Input axis equal to splitting axis.')
        zcomponent=None
        zpre=list(z.letters[:1])
        zsuf=list(z.letters[1:])
        # w is cyclically reduced, so the next computations are valid
        while zpre==((w.letters)*(1+len(zpre)//len(w)))[:len(zpre)] or zpre==([-x for x in reversed(w.letters)]*(1+len(zpre)//len(w)))[:len(zpre)]: # while zpre still lies on the w axis
            try:
                nextletter=zsuf.pop(0)
            except IndexError: # if zsuf is empty then reload zsuf with a copy of w
                zsuf=[l for l in w.letters]
                nextletter=zsuf.pop(0)
            try:
                if zpre[-1]==-nextletter:
                    zpre.pop()
                else:
                    zpre.append(nextletter)
            except IndexError:
                zpre.append(nextletter)
                    
        # zpre is the shortest prefix of z*w**n, n>=0  that does not lie on w axis
        # now want to write zpre as w**m*y with m maximal so that y and w**(-1) do not have a common prefix and w is not a prefix of y
        assert(len(zpre))
        if zpre[0]==-w.letters[-1]:
            m=(-len(zpre)//len(w))
        else:
            m=(len(zpre)-1)//len(w) # if len(zpre)==len(w**k) we want k-1
        y=w**(-m)*(F.word(zpre))
        for i in range(len(components)):
            if tuple(y.letters) in components[i]:
                currentpart=componentstopartition[i]
                break
        else:
            raise KeyError("Failed to find component containing input "+str(z)+". y="+str(y)+", zpre="+str(zpre)+", m="+str(m))
        # currentpart is the piece of the partition containing y
        # to get the part containing z we need to apply splicemap m times to get back to fundamental domain
        while m<0:
            currentpart=partsplicemap[currentpart]
            m+=1
        while m>0:
            currentpart=partsplicemap.index(currentpart)
            m-=1
        return currentpart 
        
    # Figure out the w orbits of the components
    # Within each component figure out which elements are w minimal.
    # This is essentially a sorting problem, but with only a partial order.
    # The comparison operation wlessthan is a little slow.
    # Want to make sure we do not spend too much time comparing elements that are not comparable in the partial order.
    # For this reason for each element we compute the projection of the axis of the element to the w axis.
    # We make a dictionary elementsincomplementiwhoseaxisprojectstowj[k]=indices of the elements whose axis projects to point k on the w axis
    # Elements with disjoint projection are not comparable, so we need only pick out the minimal elements for each k.
    worbitofcomplement=dict()
    worbits=dict()
    notinanorbityet=set(complementsofw.keys())
    for thiscomp in complementsofw:
        if thiscomp in notinanorbityet:
            notinanorbityet.remove(thiscomp)
            worbits[thiscomp]=[thiscomp]
            worbitofcomplement[thiscomp]=thiscomp
        nextpart=partsplicemap[partitionmap[thiscomp]]
        ####### nextpart=partsplicemap[worbitofcomplement[thiscomp]]
        while nextpart not in complementsofw: # it may be that nextpart is a part that only has an edge that goes all the way through w without hitting a vertex, so it doesn't correspond to a complementofw
            nextpart=partsplicemap[partitionmap[nextpart]]
            assert(False) ###### i think we've changed setup so that we will never be in this case, and this loop can be deleted. test.
        prefix=F.word([])
        while nextpart!=thiscomp:
            if worbitofcomplement[thiscomp]==thiscomp: # if we're starting from the orbit rep then every nextpart we see is new, otherwise, not.
                #worbits[thiscomp].append(nextpart)
                worbits[thiscomp].insert(0,nextpart)
                notinanorbityet.remove(nextpart)
                worbitofcomplement[nextpart]=thiscomp
            prefix=prefix*w # prefix becomes one higher power of w
            for e in newP2.parts[nextpart]:
                v=G.opposite_end(G.edge_order(v2)[e],v2)
                if v==v1:    
                    pass # don't get any vertices from this edge
                else:
                    for comp in bigcomponents:
                        if v in comp:
                            complementsofw[thiscomp][1].extend([prefix*F.word(v) for v in comp])
                            break
                    else:
                        raise KeyError('did not find vertex in a complementary component')
            nextpart=partsplicemap[partitionmap[nextpart]]
            while nextpart not in complementsofw: # it may be that nextpart is a part that only has an edge that goes all the way through w without hitting a vertex, so it doesn't correspond to a complementofw
                nextpart=partsplicemap[partitionmap[nextpart]]
                assert(False) ###### i think we've changed setup so that we will never be in this case, and this loop can be deleted. test.
        complementsofw[thiscomp][0]=len(worbits[worbitofcomplement[thiscomp]])
    for thisorbit in set(worbitofcomplement.values()):
        worbits[thisorbit]=[worbits[thisorbit][-1]]+worbits[thisorbit][:-1] # It will be convenient to have the orbitrep at index 0
        # right now complementsofw[i][1] contains vertices in the tree adjacent to axis of w in complementary component i
        # the difference from bigcomponents is that complementsofw[i][1] contains all the vertices in the component over the region w**complementsofw[i][0], which records the minimal power of w that stabilizes the component i (which may be different for different i)
    nearbyaxesincomplementsofw=dict()
    elementsincomplementiwhoseaxisprojectstowj=dict()
    for thiscomp in complementsofw:
        nearbyaxesincomplementsofw[thiscomp]=complementsofw[thiscomp]
        elementsincomplementiwhoseaxisprojectstowj[thiscomp]=dict([(j,set([])) for j in range((complementsofw[thiscomp][0])*len(w))])
        newguys=[]
        for thisindex in range(len(complementsofw[thiscomp][1])):
            g=nearbyaxesincomplementsofw[thiscomp][1][thisindex]
            # find all the axes of conjugates of w that go through g
            for l in range(0,len(w)):
                newguy=F.word(g.letters+[-x for x in reversed(w.letters[:l])])
                neww=w.cycle(l) # The w axis through newguy is the same as the neww axis through g.
                if l: # newguy != g
                    newindex=len(nearbyaxesincomplementsofw[thiscomp][1])+len(newguys)
                    newguys.append(newguy)
                else:
                    newindex=thisindex
                elementsincomplementiwhoseaxisprojectstowj[thiscomp][len(g)-1].add(newindex)
                if neww.letters[0]==-g.letters[-1]:# the axis of w and the forward half of the axis of gwG overlap
                    forwardoverlap=0# how far the foward axis of gwG overlaps the axis of w moving forward from the closest point to g on the w axis
                    while len(g)+forwardoverlap < nearbyaxesincomplementsofw[thiscomp][0]*len(neww) and neww.letters[(1+forwardoverlap)%len(neww)]==w.letters[(len(g)-1+forwardoverlap)%len(neww)]:
                        forwardoverlap+=1
                        elementsincomplementiwhoseaxisprojectstowj[thiscomp][len(g)-1+forwardoverlap].add(newindex)
                    backwardoverlap=0# how far the forward axis of gwG overlaps the axis of w moving backward from the closest point to g on the w axis
                    while len(g)-2-backwardoverlap >=0 and neww.letters[(1+backwardoverlap)%len(neww)]==-w.letters[(len(g)-2-backwardoverlap)%len(neww)]:
                        backwardoverlap+=1
                        elementsincomplementiwhoseaxisprojectstowj[thiscomp][len(g)-1-backwardoverlap].add(newindex)
                elif neww.letters[-1]==g.letters[-1]: # same as before but now its the backwards half of the gwG axis overlapping the axis of w
                    forwardoverlap=0
                    while len(g)+forwardoverlap < nearbyaxesincomplementsofw[thiscomp][0]*len(neww) and neww.letters[(-2-forwardoverlap)%len(neww)]==-w.letters[(len(g)-1+forwardoverlap)%len(neww)]:
                        forwardoverlap+=1
                        elementsincomplementiwhoseaxisprojectstowj[thiscomp][len(g)-1+forwardoverlap].add(newindex)
                    backwardoverlap=0
                    while len(g)-2-backwardoverlap >=0 and neww.letters[(-2-backwardoverlap)%len(neww)]==w.letters[(len(g)-2-backwardoverlap)%len(neww)]:
                        backwardoverlap+=1
                        elementsincomplementiwhoseaxisprojectstowj[thiscomp][len(g)-1-backwardoverlap].add(newindex)
        nearbyaxesincomplementsofw[thiscomp][1].extend(newguys)
        
        
    # now for each complement i in nearbyaxesincomplementsofw[i][1] a list of group elements g such that the axis of gwG passes within distance 1 of the axis of w and is in component i
    # reduce to subset of these that is minimal with respect to separation from w, ie, axes that are not separated from the axis of w by anything else in the collection.
    # These will be globablly minimal since no other axis could possibly separate if it didn't come at least as close to axis of w.
    #####   This seems to be a major bottleneck.
    ##### Is there a faster way to pick out the minimal elements of a partially ordered set?

    def wlessthan(inputx,inputy):
        """
        True if axis of xwX separates axis of w from axis of ywY.
        """
        x=F.word(inputx)
        y=F.word(inputy)
        X=x**(-1)
        z=X*y
        e=F.word([])
        if F.is_power(x,w) or F.is_power(y,w) or F.is_power(X*y,w): # degenerate cases, not three distinct axes
            return False

        # for convenience, translate by w until X and z minimally do not begin with winverse
        while X.letters[0]==w.letters[0] and z.letters[0]==w.letters[0]:
            X=w**(-1)*X
            z=w**(-1)*z
        while X.letters[0]==-w.letters[-1] or z.letters[0]==-w.letters[-1]:
            X=w*X
            z=w*z
        return which_complement(X)!=which_complement(z)

    def run_the_gauntlet(groupelementlist,champions,challenger):
        """
        Compare challenger to each champion. If challenger < champion then champion is removed from champions list. If challenger > some champion then stop and return False.
        If challenger not > champion for all champions then return True.
        """
        newchampion=False
        if not champions:
            newchampion=True
        else:
            fallen=[]
            # challenger tested against each champion. 
            for k in range(len(champions)-1,-1,-1):
                if groupelementlist[champions[k]]==groupelementlist[challenger]: # this list may have had repetitions
                    break
                elif wlessthan(groupelementlist[champions[k]],groupelementlist[challenger]): # challenger defeated
                    break
                elif wlessthan(groupelementlist[challenger], groupelementlist[champions[k]]): # challenger defeats champion and continues the trials. This challenger is guaranteed to be a champion now, but may defeat other champions too.
                    fallen.append(k)
                else:
                    pass # they are not comparable in the partial order, continue the trials
            else:
                newchampion=True # challenger went undefeated, so is a new champion
                while fallen:
                    x=fallen.pop(0)
                    champions.pop(x) # remove any fallen champions
        return newchampion

    def tournament_of_champions(groupelementlist,challengers):
        divideandcounquerparameter=16 # This is to balance the search savings of 
        if len(challengers)<2:
            return challengers
        elif len(challengers)<divideandcounquerparameter:
            champions=[challengers.pop(0)]
            while challengers:
                challenger=challengers.pop(0)
                if run_the_gauntlet(groupelementlist,champions, challenger):
                    champions.append(challenger)
            return champions
        else:
            pool1=challengers[:len(challengers)//2]
            pool2=challengers[len(challengers)//2:]
            champions1=tournament_of_champions(groupelementlist,pool1)
            champions2=tournament_of_champions(groupelementlist,pool2)
            champions3=[]
            while champions2:
                challenger=champions2.pop(0)
                if run_the_gauntlet(groupelementlist,champions1,challenger):
                    champions3.append(challenger)
            return champions1+champions3
            
    indexgiveswminmialelement=dict()
    if verbose or maxnumberof2componentcutstoconsider:
        numberofaxes=0
        for i in range(len(nearbyaxesincomplementsofw)):
            numberofaxes+=len(nearbyaxesincomplementsofw[i][1])
    if maxnumberof2componentcutstoconsider:
        if numberofaxes>maxnumberof2componentcutstoconsider:
            raise TooBigError(str(numberofaxes)+" potential cut pairs is beyond limit set by 'maxnumberof2componentcutstoconsider'",numberofaxes)
    if verbose:
        print "Splitting has "+str(len(nearbyaxesincomplementsofw))+" components and "+str(numberofaxes)+" axes. Finding minimal axes in each component."
    minimalaxesincomplementsofw=dict()
    for i in complementsofw:
        minimalaxesincomplementsofw[i]=[nearbyaxesincomplementsofw[i][0],[]]
        indexgiveswminmialelement[i]=dict([(k,True) for k in range(len(nearbyaxesincomplementsofw[i][1]))])
        if verbose:
            print "Component "+str(1+i)+" has "+str(len(nearbyaxesincomplementsofw[i][1]))+" axes in "+str(len(elementsincomplementiwhoseaxisprojectstowj[i]))+" groups."
        for j in elementsincomplementiwhoseaxisprojectstowj[i]:
            challengers=list(elementsincomplementiwhoseaxisprojectstowj[i][j])
            champions=tournament_of_champions(nearbyaxesincomplementsofw[i][1], challengers) # this is a list of indices in the list nearbyaxesincomplementsofw[i][1] such that the axis of the corresponding element projects to the point at distance j along the w axis, and such that the element is w minmal among all such elements
                                                         # the point here is that elements whose axes have disjoint projections to the w axis are not w comparable, so we should not waste time trying to compare them
            if verbose:
                print "Of "+str(len(elementsincomplementiwhoseaxisprojectstowj[i][j]))+" axes considered in group "+str(1+j)+" of component "+str(1+i)+", "+str(len(champions))+" are minimal."
            losers=set(elementsincomplementiwhoseaxisprojectstowj[i][j])-set(champions)
            for loser in losers:
                indexgiveswminmialelement[i][loser]=False
        minimalaxesincomplementsofw[i][1]=[(nearbyaxesincomplementsofw[i][1][k],which_complement([-x for x in reversed(nearbyaxesincomplementsofw[i][1][k].letters)])) for k in range(len(nearbyaxesincomplementsofw[i][1])) if indexgiveswminmialelement[i][k]] 

    # minimalaxesincomplementsofw[i][0] is the minimal power of w stabilizing component i
    # and minimalaxesincomplementsofw[i][1] is a list of pairs (g,c) where  g is a group element in component i so that the axis of gwG is w minimal and passes within distance 1 of the axis of w and c is the component containing g^{-1}.
    # In Otal's construction of Bass-Serre tree for the splitting of F over <w> we have a bipartite tree where one class of vertex is stabilized by conjugates of <w> and the other class corresponds to collections of unseperable axes
    # minimalaxesincomplementsofw records these collections of unseperable axes. The axes of gwG for g in minimalaxesincomplementsofw[i][1] are representatives of the w orbits of the unseparable axes of component i. The c tells how to transport the numbering of the complements of w to the complements of gwG. Thus, we have encoded the Bass-Serre tree of the splitting.
    # Now we have to figure out the quotient graph.

    # ----- for debugging
    componentconnections=dict()
    for i in minimalaxesincomplementsofw:
         componentconnections[i]=set([axis[1] for axis in minimalaxesincomplementsofw[i][1]])
         # assert(all([(y in componentconnections[x])==(x in componentconnections[y]) for (x,y) in [(x,y) for x in minimalaxesincomplementsofw for y in minimalaxesincomplementsofw]]))
         # If this assertion fails it may mean some non-minimal axis did not get weeded out. This will result in a too large vertex stabilizer. But this is not always the case.
    # -----

    if verbose:
        print "Computing quotient graph."
    quotientgraph=nx.MultiDiGraph()
    quotientgraph.add_node(smash(nameprefix,w()), {'stabilizer':set([tuple(w.letters)])}) # one vertex with a cyclic stabilizer from the splitting word
    edgestobe=set(worbits.keys()) # the set of w orbits of complementary components of w. quotient graph of group has one edge per w orbit. Has one vertex per F orbit.
    while edgestobe:
        thisedge=edgestobe.pop()
        quotientgraph.add_node(smash(nameprefix,thisedge), stabilizer=set([tuple((w**(minimalaxesincomplementsofw[thisedge][0])).letters)]))
        quotientgraph.add_edge(smash(nameprefix,thisedge),smash(nameprefix,w()),smash(nameprefix,thisedge),label=F.word([]),headstabilizer=set([w**(minimalaxesincomplementsofw[thisedge][0])]), tailstabilizer=set([w**(minimalaxesincomplementsofw[thisedge][0])]))
        minimaladjacentaxes=set(minimalaxesincomplementsofw[thisedge][1])
        while minimaladjacentaxes:
            nextaxis=minimaladjacentaxes.pop()
            g=nextaxis[0] # the group element that takes axis of w to nextaxis, ie, nextaxis=g*axis(w)
            c=nextaxis[1] # the component containing g inverse.
                          # from point of view of nextaxis, axis(w) is in component c
            x = worbits[worbitofcomplement[c]].index(c) # x is minimal non-negative integer such that w**x takes the w-orbit rep for the w-orbit of c to c
            if worbitofcomplement[c] in edgestobe: # if this is the first time we've seen this orbit of component
                edgestobe.remove(worbitofcomplement[c])
                # This is the first time we've seen this component. Add the stabilizer of nextaxis to the vertex stabilizer of the quotient.
                quotientgraph.node[smash(nameprefix,thisedge)]['stabilizer'].add(tuple((g*(w**(minimalaxesincomplementsofw[c][0]))*g**(-1)).letters))
                quotientgraph.add_edge(smash(nameprefix,thisedge),smash(nameprefix,w()),smash(nameprefix,worbitofcomplement[c]),label=g*(w**(x)),headstabilizer=set([w**minimalaxesincomplementsofw[worbitofcomplement[c]][0]]), tailstabilizer=set([g*w**minimalaxesincomplementsofw[c][0]*g**(-1)]))
                for (h,d) in minimalaxesincomplementsofw[worbitofcomplement[c]][1]:
                    minimaladjacentaxes.add((g*w**x*h,d)) # need to do this because so far we only know about the minmimaladjacent axes that are distance at most 1 from axis of w. This adds the group elements we need to get minmial adjacent axes that are distance at most 1 from those.
            else:
                # we've been here before, so an edge has already been added to the quotientgraph. Just need to add to the stabilizer.
                h=quotientgraph[smash(nameprefix,thisedge)][smash(nameprefix,w())][smash(nameprefix,worbitofcomplement[c])]['label']
                quotientgraph.node[smash(nameprefix,thisedge)]['stabilizer'].add(tuple((g*w**x*(h)**(-1)).letters))

    # We've got the quotient graph along with a set of generators for each stabilizers. Make it into a graph of groups.
    # We need to do this because for quotientgraph we built up stabilizers one generator at a time, but for a graph of groups we need to know the whole stabilizer subgroup.
    qgog=gog.FPGraphOfGroups()
    for v in quotientgraph.nodes():
        thisstabwl=[F.word(s) for s in quotientgraph.node[v]['stabilizer']]
        thisstabilizer=freegroup.FGSubgroupOfFreeFrom(F,thisstabwl, generatorbasename=str(v), displaystyle=list)
        qgog.add_vertex(v,vertgroup=thisstabilizer)
    for e in quotientgraph.edges(keys=True):
        ogroup=qgog.node[e[0]]['group']
        tgroup=qgog.node[e[1]]['group']
        originword=F.word(quotientgraph[e[0]][e[1]][e[2]]['tailstabilizer'].pop())
        terminusword=F.word(quotientgraph[e[0]][e[1]][e[2]]['headstabilizer'].pop())
        assert(len(originword))
        assert(len(terminusword))
        edgegroup=freegroup.FGSubgroupOfFree(F,[originword],generatorbasename='e'+str(e[2]), displaystyle=list)
        orestricted=ogroup.restrict_word(originword)
        trestricted=tgroup.restrict_word(terminusword)
        omap=group.Homomorphism(edgegroup,ogroup,dict([(1,orestricted)]))
        tmap=group.Homomorphism(edgegroup,tgroup,dict([(1,trestricted)]))
        qgog.add_edge(e[0],e[1],edgegroup,omap,tmap,label=quotientgraph[e[0]][e[1]][e[2]]['label'],key=e[2]) 

    # figure out where words of the original wordlist go in qgog
    ##### figure out where generators of the free group go in qgog
    wordmap=[]
    for thisword in wordlist:
        for vert in qgog:
            conjugateword=qgog.localgroup(vert).find_conjugate_in_subgroup(thisword)
            if conjugateword is not None:
                wordmap.append((vert,conjugateword,1))
                break
        else:
            raise KeyError(thisword()+" is not elliptic")
    if extrawordlist:
        extrawordmap=[]
        for thisword in extrawordlist:
            for vert in qgog:
                conjugateword=qgog.localgroup(vert).find_conjugate_in_subgroup(thisword)
                if conjugateword is not None:
                    extrawordmap.append((vert,conjugateword,1))
                    break
            else:
                raise KeyError(thisword()+" is not elliptic")
        return qgog, wordmap, extrawordmap
    return qgog, wordmap

freegroup.FGFreeGroup.get_relative_cyclic_splitting_over=get_relative_cyclic_splitting_over 
                




def get_induced_multiword(thisgog, thisvert, multiwordmap, simplifyandminimize=False, blind=False,withedgemap=False):
    """
    Find the induced multiword in the vertex group of vert.
    multiwordmap is a list [(v,w),...] so that v is a vertex of thisgog and w is a word in thisgog.localgroup(v).
    If simplifyandminimize=True then the resulting inducedmultiword is simplified and minimized.
    """
    inducedmultiword=[]
    edgeupdate={}
    thisgroup=thisgog.localgroup(thisvert)
    outedges=set(thisgog.out_edges(thisvert,keys=True))
    inedges=set(thisgog.in_edges(thisvert,keys=True))
    loops=outedges&inedges
    outedges=outedges-loops
    inedges=inedges-loops
    images=dict()

    # build the induced multiword
    for e in outedges:
        theindex=len(inducedmultiword)
        try:
            inducedmultiword.append(thisgog.getomap(e)(thisgog.localgroup(e).word([1])))
        except IndexError:
            inducedmultiword.append(thisgroup.word([]))
        images[(e,0)]=(theindex,1)
    for e in inedges:
        theindex=len(inducedmultiword)
        try:
            inducedmultiword.append(thisgog.gettmap(e)(thisgog.localgroup(e).word([1])))
        except IndexError:
            inducedmultiword.append(thisgroup.word([]))
        images[(e,1)]=(theindex,1)
    for e in loops:
        theindex=len(inducedmultiword)
        try:
            inducedmultiword.append(thisgog.getomap(e)(thisgog.localgroup(e).word([1])))
        except IndexError:
            inducedmultiword.append(thisgroup.word([]))
        try:
            inducedmultiword.append(thisgog.gettmap(e)(thisgog.localgroup(e).word([1])))
        except IndexError:
            inducedmultiword.append(thisgroup.word([]))
        images[(e,0)]=(theindex,1)
        images[(e,1)]=(1+theindex,1)
    for i in range(len(multiwordmap)):
        if multiwordmap[i][0]==thisvert:
            theindex=len(inducedmultiword)
            inducedmultiword.append(multiwordmap[i][1])
            images[i]=(theindex,1)
            
    if not simplifyandminimize:
        if not withedgemap:
            return inducedmultiword
        else:
            return inducedmultiword, images
    else:
        smwl=wg.simplify_and_minimize_wordlist(thisgroup, inducedmultiword, blind=blind)
        mswordmap=smwl['wordmap']
        if not withedgemap:
            return smwl['wordlist']
        else:
            for k in images:
                images[k]=(mswordmap[images[k][0]][0], mswordmap[images[k][0]][1] * images[k][1])
            if blind:
                return smwl['wordlist'], images
            else:
                return smwl['wordlist'], images, smwl['minimizingautomorphism']
        
def is_RJSJ(F,wlmap,thisgog, verbose=False):
    """
    Check if the given graph of groups thisgog is the RJSJ for F rel a wordlist.
    wlmap is a list with entries (vertex of thisgog, word in the vertex group, power)
    """
    rigidverts=[]
    circleverts=[]
    degree2verts=[]
    for vert in thisgog.nodes():
        if thisgog.localgroup(vert).rank>1:
            iwl=get_induced_multiword(thisgog,vert, wlmap,simplifyandminimize=True)
            if thisgog.localgroup(vert).isCircle(iwl):
                circleverts.append(vert)
            elif thisgog.localgroup(vert).is_rigid_rel(iwl):
                rigidverts.append(vert)
            else:
                if verbose:
                    print "vertex "+str(vert)+" is neither rigid nor circle."
                return False
        elif thisgog.localgroup(vert).rank==1:
            totaldegree=0
            for e in thisgog.out_edges(vert,keys=True):
                if thisgog.localgroup(e).rank:
                    totaldegree+=thisgog.localgroup(vert).degree(thisgog.getomap(e)(thisgog.localgroup(e).word([1])))
            for e in thisgog.in_edges(vert,keys=True):
                if thisgog.localgroup(e).rank:
                    totaldegree+=thisgog.localgroup(vert).degree(thisgog.gettmap(e)(thisgog.localgroup(e).word([1])))
            if totaldegree==2:
                degree2verts.append(vert)

    # if still going then all nonabelian vertex groups are rigid or circles. This is not the RJSJ if there is a degree two vertex adjacent only to circles, unless it is a cutpoint.
    toosplit=False
    for vert in degree2verts:
        if not F.is_conjugate_into(thisgog.localgroup(vert).word([1]), *[w for (v,w,p) in wlmap]):
            if all([v in circleverts for v in thisgog.neighbors(vert)]):
                if verbose:
                    print "vertex "+str(vert)+" should be inside a larger circle"
                toosplit=True
                break
    return not toosplit
freegroup.FGFreeGroup.is_RJSJ=is_RJSJ

def simplify_GOG(thisgog, thisvert, multiwordmap):
    # What is this function doing here? This stuff was done directly in get_RJSJ. Factor it out?
    """
    thisgog is a graph of free groups with cyclic edge groups. multiwordmap is a list [(v,w),...] where v is a vertex of thisgog and w is a word in thisgog.localgroup(v). Change the graph of groups structure at vert so that the induced multiword is simplified and Whitehead minimized.
    That
    First change the graph of groups structure so that edge groups that are conjugate into the same cyclic subgroup actually map into the same cyclic subgroup. Then Whitehead minimize the inducedmultiword at the vertex and update the vertex group and edge maps. Update the multiwordmap to reflect the Whitehead minimization.
    """
    inducedmultiword, images, minimizingautomorphism, inverseminimizingautomorphism = getInducedMultiword(thisgog, thisvert, multiwordmap, simplifyandminimize=True, withedgemap=True)
    # now update the multiwordmap
    for i in range(len(multiwordmap)):
        if multiwordmap[i][0]==thisvert:
            multiwordmap[i]=(thisvert, inducedmultiword[images[i][0]], multiwordmap[i][1]*images[i][1])
    thisgog.change_vertex_group_by_automorphism(thisvert,minimizingautomorphism)
    edges=[(e,0) for e in thisgog.out_edges(thisvert,keys=True)]+[(e,1) for e in thisgog.in_edges(thisvert,keys=True)]
    for edge in edges:
        # the image of the edge generator is conjugate to a power of images[edge]. Need to find the conjugator and adjust the edgemap.
        if edge[1]==0:
            edgeimage=thisgog.getomap(edge[0]).images[0]
        elif edge[1]==1:
            edgeimage=thisgog.gettmap(edge[0]).images[0]
        else:
            raise KeyError() 
        conjugator=thisgroup.is_conjugate_into(edgeimage,inducedmultiword[images[edge][0]])['conjugator']
        thisgog.change_edge_map(edge[0],edge[1], conjugator)
    

def get_RJSJ(F,whiteheadgraphorwordlist,withmap=False, printresult=False, nameprefix='', blind=True, simplified=False, minimized=False, verbose=False,cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None):
    """
    Find the JSJ splitting for F relative to the worldlist represented as whiteheadgraphorwordlist.
    """
    # Algorithm as follows:
    # 0. Find all the splitting elements.
    # 1. Create a graph of groups with a single vertex stabilized by F, and the original wordlist as inducedmultiword.
    # 2. Terminate if there are no unused splitting elements.
    # 3. Take an unused splitting element. Figure out which vertex group it is conjugate into.
    # 4. Split that vertex group over the splitting element.
    # 5. Refine the graph of groups by inserting the splitting from 4 into the graph in place of the vertex.
    # 6. goto 2

    assert(blind) # blind=False not implemented yet - have to set up tracking of hyperbolic elements
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, simplified=simplified, blind=blind, verbose=verbose,minimized=minimized, simplifyandminimize=True)
    W=wgp['WhiteheadGraph']
    wordlist=wgp['wordlist']
    wordmap=wgp['wordmap']
    assert(W.is_connected()) # if it's not connected must first pass to free factors

    rjsj=gog.FPGraphOfGroups()
    thisvert=smash(nameprefix,0)
    rjsj.add_vertex(thisvert,F)
    wheredidmywordsgo=[]
    for i in range(len(wordmap)):
        wheredidmywordsgo.append((thisvert,wordlist[wordmap[i][0]],wordmap[i][1]))
    universal_split_vertex(rjsj,thisvert,wheredidmywordsgo,MinNumComponents=3,verbose=verbose, cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider) # performs all splittings corresponding to cut points or uncrossed cut pairs with at least 3 components
    firstroundverts=[n for n in rjsj.nodes()]
    for thisvert in firstroundverts: # now for each higher rank vertex try to split it over uncrossed cut pairs with 2 components
        if rjsj.localgroup(thisvert).rank>1:
            universal_split_vertex(rjsj,thisvert,wheredidmywordsgo,MinNumComponents=2,verbose=verbose, cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider)
    if printresult:
        print rjsj
        if withmap and type(whiteheadgraphorwordlist)==list:
            for i in range(len(wheredidmywordsgo)):
                print "The word "+(whiteheadgraphorwordlist[i])()+" corresponds to the word "+(wheredidmywordsgo[i][1])()+" in vertex "+str(wheredidmywordsgo[i][0])
    if withmap:
        return rjsj, wheredidmywordsgo
    else:
        return rjsj
freegroup.FGFreeGroup.get_RJSJ=get_RJSJ

def universal_split_vertex(thisgog,thisvert,wheredidmywordsgo, MinNumComponents=2,verbose=False, cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None):
    """
    Refine thisgog by performing universal splittings of thisvert.
    """
    inducedmultiword,images,inducedmultiwordminimizingautomorphism =get_induced_multiword(thisgog,thisvert,wheredidmywordsgo, blind=False, simplifyandminimize=True,withedgemap=True)
    thisgog.change_vertex_group_by_automorphism(thisvert,inducedmultiwordminimizingautomorphism)
    for i in range(len(wheredidmywordsgo)):
        if wheredidmywordsgo[i][0]==thisvert:
            wheredidmywordsgo[i]=(thisvert, inducedmultiwordminimizingautomorphism(wheredidmywordsgo[i][1]),wheredidmywordsgo[i][2])
    W=wg.WGraph(inducedmultiword, simplified=True, autominimize=False)
    if W.is_circle(): # we're done, no universal splittings
        pass
    else:
        if verbose:
             print "Searching for cut pairs."
        if MinNumComponents>2:
            DoNotVerifyTwoComponentWords=True
        else:
            DoNotVerifyTwoComponentWords=False
        cuts,surethatsall=find_universal_splitting_words(thisgog.localgroup(thisvert), W, inducedmultiword, DoNotVerifyTwoComponentWords=DoNotVerifyTwoComponentWords, MinNumComponents=MinNumComponents, simplified=True,minimized=True,verbose=verbose,cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider)
        if not surethatsall:
            raise TooBigError("Could not determine with 'cutpairsearchrecursionlimit'="+str(cutpairsearchrecursionlimit))
        splittingelements=cuts['cutpoints']|cuts['uncrossed']
        splittingelementsbyvertex=dict() # At first all of the splitting words are in thisvert, but each time we split the splitting words will be in one of the new vertex stabilizers, so we use a dict to keep track of where the splitting elements are at each step
        if splittingelements:
            splittingelementsbyvertex[thisvert]=splittingelements
        while splittingelementsbyvertex: 
            if verbose:
                print "Performing next splitting."
            newvert=splittingelementsbyvertex.keys()[0]
            thiscut=splittingelementsbyvertex[newvert].pop()
            if not splittingelementsbyvertex[newvert]:
                del splittingelementsbyvertex[newvert]
            split_and_refine(thisgog,newvert,thiscut,wheredidmywordsgo, splittingelementsbyvertex, verbose=verbose, cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider)
            # function refines thisgog, no return

def split_and_refine(thisgog, thisvert, thiscut, wheredidmywordsgo, splittingelementsbyvertex, verbose=False, cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None):
    """
    Refine thisgog by splitting the local group of thisvert over the word thiscut relative to incident edges.
    """
    # split thisvert over thiscut relative to multiword and incident edges
    thisgroup=thisgog.localgroup(thisvert)
    assert(thisgroup is thiscut.group)
    outedges=set(thisgog.out_edges(thisvert,keys=True))
    inedges=set(thisgog.in_edges(thisvert,keys=True))
    loops=inedges&outedges
    outedges=outedges-loops
    inedges=inedges-loops

    inducedmultiword,images,inducedmultiwordminimizingautomorphism =get_induced_multiword(thisgog,thisvert,wheredidmywordsgo, blind=False, simplifyandminimize=True,withedgemap=True)
    # we want the inducedmultiword to be simplified and minimized, which means we may need to apply an automorphism to the vertex group at thisvertex. We should then change the incident edge maps and the splitting elements at this vertex.
    thisgog.change_vertex_group_by_automorphism(thisvert,inducedmultiwordminimizingautomorphism)
    thiscutconjugateword=thisgroup.cyclic_reduce(inducedmultiwordminimizingautomorphism(thiscut))       
    if thisvert in splittingelementsbyvertex:
        newcutsinthisvert=set([inducedmultiwordminimizingautomorphism(w) for w in splittingelementsbyvertex[thisvert]])
        splittingelementsbyvertex[thisvert]=newcutsinthisvert
            
    inducedW=wg.WGraph(inducedmultiword, simplified=True, autominimize=False)
        # images is a dict whose keys are  of the form (edge,0) or (edge,1) or i and whose values are interpreted:
        # images[(edge,0)]=(i,p) means the omap for edge e maps the generator of the edge group to the word inducedmultiword[i]**p
        # images[j]=(i,p) means the word from wheredidmywordsgo[i][0] is conjugate to inducedmultiword[i]**p
        
        # find the splitting of the vertex group relative to the inducedmultiword
        # refinedvert is a graph of groups decomposition of the thisvert
        # emap is a list so that emap[i]=(v,w,p) means that word i of the inducedmultiword is conjugate to the word w**p in the group of vertex v in refinedvert


    if thisvert in splittingelementsbyvertex: # there are more splitting words in thisgroup so make sure to track where they go when we refine thisgog
        exwordlist=list(splittingelementsbyvertex[thisvert])
        refinedvert, emap, extrawordmap=thisgroup.get_relative_cyclic_splitting_over(inducedW,inducedmultiword, thiscutconjugateword, thisvert, extrawordlist=exwordlist,verbose=verbose, cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider, simplified=True, minimized=True)
        for (newvert,newsplittingelement,newpower) in extrawordmap:
            try:
                splittingelementsbyvertex[newvert].add(newsplittingelement)
            except KeyError:
                splittingelementsbyvertex[newvert]=set([newsplittingelement])
        del splittingelementsbyvertex[thisvert] # any splitting elements that were in thisvert are now contained in one of the refinement vertices
    else: #there are no further splitting words in thisgroup so this is the only refinement of thisvert
        refinedvert, emap=thisgroup.get_relative_cyclic_splitting_over(inducedW,inducedmultiword, thiscutconjugateword, thisvert, verbose=verbose, cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider, simplified=True, minimized=True)

    for i in range(len(wheredidmywordsgo)):
        if wheredidmywordsgo[i][0]==thisvert: # words conjugate to thisvert are now conjugate into one of the verts of refinedvert
            wheredidmywordsgo[i]=(emap[images[i][0]][0], emap[images[i][0]][1], emap[images[i][0]][2]*wheredidmywordsgo[i][2])
        
    # figure out how the edges of the original graph are going to connect to the refined vertex
    edgeupdate={}
    for e in inedges:
        thevertex=emap[images[(e,1)][0]][0]
        theword=emap[images[(e,1)][0]][1]
        thepower=emap[images[(e,1)][0]][2]*images[(e,1)][1]
        edgeupdate[e]=(thevertex, group.Homomorphism(thisgog.localgroup(e),refinedvert.localgroup(thevertex),dict([(1,(theword)**(thepower))])))
    for e in outedges:
        thevertex=emap[images[(e,0)][0]][0]
        theword=emap[images[(e,0)][0]][1]
        thepower=emap[images[(e,0)][0]][2]*images[(e,0)][1]
        edgeupdate[e]=(thevertex, group.Homomorphism(thisgog.localgroup(e),refinedvert.localgroup(thevertex),dict([(1,(theword)**(thepower))])))
    for e in loops:
        thevertex=emap[images[(e,0)][0]][0]
        theword0=emap[images[(e,0)][0]][1]
        theword1=emap[images[(e,1)][0]][1]
        thepower0=emap[images[(e,0)][0]][2]*images[(e,0)][1]
        thepower1=emap[images[(e,1)][0]][2]*images[(e,0)][1]
        edgeupdate[e]=((thevertex,thevertex),(group.Homomorphism(thisgog.localgroup(e),refinedvert.localgroup(thevertex),dict([(1,(theword0)**(thepower0))])), group.Homomorphism(thisgog.localgroup(e),refinedvert.localgroup(thevertex),dict([(1,(theword1)**(thepower1))]))))
            
    if verbose:
        print "Refining the splitting."
    # replace thisvert by a graph of groups compatible with previous splitting
    thisgog.refine_vertex(thisvert, refinedvert, edgeupdate)
    # the function modifies thisgog, wheredidmywordsgo, and splittingelementsbyvertex
    # does not return anything


def get_max_free_and_cyclic_splitting_rel(F, whiteheadgraphorwordlist, withmap=False, printresult=False, verbose=False, cutpairsearchrecursionlimit=None, maxnumberof2componentcutstoconsider=None, simplified=False, minimized=False, blind=True):
    assert(blind) # blind=False not implemented yet
    wgp=wg.wgparse(F,whiteheadgraphorwordlist, blind=blind,simplified=simplified, minimized=minimized,simplifyandminimize=True, verbose=verbose)
    W=wgp['WhiteheadGraph']
    wordlist=wgp['wordlist']
    wordmap=wgp['wordmap']
    if verbose:
        print "Looking for free splittings."
    freesplitting,wmap=F.get_free_splitting_rel(wordlist, withwordmap=True, minimized=True, simplified=True, blind=blind)
    wheredidmywordsgo=[(wmap[wordmap[i][0]][0], wmap[wordmap[i][0]][1],wmap[wordmap[i][0]][2]*wordmap[i][1]) for i in range(len(wordmap))]
    higherrankvertices=[v for v in freesplitting.nodes() if freesplitting.localgroup(v).rank>1]
    if verbose:
        print "Found a free splitting with "+str(len(higherrankvertices))+" higher rank vertices."
    for thisvert in higherrankvertices: # find cyclic splittings of the vertex groups
        if verbose:
            print "Finding cyclic splittings of vertex "+str(1+higherrankvertices.index(thisvert))+"."
        thisgroup=freesplitting.localgroup(thisvert)
        indexofthiswordlistintomainwordlist=[]
        thiswordlist=[]
        for i in range(len(wheredidmywordsgo)):
            if wheredidmywordsgo[i][0]==thisvert:
                indexofthiswordlistintomainwordlist.append(i)
                thiswordlist.append(wheredidmywordsgo[i][1])
        thisrjsj,thiswmap=thisgroup.get_RJSJ(thiswordlist,withmap=True, nameprefix=thisvert, verbose=verbose, simplified=True, minimized=True, blind=blind, cutpairsearchrecursionlimit=cutpairsearchrecursionlimit, maxnumberof2componentcutstoconsider=maxnumberof2componentcutstoconsider)
        for j in range(len(thiswordlist)):
            wheredidmywordsgo[indexofthiswordlistintomainwordlist[j]]=thiswmap[j]
        # To refine the freesplitting we need to know how to attach edges from freesplitting incident to thisvert to thisrjsj.
        # Since the stabilizers of edges in freesplitting are trivial it doesn't matter how we attach the edges.
        # Take any vertex somevertex in thisrjsj and attach all the edges there, with trivial edge maps.
        somevertex=thisrjsj.nodes()[0]
        edgeupdate=dict()
        outedges=set(freesplitting.out_edges(thisvert,keys=True))
        inedges=set(freesplitting.in_edges(thisvert,keys=True))
        loops=outedges&inedges
        nonloops=(inedges|outedges)-loops
        edgeupdate=dict()
        for e in nonloops:
            edgeupdate[e]=(somevertex, group.Homomorphism(freesplitting.localgroup(e), thisrjsj.localgroup(somevertex)))
        for e in loops:
            edgeupdate[e]=((somevertex,somevertex), (group.Homomorphism(freesplitting.localgroup(e), thisrjsj.localgroup(somevertex)),group.Homomorphism(freesplitting.localgroup(e), thisrjsj.localgroup(somevertex))))
        freesplitting.refineVertex(thisvert, thisrjsj, edgeupdate)
    if verbose:
        print ""
    if printresult:
        print freesplitting
        if withmap:
            for i in range(len(wheredidmywordsgo)):
                    print "The word "+(wordlist[i])()+" corresponds to a power of "+(wheredidmywordsgo[i][1])()+" in vertex "+str(wheredidmywordsgo[i][0])
    if withmap:
        return freesplitting, wheredidmywordsgo
    else:
        return freesplitting
freegroup.FGFreeGroup.get_max_free_and_cyclic_splitting_rel=get_max_free_and_cyclic_splitting_rel

    
        

def missing_3_letter_subwords(*wordlist):
    """
    Return the set of 3 letter words in F that do not occur as subwords in w or w**(-1) as cyclic words.
    """
    thegroup=wordlist[0].group
    wordgenerator=enumeratewords.generate_words(thegroup,3,3)
    missingWords=set([tuple(x.letters) for x in wordgenerator])
    for theword in wordlist:
        reducedword=thegroup.cyclic_reduce(theword)
        for i in range(len(reducedword)):
            missingWords.discard((reducedword.letters[i],reducedword.letters[(i+1)%len(reducedword)],reducedword.letters[(i+2)%len(reducedword)]))
            missingWords.discard((-reducedword.letters[(i+2)%len(reducedword)],-reducedword.letters[(i+1)%len(reducedword)],-reducedword.letters[i]))
    return [thegroup.word(lets) for lets in missingWords]
    
def contains_all_3_letter_subwords(*wordlist):
    thegroup=wordlist[0].group
    total3LetterSubwords=2*thegroup.rank*(2*thegroup.rank-1)**2
    the3LetterSubwords=set()
    for theword in wordlist:
        reducedword=thegroup.cyclic_reduce(theword)
        for i in range(len(reducedword)):
            the3LetterSubwords.add((reducedword.letters[i],reducedword.letters[(i+1)%len(reducedword)],reducedword.letters[(i+2)%len(reducedword)]))
            the3LetterSubwords.add((-reducedword.letters[(i+2)%len(reducedword)],-reducedword.letters[(i+1)%len(reducedword)],-reducedword.letters[i]))
    assert(len(the3LetterSubwords)<=total3LetterSubwords)
    return len(the3LetterSubwords)==total3LetterSubwords
        
