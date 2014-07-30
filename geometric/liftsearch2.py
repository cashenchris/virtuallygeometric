import re, copy
import heegaard
import gengraphs
import subgroup


def look_for_good_cover(wordlist,rank,degree,Heegaardwaittime=10,verbose=True):
    """
    Try to find a cover in which lifts of things in wordlist are embedded.
    """
    working = []
    for graph in gengraphs.graphs(rank,degree):
        lifts = []
        H = subgroup.subgroup(graph)
        for word in wordlist:
            lifts = lifts+H.lifts(word)
        strlifts = map(lambda x: x.alpha(), lifts)
        (yesno, fulltext) = heegaard.is_realizable(strlifts,full_answer=True, maxtime=Heegaardwaittime)
        if yesno:
            if re.search("Unable to determine",fulltext)==None:
                if verbose:
                    print H.graph.edges()
                    print '\n Found a subgroup which seems to work.'
                tempH = copy.copy(H.graph)
                working.append(tempH)
            else:
                if verbose:
                    print 'Unable to determine...'
    return working

def alldegrees(wordlist,rank):
    """
    Do look_for_good_cover for degree 2, then degree 3, etc.  until you find
    something.
    """
    ans = [];d=2
    while ans==[]:
        print "Searching degree",d,"covers."
        ans = look_for_good_cover(wordlist,rank,d)
        d=d+1
