import freegroup
import group
import itertools
import whiteheadgraph.build.whiteheadreduce as wr

def advance_counter(thecounter,place,resetval):
    """
    Advance a little endian counter. thecounter is a list of non-negative integers. place is index to be incremented by 1. resetval is number at which the place should rollover and the next place should be incremented.
    Return value is 0 if the result has the same number of places as the original counter, or 1 if the length of the counter has increased.
    """
    thecounter[place]=(thecounter[place]+1)%resetval
    if thecounter[place]==0:
        if place+1 not in thecounter:
            thecounter[place+1]=0
            return 1
        else:
            result=advance_counter(thecounter,place+1,resetval)
            return result
    else:
        return 0
    

def generate_words(F,maxlen,startlength=1):
    """
    Generator of unique, non-trivial words in a free group F up to length maxlen.
    """
    rank=F.rank
    letters=[x for x in range(1,1+rank)]+[-x for x in range(1,1+rank)]
    counters=dict()
    for i in range(startlength):
        counters[i]=0
    maxindex=startlength-1
    theletters=[letters[counters[i]] for i in range(maxindex+1)]
    theword=F.word(theletters)
    while maxindex<maxlen:    
        yield theword
        maxindex=maxindex+advance_counter(counters,0,2*rank)
        theletters=[letters[counters[i]] for i in range(maxindex+1)]
        theword=F.word(theletters)
        while len(theword)<1+maxindex:# while there is some free reduction in the word, discard it and get a new one
            maxindex=maxindex+advance_counter(counters,0,2*rank)
            theletters=[letters[counters[i]] for i in range(maxindex+1)]
            theword=F.word(theletters)

def generate_words_in_commutator_subgroup(F,maxcommutatorlength,maxsinglewordlength):
    """
    Generator of words in commutator subgroup of F.
    """
    wordgen=generate_words(F,maxsinglewordlength)
    while True:
        w=wordgen.next()
        if freegroup.ishomologicallytrivial(w):
            yield w
    #wordgen=generatewords(F,maxsinglewordlength)
    #multiwordgen=itertools.product(wordgen,repeat=2*maxcommutatorlength)
    #while True:
    #    multiword=multiwordgen.next()
    #    pairlist=[(multiword[2*i],multiword[2*i+1]) for i in range(maxcommutatorlength)]
    #    commutatorlist=[p[0]*p[1]*p[0]**(-1)*p[1]**(-1) for p in pairlist]
    #    commutator=group.product(*commutatorlist)       
    #    if len(commutator)>0:
    #       yield commutator

def is_minimally_ordered(thetuple):
    if len(thetuple)>1:
        if len(thetuple)%2>0:
            return False
    maxindices=[]
    for i in range(1,len(thetuple)):
        if abs(thetuple[i])>thetuple[0]:
            return False
        elif abs(thetuple[i])==thetuple[0]:
            maxindices.append(i)
    while maxindices:
        theindex=maxindices.pop()
        newword=thetuple[theindex:]+thetuple[:theindex]
        if newword[0]<0:
            newnewword=tuple((-1)**(i+1)*newword[i] for i in range(len(newword)))
        else:
            newnewword=newword
        if len(newnewword)>1:
            if newnewword[1]<0:
                newnewnewword=tuple((-1)**(i+1)*newnewword[i] for i in range(len(newnewword)))
            else:
                newnewnewword=newnewword
        else:
            newnewnewword=newnewword
        if thetuple<newnewnewword:
            return False
    return True
        
            
            

def aut2forbits(F,maxlength):
    reps=[]
    lastroundadded=[]
    newwords=[]
    for i in range(1,1+maxlength):
        reps.append((i,))
        lastroundadded.append((i,))
    while lastroundadded:
        oldword=lastroundadded.pop()
        wordlength=sum(abs(l) for l in oldword)
        if wordlength<maxlength:
            for k in range(1,1+min(maxlength-wordlength,oldword[0])):
                newword=oldword+(k,)
                if is_minimally_ordered(newword):
                    newwords.append(newword)
                    reps.append(newword)
    lastroundadded=newwords
    newwords=[]
    while newwords+lastroundadded:
        while lastroundadded:
            oldword=lastroundadded.pop()
            wordlength=sum(abs(l) for l in oldword)
            if wordlength<maxlength:
                for k in range(1,1+min(oldword[0],maxlength-wordlength))+range(max(-oldword[0],-maxlength+wordlength+1),0):
                    for h in range(1,1+maxlength-wordlength-k)+range(-1-maxlength+wordlength+k,-1):
                        newword=oldword+(k,h)
                        if is_minimally_ordered(newword):
                            newwords.append(newword)
                            reps.append(newword)
        lastroundadded=newwords
        newwords=[]
    wordreps=[]
    for rep in reps:
        theletters=[]
        for i in range(len(rep)):
            if i%2:
                if rep[i]<0:
                    theletters.extend([-2]*abs(rep[i]))
                else:
                    theletters.extend([2]*rep[i])
            else:
                if rep[i]<0:
                    theletters.extend([-1]*abs(rep[i]))
                else:
                    theletters.extend([1]*rep[i])
        theword=F.word(theletters)
        theroot=F.max_root(theword,withpower=False)
        if len(theroot)>1:
            minimal=wr.whitehead_minimal(F,[theroot],simplified=True,blind=True)
            if len(minimal['wordlist'][0])<len(theroot):
                pass
            else:
                wordreps.append(theword)
        else:
            wordreps.append(theword)
    return wordreps
