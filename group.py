import copy
import random
from numpy import log, ceil, sqrt, sign
     
 
class FGGroup(object):
    """
    A finitely generated group with a fixed generating set. Internally the generating set is integers 1..n.
    gens = optional list of strings containing names of the generators, ie ['x','y','z']
    optional keyword arguments:
    inverses = list of strings containing names of the inverses of gens, ie inverses=['X','Y','Z']
        If gens are supplied but no inverses the inverses will be named like ['x**(-1)', 'y**(-1)', 'z**(-1)']
        This is ignored if no gens are supplied. 
    numgens =  number of generators, to be given generic names. This is ignored if gens is supplied.
    identity = string representing the identity element, ie identity='' or identity='e'. Defaults to ''.
    generatorbasename = prefix to give generator names if no gens is supplied, ie gens=[] and generatorbasename='a' will give generators 'a1', 'a2',...
    displaystyle =  choice of list or string style for displaying group elements.
        In the above examples, with inveses supplied and displaystyle=str the commutator word [1,2,-1,-2] would print 'xyXY'.
        With inverses not supplied and displaystyle=list it would print ['x','y','x**(-1)','y**(-1)'].
        Default:  str if gens=[] and numgens <= 26 else list
    
    """
    # Note at this point we don't do anything with relations, so this is really just the free group on the given generators.

    def __init__(self,gens=[],**kwargs):
        try:
            self.identity=kwargs['identity']
        except KeyError:
            self.identity=''
        if not gens: # no generator names specified
            try:
                numgens=kwargs['numgens']
            except KeyError:
                numgens=0
                
            if numgens==0: # group is the trivial group
                self.gens=[]
                self.lettering=[self.identity]
                self.displaystyle=list
            elif numgens<0:
                raise ValueError("numgens must be non-negative")
            else:
                try:
                    generatorbasename=kwargs['generatorbasename']
                except KeyError:
                    generatorbasename=None
                if numgens<=26 and generatorbasename is None: # default to alphabetic representation
                    self.gens=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'][:numgens]
                    self.lettering=['Z','Y','X','W','V','U','T','S','R','Q','P','O','N','M','L','K','J','I','H','G','F','E','D','C','B','A',self.identity,'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'][26-numgens:27+numgens]
                    self.displaystyle=str
                else:
                    self.displaystyle=list
                    if generatorbasename is None:
                        generatorbasename='g'
                    self.gens=[generatorbasename+'_'+str(i) for i in range(1,numgens+1)]
                    self.lettering=['('+generatorbasename+'_'+str(i)+')**(-1)' for i in reversed(range(1,numgens+1))]+[self.identity]+[generatorbasename+'_'+str(i) for i in range(1,numgens+1)]
        else: # generator names were given
            self.gens=gens
            try:
                inverses=kwargs['inverses']
            except KeyError:
                inverses=[]
            if inverses:
                self.lettering=[x for x in reversed(inverses)]+[self.identity]+gens
                self.displaystyle=str
            else:
                self.lettering=['('+x+')**(-1)' for x in reversed(gens)]+[self.identity]+gens
                self.displaystyle=list
        try: # Check to see if automatically set self.displaystyle is overridden by an argument
            self.displaystyle=kwargs['displaystyle']
        except KeyError:
            pass
                

    def __repr__(self):
        return "< "+", ".join([str(g) for g in self.gens])+" >"

    def generators(self):
        return self.gens

    def free_reduce(self,w):
        return self.word(freereduce(w.letters))

    # comparison operators for words in the group, to be overriden in subclasses if we know a good way to compare
    def word__ne__(self,v, w):
        assert(False) # do not compare words in a generic group
        return v.letters!=w.letters
        
    def word__eq__(self,v, w):
        assert(False) # do not compare words in a generic group
        return v.letters==w.letters

    def word__cmp__(self,v, w):
        assert(False) # do not compare words in a generic group
        return v.letters<w.letters

    def word__hash__(self,w):
        return None
       

    def cyclic_reduce(self,w):
        w1=copy.copy(w.letters)
        w1=freereduce(w1)
        while len(w1) > 2 and w1[0]+w1[-1]==0:
            w1=w1[1:-1]
        return self.word(w1)

    def cyclic_reducer(self,w):
        """
        return w0,w1 such that w1 is cyclically reduced and w0**(-1)w1w0=wi.
        """
        w1=copy.copy(w.letters)
        w0=[]
        while len(w1) > 2 and w1[0]+w1[-1]==0:
            w0.insert(0,w1[-1])
            w1=w1[1:-1]
        return self.word(w0), self.word(w1)

    def word(self,letters):
        """
        Construct a word in the group.
        """
        return Word(letters, self)

    def is_subgroup(self,G):
        """
        Check if group is a subgroup of G.
        """
        # No group theory here, subgroups are explicitly declared. This just looks for a chain of declarations leading to G.
        currentgroup=self
        while currentgroup!=G:
             try:
                  currentgroup=currentgroup.supergroup
             except AttributeError:
                  return False
        return True
     
    def get_inclusion(self,G):
        """
        If self is a subgroup of G return the inclusion homomorphism, else error.
        """
        # As in isSubgroup, nothing intelligent here. Subgroups are defined with an inclusion into a supergroup.
        currentgroup=self
        inclusionchain=Automorphism(self)
        while currentgroup is not G:
            try:
                inclusionchain=compose(currentgroup.inclusion,inclusionchain)
                currentgroup=currentgroup.supergroup
            except AttributeError:
                 raise ValueError("not a known supergroup")
        return inclusionchain

    def randomword(self,length):
        """
        Word that is the result of a random walk without backtracking of given length in the generators and inverses.
        """
        numgens=len(self.gens)
        letterlist = range(1,numgens+1)+range(-1,-(numgens+1),-1)
        letters = []
        for n in range(length):
            nextletter=random.choice(letterlist)
            if len(letters):
                while nextletter==-letters[-1]:
                    nextletter=random.choice(letterlist)
            letters.append(nextletter)
        return self.word(letters)

    def random_word(self,length):
        """
        same as randomword
        """
        return self.randomword(length)

    def random_cyclically_reduced_word(self,length):
        """
        Word that is the result of a random walk without backtracking of given length in the generators and inverses, and such that last letter is not inverse of first letter.
        """
        numgens=len(self.gens)
        letterlist = range(1,numgens+1)+range(-1,-(numgens+1),-1)
        letters = []
        if length==0:
            pass
        elif length==1:
            nextletter=random.choice(letterlist)
            letters.append(nextletter)
        else:
            for n in range(length-1):
                nextletter=random.choice(letterlist)
                if len(letters):
                    while nextletter==-letters[-1]:
                        nextletter=random.choice(letterlist)
                letters.append(nextletter)
            nextletter=random.choice(letterlist)
            while nextletter==-letters[-1] or nextletter==-letters[0]:
                nextletter=random.choice(letterlist)
            letters.append(nextletter)
        return self.word(letters)

    def random_walk(self,length):
        """
        Word that is the free reduction of a random walk of length n in the generators and inverses.
        """
        numgens=len(self.gens)
        letterlist = range(1,numgens+1)+range(-1,-(numgens+1),-1)
        letters = []
        for n in range(length):
            letters.append(random.choice(letterlist))
        return self.word(letters)
        

    def random_multiword(self,length):
        """
        Generate a list of random words of random length with total characters of about length.
        """
        if length==0:
            return [self.word([])]
        numgens=len(self.gens)
        wordlist=[]
        totallength=0
        numberofwords=random.choice(range(1,max(2,1+int(sqrt(length)))))
        for i in range(1,numberofwords+1):
            if totallength<length:
                ilen=random.choice(range(1,1+length-totallength))
                w=self.randomword(ilen)
                wordlist.append(w)
                totallength+=len(w.letters)
            else:
                wordlist.append(self.word([]))
        return wordlist
      
                  
    
class FPGroup(FGGroup):
    """
    A finitely presented group.
    """
    def __init__(self,gens=[], rels=[], **kwargs):
        FGGroup.__init__(self,gens, **kwargs)
        self.rels=rels

    def __repr__(self):
        return "< "+", ".join([str(g) for g in self.gens])+" | "+", ".join([self.word(r)() for r in self.rels])+" >"

    def relators(self):
        return self.rels
    





class FGSubgroup(FGGroup):
    """
    A subgroup defined by an inclusion homomorphism into a supergroup, with optional names for the generators of the subgroup.
    """
    def __init__(self,supergroup, inclusionlist, **kwargs):
        self.supergroup=supergroup
        # To initialize the subgroup the default arguments become those of the supergroup
        try:
            idelt=kwargs['identity']
        except KeyError:
            idelt=supergroup.identity
        try:
            dstyle=kwargs['displaystyle']
        except KeyError:
            dstyle=supergroup.displaystyle
        try:
            generators=kwargs['gens']
        except KeyError:
            generators=None
        try:
            invs=kwargs['inverses']
        except KeyError:
            invs=None
        try:
            gbn=kwargs['generatorbasename']
        except KeyError:
            gbn=None
            
        if (not generators) and (not gbn): # if we haven't said how to name subgroup elements then we'll just call them as we do in the supergroup
            generators=[w() for w in inclusionlist]
        
        FGGroup.__init__(self, generators, inverses=invs, generatorbasename=gbn, numgens=len(inclusionlist), displaystyle=dstyle, identity=idelt)
        self.inclusion=Homomorphism(self, supergroup, dict([(i,inclusionlist[i-1]) for i in range(1,1+len(inclusionlist))]))         

    def __repr__(self):
        if self.gens==[]: #trivial group
            return "{"+self.identity+"}"
        elif str(self.word([1])(self))==str(self.word([1])(self.supergroup)): # check if the subgroup generators were given new names
            return "< "+", ".join([str(self.word([i])(self)) for i in range(1,1+len(self.gens))])+" >" # if not, just output names of the generators
        else:
            return "< "+", ".join([str(self.word([i])(self))+"="+str(self.word([i])(self.supergroup)) for i in range(1,1+len(self.gens))])+" >" # if so, output new name for generator = word in supergroup

    def __call__(self,ancestor=None):
        """
        return string < w1, ... > where g1 is a generator of the subgroup and w1 is inclusion of g1 into the supergroup ancestor.
        """
        if ancestor is None:
            ancestor=self.supergroup
        G=ancestor
        selfinG=self.get_inclusion(G)
        return "< "+", ".join([selfinG(self.word([i]))() for i in range(1,1+len(self.gens)) ])+" >"
        
        








def common_ancestor(G,H):
        """
        Return a (smallish) group containing both G and H.
        """
        Gancestors=[G]
        Hancestors=[H]
        while all([Gancestors[-1] is not y for y in Hancestors]) and all([Gancestors[-1] is not y for y in Hancestors]):
            if (not hasattr(Gancestors[-1],'supergroup')) and (not hasattr(Hancestors[-1],'supergroup')):
                return None
            else:
                try:
                    Gancestors.append(Gancestors[-1].supergroup)
                except AttributeError:
                    pass
                try:
                    Hancestors.append(Hancestors[-1].supergroup)
                except AttributeError:
                    pass
                
        if any([Gancestors[-1] is y for y in Hancestors]):
            return Gancestors[-1]
        elif any([Hancestors[-1] is y for y in Gancestors]):
            return Hancestors[-1]
        else:
            raise RuntimeError()
     

#----------  words
    
def freereduce(l):
    """
    Free reduces a list of integers.
    """
    if len(l)<2:
        return l
    else:
        reduction=[l[0]]
        for i in range(len(l)-1):
            if reduction==[]:
                reduction=[l[i+1]]
            else:
                if reduction[-1]==-l[i+1]:
                    reduction.pop()
                else:
                    reduction.append(l[i+1])
    return reduction

def stringtolist(letters, lettering):
    """
    Convert some alphabetic word to a list of numbers.
    """
    numbers=[]
    prefix=''
    suffix=copy.copy(letters)
    while suffix:
        prefix+=suffix[0]
        suffix=suffix[1:]
        if prefix in lettering:
             numbers.append(lettering.index(prefix)-len(lettering)//2)
             prefix=''
    if prefix!='':
         raise ValueError(prefix+" was not found in "+str(lettering))
    return numbers
         


class Word(object):
    """
    Word in a group.
    The word itself is just a list of nonzero integers. Expressing the word as a group element is handled by the group.
    Letters can be either a list of nonzero integers or a list or string consisting of generators or inverse generators of the group.

        F=FGFreeGroup(gens=['x','y'], inverses=['X','Y'])
        w1=Word([1,2,-1,-2],F)
        w2=Word('xyXY',F)
        
        produce the same word, as does
        F.word([1,2,-1,-2]), or
        F.word('xyXY')
        
    Calling a word gives its representation in the group.
    w1() = 'xyXY'

    Calling a word with argument a supergroup of the defining group gives the representation of the word in the supergroup.
        G=FGFreeSubgroup([w1], F, gens=[a], ivnerses=['A']) is the cyclic subgroup of F generated by a=w1.
        w3=G.word([1,1,1])
        w3()='aaa'
        w3(G)='aaa'
        w3(F)='xyXYxyXYxyXY'
    """   
    def __init__(self,letters, group):
        self.group=group
        if hasattr(letters,'group'):# letters is already a word in some group
            self.letters=[x for x in letters.letters]
        elif hasattr(letters,'isalpha'): #letters is a string
            self.letters=freereduce(stringtolist(letters,self.group.lettering))
        else:
            if len(letters)==0:
                self.letters=[]
            elif letters==['']:
                self.letters=[]
            elif hasattr(letters[0],'isalpha'):
                self.letters=freereduce([group.lettering.index(x)-(len(group.lettering)//2) for x in letters])
            else:
                self.letters=freereduce([x for x in letters])


    def __len__(self):
        return len(self.letters)

    def __repr__(self):
        return str(self.letters)      

    def __mul__(self,other):
        G=common_ancestor(self.group,other.group)
        if G is None:
            raise TypeError("can not multiply words that are not in a common group")
        else:
            selfinG=self.group.get_inclusion(G)(self)
            otherinG=other.group.get_inclusion(G)(other)
            return G.word(selfinG.letters+otherinG.letters)
    
    def __ne__(self, other):
        G=common_ancestor(self.group,other.group)
        if G is None:
            return True
        else:
            return G.word__ne__(self.group.get_inclusion(G)(self),other.group.get_inclusion(G)(other))

    def __eq__(self, other):
        G=common_ancestor(self.group,other.group)
        if G is None:
            return False
        else:
            return G.word__eq__(self.group.get_inclusion(G)(self),other.group.get_inclusion(G)(other))

    def __cmp__(self, other):
        G=common_ancestor(self.group,other.group)
        if G is None:
            return False
        else:
            return G.word__cmp__(self.group.get_inclusion(G)(self),other.group.get_inclusion(G)(other))

    def __hash__(self):
        G=self.group
        while hasattr(G,'supergroup'):
            G=G.supergroup
        return G.word__hash__(self.group.get_inclusion(G)(self))
        


    def __pow__(self,n):
        # take n'th power
        result = Word([], self.group)
        if n == 0:
            return result
        elif n>0:
            for i in range(n):
                result = result * self
            return result
        else:
            inverse = Word([-i for i in reversed(self.letters)], self.group)
            for i in range(-n):
                result = result * inverse
            return result
        
    def is_element(self,G):
        """
        True if word belongs to a subgroup of G.
        """
        return self.group.is_subgroup(G)

    def __call__(self,supergroup=None):
        """
        Write out the word in terms of the generators of a supergroup.
        """
        if supergroup is None:
             supergroup=self.group
        inclusionchain=Automorphism(self.group) # start with the identity automorphism
        currentgroup=self.group
        while currentgroup!=supergroup:
            try:
                inclusionchain=compose(currentgroup.inclusion,inclusionchain)
                currentgroup=currentgroup.supergroup
            except AttributeError:
                print str(self)+" is not a recognized word in "+supergroup
        if len(self.letters)==0:
            return supergroup.identity
        else:
            imageword=inclusionchain(self)
            if supergroup.displaystyle==str:
                answer = ''
                for letter in imageword.letters:
                    answer+=supergroup.lettering[letter+len(supergroup.lettering)//2]
            else:
                answer=[]
                for letter in imageword.letters:
                    answer.append(supergroup.lettering[letter+len(supergroup.lettering)//2])
            return str(answer)

    def cycle(self,steps=1):
        """
        Return cyclic permutation by some number of steps. abc -> bca
        """
        cycledletters=copy.copy(self.letters)
        if len(self)>0:
            remainingsteps=steps%len(self)
            while remainingsteps:
                cycledletters=cycledletters[1:]+cycledletters[:1]
                remainingsteps-=1
        return self.group.word(cycledletters)
                
    
                  
         
       
        
    def pop(self):
        """
        Return first letter (as a number!), and shorten word.
        """
        first = self.letters[0]
        self.letters = self.letters[1:]
        return first

    def alpha(self):
        """
        Return word as string.
        """
        key = 'ZYXWVUTSRQPONMLKJIHGFEDCBA abcdefghijklmnopqrstuvwxyz'
        key_offset = 26
        strout = ''
        for letter in self.letters:
            strout = strout+key[letter+key_offset]
        return strout



def guess_rank(*wordlist):
    """
    The highest generator appearing in wordlist.
    """
    rank=0
    for w in wordlist:
        if w.letters!=[]:
            rank=max(rank,max(abs(i) for i in w.letters))
    return rank



#---------- homomorphisms

class PDHomo(object):
    """
    Partially defined group homomorphism. May only be defined on some of the generators of the domain.
    """
    def __init__(self,domain, codomain, generatorimagedict=None): 
        # listofimagesofgens[i]=word in codomain that is image of i st generator of domain
        self.domain=domain
        self.codomain=codomain
        if generatorimagedict is None:
             self.images=dict()
        else:
             self.images=generatorimagedict

    def variant_generators(self):
        return self.images.keys()
    
    def __repr__(self):
        return str(self.domain)+" -> "+str(self.codomain)+":"+str(self.images)
    
    def __str__(self):
        imagelist=[]
        vG=list(set([abs(j) for j in self.variant_generators()]))
        vG.sort()
        for i in vG:
            imagelist+=[self.domain.word([i])()+" -> "+self(self.domain.word([i]))()]
        if len(vG)<len(self.domain.gens):
            imagelist+=["else undefined"]
        return str(self.domain)+" -> "+str(self.codomain)+":\n"+"\n".join(imagelist)

    def __call__(self,w): #evaluate the homomorphism on the word w and return a word in codomain
        theletters=copy.copy(w.letters)
        imagewordletters=[]
        while theletters:
            nextletter=theletters.pop(0)
            try:
                nextwordletters=self.images[nextletter].letters
            except KeyError:
                try:
                    nextwordletters=[-i for i in reversed(self.images[-nextletter])]
                except KeyError:
                    raise KeyError('The map is not defined on generator '+str(nextletter)+' of the domain.')
            imagewordletters+=nextwordletters
        return self.codomain.word(imagewordletters)

    def alpha(self):
        print str(domain)+" -> "+str(codomain)+":\n"+"\n".join([domain.word([i]).alpha()+" -> "+self(self.domain.word([i])).alpha() for i in self.variant_generators()])

class Homomorphism(PDHomo):
    """
    Homomorphism between groups. Generators not  in the generatorimagedict are sent to the trivial word.
    """
    def __init__(self,domain,codomain,generatorimagedict=None):
        PDHomo.__init__(self,domain,codomain,generatorimagedict)
        
    def __call__(self,w): #evaluate the homomorphism on the word w and return a word in codomain
        theletters=copy.copy(w.letters)
        imagewordletters=[]
        while theletters:
            nextletter=theletters.pop(0)
            try:
                nextwordletters=self.images[nextletter].letters
            except KeyError:
                try:
                    nextwordletters=[-i for i in reversed(self.images[-nextletter].letters)]
                except KeyError:
                    nextwordletters=[] # generators not in generatorimagedict are sent to trivial word
            imagewordletters+=nextwordletters
        return self.codomain.word(imagewordletters)

    def __str__(self):
        imagelist=[]
        vG=list(set([abs(j) for j in self.variant_generators()]))
        vG.sort()
        for i in vG:
            imagelist+=[self.domain.word([i])()+" -> "+self(self.domain.word([i]))()]
        if len(vG)<len(self.domain.gens):
            imagelist+=["else"+" -> "+str(self.codomain.word([]))]
        return str(self.domain)+" -> "+str(self.codomain)+":\n"+"\n".join(imagelist)


class PDEndo(PDHomo):
    """
    Partially defined endomorphism of a group.
    """
    def __init__(self, domain, generatorimagedict=None):
        PDHomo.__init__(self, domain, domain, generatorimagedict)

class Endomorphism(PDEndo, Homomorphism):
    """
    Endomorphism of a group. Generators not in the generatorimagedict are sent to the trivial element.
    """
    def __init__(self, domain, generatorimagedict=None):
        Homomorphism.__init__(self, domain, domain, generatorimagedict)

class Automorphism(Endomorphism):
    """
    Automorphism of a group. Generators not in the generatorimagedict are assumed to be fixed.
    """
    # Note, the class doesn't check that the input actually defines an automorphism.
    def __init__(self, domain, generatorimagedict=None):
        Endomorphism.__init__(self,domain,generatorimagedict)
    
    def __mul__(self,other):
        assert(self.domain is other.domain)
        return Automorphism(self.domain, dict([(i,self(other(other.domain.word([i])))) for i in set(other.variant_generators())|set(self.variant_generators())]))

    def __str__(self):
        imagelist=[]
        vG=list(set([abs(j) for j in self.variant_generators()]))
        vG.sort()
        for i in vG:
            imagelist+=[self.domain.word([i])()+" -> "+self(self.domain.word([i]))()]
        if len(vG)<len(self.domain.gens):
            imagelist+=["else invariant"]
        return str(self.domain)+" -> "+str(self.codomain)+":\n"+"\n".join(imagelist)

    def inverse(self):
         raise RuntimeError("inverse not implemented for arbitrary automorphisms")

    def __pow__(self,n):
        # for large n would recursive definition be faster or slower?
        result=Automorphism(self.domain) # the identity automorphism
        if(n)>=0:
            for i in range(n):
                result=result * self
        else:
            inverse=self.inverse()
            for i in range(-n):
                result=result * inverse
        return result
        
    def __call__(self,w): #evaluate the automorphism on the word w and return a word in codomain
        theletters=copy.copy(w.letters)
        imagewordletters=[]
        while theletters:
            nextletter=theletters.pop(0)
            try:
                nextwordletters=self.images[nextletter].letters
            except KeyError:
                try:
                    nextwordletters=[-i for i in reversed(self.images[-nextletter].letters)]
                except KeyError:
                    nextwordletters=[nextletter] # if neither the generator nor its inverse are in the dict then it is fixed by the automorphism
            imagewordletters+=nextwordletters
        return self.codomain.word(imagewordletters)

class InnerAutomorphism(Automorphism):
    """
    Inner automorphism.
    """
    def __init__(self, domain, gpelement):
        self.gpelement=gpelement
        self.domain=domain
        self.codomain=domain

    def variant_generators(self):
        return range(1,1+len(self.domain.gens))
        
    def __call__(self,w):
        return self.codomain.word([-i for i in reversed(gpelement.letters)]+w.letters+gpelement.letters)
        
    def __mul__(self,other):
        if type(self)==type(other) and self.domain is other.domain:
            return InnerAutomorphism(self.domain,self.gpelement*other.gpelement)
        else:
            return Automorphism.__mul__(self,other)

    def __pow__(self,n):
        return InnerAutomorphism(self.domain, self.gpelement**(n))

    def inverse(self):
        return self**(-1)


def PDcompose(alpha,beta):
    return PDHomo(beta.domain,alpha.codomain,dict([(i,alpha(beta(beta.domain.word([i])))) for i in beta.variant_generators()]))


def compose(alpha,beta):
    """
    Compose two homomorphisms.
    """
    # special cases because an non-variant generator for an autmorphism maps to itself, while for a homomorphism it maps to the identity element
    if isinstance(beta,Automorphism) and isinstance(alpha,Automorphism):
        return alpha*beta
    elif isinstance(beta,Automorphism):
        return Homomorphism(beta.domain,alpha.codomain,dict([(i,alpha(beta(beta.domain.word([i])))) for i in set(beta.variant_generators()+alpha.variant_generators())]))
    else:
        return Homomorphism(beta.domain,alpha.codomain,dict([(i,alpha(beta(beta.domain.word([i])))) for i in beta.variant_generators()]))


def product(*args):
    def product2(alpha,beta):
        return alpha*beta
    if len(args)>1:
        return reduce(product2,args)
    elif len(args)>0:
        return args[0]
    else:
        raise TypeError("Don't know how to take product of 0 things.")
