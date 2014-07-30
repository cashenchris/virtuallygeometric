from group import *
import copy
import random



class WhiteheadAuto(Automorphism):
    """
    Whitehead automorphism defined by generator or inverse generator x and set or list Z of generators and inverses including x and not including X defined by:
    x -> x
    y -> xy if y in Z and Y not in Z
    y -> xyX if y and Y both in Z
    y -> y if neither y nor Y in Z
    """
    def __init__(self,F,xinput,Zinput):
        assert(F.rank)
        self.codomain=F
        self.domain=F
        if xinput or Zinput:
            self.x=xinput 
            self.Z=set(Zinput)
        else: # identity automorphism
            self.x=1
            self.Z=set([1])
        assert((self.x in self.Z) and not(-self.x in self.Z))          

    def variant_generators(self):
        return [i for i in self.Z if i!=self.x]
        
    def __call__(self,w): #evaluate the automorphism on the word w and return a word in codomain
        theletters=copy.copy(w.letters)
        imagewordletters=[]
        while theletters:
            nextletter=theletters.pop(0)
            if abs(nextletter)==abs(self.x):
                nextwordletters=[nextletter]
            elif nextletter in self.Z and -nextletter not in self.Z:
                nextwordletters=[self.x,nextletter]
            elif nextletter not in self.Z and -nextletter in self.Z:
                nextwordletters=[nextletter,-self.x]
            elif nextletter in self.Z and -nextletter in self.Z:
                nextwordletters=[self.x,nextletter,-self.x]
            else:
                nextwordletters=[nextletter]
            imagewordletters+=nextwordletters
        return self.codomain.word(imagewordletters)

    def inverse(self):
        if self.x and self.Z:
            newZ=copy.copy(self.Z)
            newZ.remove(self.x)
            newZ.add(-self.x)
            return WhiteheadAuto(self.domain,-self.x,newZ)
        else:
            return self

    def __mul__(self,other): # In some special cases the result is still a Whithead automorphism
        if type(other)==type(self):
            if self.x==other.x and not self.Z&other.Z:
                return WhiteheadAuto(self.domain,self.x,self.Z|other.Z)
            elif self.x==-other.x and not other.Z-self.Z:
                return WhiteheadAuto(self.domain,self.x,self.Z-other.Z)
            elif self.x==-other.x and not self.Z-other.Z:
                return WhiteheadAuto(self.domain,-self.x,other.Z-self.Z)
            else: # otherwise just multiply them as automorphisms
                return Automorphism.__mul__(self,other)
        else: # otherwise just multiply them as automorphisms
            return Automorphism.__mul__(self,other)
            
            

    def __pow__(self,n):
        result=WhiteheadAuto(self.domain,0,[]) # the identity automorphism
        # If n=0,1,-1 then result is still a Whitehead automorphism
        if n==0:
            return result
        elif n==1:
            return self
        elif n==-1:
            return self.inverse()
        elif(n)>1:
            for i in range(n):
                result=result * self
            return result
        else:
            inverse=self.inverse()
            for i in range(-n):
                result=result * inverse
            return result
            

def random_whitehead_automorphism(F):
    """
    Generate a random Whitehead automorphism of a free group F.
    """
    vertices=range(-F.rank,0)+range(1,F.rank+1)
    x=random.choice(vertices)
    Z=[x]
    vertices.remove(x)
    vertices.remove(-x)
    for v in vertices:
        if random.random()<.5:
            Z.append(v)
    return WhiteheadAuto(F,x,Z)


def random_automorphism_pair(F,length):
    """
    Generate an automorphism and its inverse by taking a product of 'length' random Whitehead automorphisms.
    """
    randomaut=Automorphism(F)
    inverse=Automorphism(F)
    for i in range(length):
        w=random_whitehead_automorphism(F)
        randomaut=w*randomaut
        inverse=inverse*(w.inverse())
    return randomaut, inverse

def random_automorphism(F, length):
    return random_automorphism_pair(F,length)[0]



def is_inner_auto_by(alpha):
    """
    Return a word w such that alpha(f)=w**(-1)*f*w for all f in alpha.domain.
    """
    F=alpha.domain
    assert(F.isFree())
    if F.rank==1:
        if alpha(F.word([1]))==F.word([1]):
            return F.word([])
        else:
            return None
    else:
        conjugators=[F.get_conjugator(F.word([i]),alpha(F.word([i]))) for i in range(1,1+F.rank)]
        if any(conj is None for conj in conjugators):
            return None
        # alpha(F.word([i]))=conjugators[i-1]**(-1)*F.word([i])*conjugators[i-1]
        else:
            # if alpha(z)=w**(-1)*z*w for all z in F
            # then for each i there is an a[i-1] so that F.word([i])**(a[i-1])*conjugators[i-1]=w
            a=[None for i in range(F.rank)]
            x=(conjugators[0])*(conjugators[1])**(-1)
            # if inner, x must be of the form F.word([2])**(a[1])*F.word([1])**(-a[0])
            if len(x)==0:
                a[0]=0
                a[1]=0
            else:
                if x.letters[0]==2 and x.letters[-1]==1:
                    change=x.letters.index(1)
                    if x.letters==[2 for i in range(change)]+[1 for i in range(len(x.letters)-change)]:
                        a[1]=change-1
                        a[0]=-(len(x.letters)-change)
                    else:
                        return None
                elif x.letters[0]==-2 and x.letters[-1]==1:
                    change=x.letters.index(1)
                    if x.letters==[-2 for i in range(change)]+[1 for i in range(len(x.letters)-change)]:
                        a[1]=-(change-1)
                        a[0]=-(len(x.letters)-change)
                    else:
                        return None
                elif x.letters[0]==-2 and x.letters[-1]==-1:
                    change=x.letters.index(1)
                    if x.letters==[-2 for i in range(change)]+[-1 for i in range(len(x.letters)-change)]:
                        a[1]=-(change-1)
                        a[0]=(len(x.letters)-change)
                    else:
                        return None
                elif x.letters[0]==2 and x.letters[-1]==1:
                    change=x.letters.index(1)
                    if x.letters==[-2 for i in range(change)]+[1 for i in range(len(x.letters)-change)]:
                        a[1]=(change-1)
                        a[0]=-(len(x.letters)-change)
                    else:
                        return None
                elif abs(x.letters[0])==2:
                    if x.letters==[x.letters[0] for i in range(len(x.letters))]:
                        a[1]=sign(x.letters[0])*len(x.letters)
                        a[0]=0
                    else:
                        return None
                elif abs(x.letters[0])==1:
                    if x.letters==[x.letters[0] for i in range(len(x.letters))]:
                        a[1]=0
                        a[0]=-sign(x.letters[0])
                else:
                    return None
            w=F.word([1])**(a[0])*conjugators[0]
            # on the first two generators alpha is conjugation by w
            for i in range(3,1+F.rank):
                # if inner then w=F.word([i])**(a[i-1])*conjugators[i-1] for some a[i-1]
                x=w*(conjugators[i-1])**(-1)
                if len(x)==0:
                    a[i-1]=0
                    i+=1
                else:
                    if abs(x.letters[0])==i:
                        if x.letters==[x.letters[0]]*len(x.letters):
                            a[i-1]=sign(x.letters[0])*len(x.letters)
                            i+=1
                        else:
                            return None
                    else:
                        return None
            return w
                
def is_inner_auto(alpha):
    w=is_inner_auto_by(alpha)
    if w is None:
        return False
    else:
        return True
