from __future__ import division
# from Numeric import sign  
import random
import copy

key = 'ZYXWVUTSRQPONMLKJIHGFEDCBA abcdefghijklmnopqrstuvwxyz'
key_offset = 26

def sign(n):
    # Returns 0 for 0, -1 if <0, +1 if >0
    if n==0:
        return n
    return (-1)**(n<0)

def freereduce(letters,n=0):
    positions = range(n,len(letters)-1)+range(0,n)
    for pos in positions:
        # If there is a reduction, recurse.
        if letters[pos]+letters[pos+1] == 0:
            return freereduce(letters[:pos]+letters[pos+2:],max([pos-1,0]))
    #If there were no reductions to do, we return the list given us.
    return letters

def stringtolist(letters):
    """
    Convert some alphabetic word to a list of numbers.
    """
    if not letters.isalpha():
        return []
    return [key.find(letters[0])-key_offset]+stringtolist(letters[1:])

    

class word(object):
    """
    Word in the free group, together with function for evaluating word in
    some representation
    """
    # Primary representation is as a list of integers
    def __init__(self,letters):
        """
        Constructor takes a list of integers or string.
        
        Example: w = word([1,2,-1,-2]) followed by w([a,b]) returns
        the commutator abAB.

        Words can also be specified as strings
        of lower and upper case letters, e.g. word('abAB') returns the
        same thing as word([1,2,-1,-2]).
        """
        if type(letters)==type([]):
            self.letters=freereduce(letters)
        elif type(letters)==type('') and letters.isalpha():
            self.letters=stringtolist(letters)
        else:
            print('WARNING:List of integers or alphabetic string expected, returning [].\n')
            self.letters=[]

    def __len__(self):
        return len(self.letters)

    def __repr__(self):
        return str(self.letters)

    def __mul__(self,other):
        return word(self.letters+other.letters)

    def __cmp__(self,other):
        " shortlex order "
        if self.letters == other.letters:
            return 0
        elif len(self.letters)!=len(other.letters):
            return len(self.letters)-len(other.letters)
        else:
            # This is probably too clever.  Also it might be better to order differently.
            return -2*int(self.letters<other.letters)+1

    def __pow__(self,n):
        # take n'th power
        result = word([])
        if n == 0:
            return result
        elif n>0:
            for i in range(n):
                result = result * self
            return result
        else:
            inverse = word(map(lambda x: -x, self.letters[::-1]))
            for i in range(-n):
                result = result * inverse
            return result

    def __call__(self,arg,one='Identity?'):
# The optional argument is so that an identity element of the target
# can be provided.
        if len(self.letters)==0:
            if type(one)==type(''):
                print 'Warning, returning string.'
            return one
        else:
            answer = arg[abs(self.letters[0])-1]**(sign(self.letters[0]))
            for letter in self.letters[1:]:
                answer = answer * arg[abs(letter)-1]**(sign(letter))
            return answer
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
        strout = ''
        for letter in self.letters:
            strout = strout+key[letter+key_offset]
        return strout

def randomwalk(numgens,length):
    # This is the old randomword, which is really a random walk
    letterlist = range(1,numgens+1)+range(-1,-(numgens+1),-1)
    letters = []
    for n in range(length):
        letters.append(random.choice(letterlist))
    return word(letters)

def randomword(numgens,length):
    # random word not random walk
    letterlist = range(1,numgens+1)+range(-1,-(numgens+1),-1)
    letters = [random.choice(letterlist)] #first letter can be anything
    for n in range(length-1):
        lastletter=letters[len(letters)-1]
        letterlist.remove(-lastletter)  # Don't reverse the last letter
        letters.append(random.choice(letterlist))
        letterlist.append(-lastletter)  # restore the list.
    return word(letters)


def ishomologicallytrivial(w):
    reorder = copy.copy(w.letters)
    reorder.sort()
    reorder = freereduce(reorder)
    if len(reorder)>0:
        return False
    else:
        return True

def randomcommutator(numgens,length,stopafter=10000):
    rndwd = randomword(numgens,length)
    for n in range(stopafter):
        if ishomologicallytrivial(rndwd):
            print('Found commutator on '+str(n)+'th try.')
            return rndwd
        else:
            rndwd = randomword(numgens,length)
    print('Failed to find commutator.')
    return []

def guessmean(numgens,length):
# Here's some experimental verification that the expected progress of
# a random walk of length n in a free group of rank k is about
# (1-1/k)n.  (Not quite right for small values of n.)
    result = []
    for n in range(100):
        result.append(randomword(numgens,length).letters)
    result = map(len,result)
    return sum(result)/100

