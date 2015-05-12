import sys
from collections import defaultdict

def canonicalize(key):
    s1,s2 = key.split('|')
    s1 = ','.join(sorted(s1.split(',')))
    s2 = ','.join(sorted(s2.split(',')))
    key = '|'.join(sorted([s1, s2]))        
    return key

def alternatetopologies(key):
    key = canonicalize(key)
    s1,s2 = key.split('|')
    s1a, s1b = s1.split(',')
    s2a, s2b = s2.split(',')
    sa = ','.join(sorted([s1a, s2a]))
    sb = ','.join(sorted([s1b, s2b]))
    k1 = '|'.join(sorted([sa, sb]))
    
    sa = ','.join(sorted([s1a, s2b]))
    sb = ','.join(sorted([s1b, s2a]))
    k2 = '|'.join(sorted([sa, sb]))
    
    return canonicalize(k1), canonicalize(k2)

def addmissingtopos(quartets, val=0.):
    toadd = set()
    for q in quartets:
        for q1 in alternatetopologies(q):
            if q1 not in quartets:
                toadd.add(q1)
    for q1 in toadd:
        quartets[q1] = 0.


class Quartet:
    def __init__(self, qt):
#        ((self.a1, self.a2), (self.b1, self.b2)) = self.canonicalize(((a1, a2), (b1, b2)))
        if type(qt) == str:
            a, b = qt.split('),(') 
            a1, a2 = a.replace('(', '').replace(')', '').replace(';', '').split(',')
            b1, b2 = b.replace('(', '').replace(')', '').replace(';', '').split(',')
            self.data = self.canonicalize(((a1, a2), (b1, b2)))
        else:
            self.data = self.canonicalize(qt)
    def canonicalize(self, ((a1, a2), (b1, b2))):
        if a1 < a2:
            a = (a1, a2)
        else:
            a = (a2, a1)
        if b1 < b2:
            b = (b1, b2)
        else:
            b = (b2, b1)
        if a < b:
            return (a, b)
        else:
            return (b, a)
                
    def __eq__(self, other):
        return self.data == other.data
    def __hash__(self):
        return hash(self.data)
    def __repr__(self):
        return str(self.data)


def complementquartets(inputq, outputq):
    taxa = set()
    quartets = defaultdict(lambda:0)
    weights = []
    for i in inputq.readlines():
        quartet, weight = i.split(':')
        q1, q2 = quartet.split('|')
        a, b = sorted(q1.split(','))
        c, d = sorted(q2.split(','))
        taxa.add(a)
        taxa.add(b)
        taxa.add(c)
        taxa.add(d)
#        weights.append(float(weight))
        quartets[Quartet(((a, b), (c, d)))] = float(weight)

    parenquartets = []

    for i in taxa:
        for j in taxa:
            for k in taxa:
                for l in taxa:
                    if i < j and j < k and k < l:
                        parenquartets.append('((' + i + ',' + j + '),(' + k + ',' + l + '); ')
                        weights.append(quartets[Quartet(((i,j),(k,l)))])
                        parenquartets.append('((' + i + ',' + k + '),(' + j + ',' + l + '); ')
                        weights.append(quartets[Quartet(((i,k),(j,l)))])
                        parenquartets.append('((' + i + ',' + l + '),(' + k + ',' + j + '); ')
                        weights.append(quartets[Quartet(((i,l),(k,j)))])

    mw = max(weights)
    qd = []
        
    for quartet, weight in zip(parenquartets, weights):
        outputq.write(quartet + str(- weight) + '\n')
#        outputq.write(quartet + str(weight) + '\n')
    return list(taxa)
#    outputt.write('\n'.join(taxa))

#inputquartets = open(sys.argv[1])
#outputquartets = open(sys.argv[2],'w')

#complementquartets(inputquartets, outputquartets)
