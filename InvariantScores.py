import dendropy
import itertools
import collections
import pprint
import string
from dendropy import treecalc
import matrixmaker
import numpy as np
import copy

#from dendropy import *
#nuclear option to not type dendropy. I think this is frowned upon
#pp = pprint.PrettyPrinter(indent=4)
#inv3 gets the score for 3-taxon trees or 4-taxon trees.  Use for u1 the
#probability that the gene tree matches the species tree, then u2, u3 those that don't
def inv3(u1, u2, u3):
    a12 = min(u1-u2,0)
    a13 = min(u1-u3,0)
    score = abs(u2-u3) + (-1)*a12 + (-1)*a13
    return score

# a version of the penalty without equality....
def inv3I(u1, u2, u3):
    a12 = min(u1-u2,0)
    a13 = min(u1-u3,0)
    score (-1)*a12 + (-1)*a13
    return score

#inv51 is for the balanced species tree 5 leaves (((a,b),c),(d,e))
def inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15):
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0)
    a57 = min(u5-u7,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a25 + (-1)*a45 + (-1)*a57
    score2 = abs(u14 - u15) + abs(u11-u15) + abs(u10-u15) + abs(u9-u12) + abs(u8-u15) + abs(u7-u15) + abs(u6-u12) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3)
    score = score1 + score2
    return score

#inv52 is for the caterpillar tree 5 leaves ((((a,b),c),d),e)
def inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15): 
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0)
    a57 = min(u5-u7,0) 
    a32 = min(u3-u2,0)
    a36 = min(u3-u6,0)
    a65 = min(u6-u5,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a25 + (-1)*a45 + (-1)*a57 + (-1)*a32 + (-1)*a36 + (-1)*a65
    score2 = abs(u14-u15) + abs(u11-u15) + abs(u10-u15) + abs(u8-u15) + abs(u7-u15) + abs(u6-u9) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3 + u9-u12)
    score = score1 + score2
    return score

#inv53 is for the pseudocaterpillar tree (((a,b),(d,e)),c)
def inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15):
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a18 = min(u1-u8,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0) 
    a85 = min(u8-u5,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a18 + (-1)*a25 + (-1)*a45 + (-1)*a85
    score2 = abs(u14-u15) + abs(u12-u15) + abs(u10-u15) + abs(u9-u15) + abs(u8-u11) + abs(u7-u15) + abs(u6-u15) + abs(u5-u15) + abs(u4-u13) + abs(u2-u3)
    score = score1 + score2
    return score

#l=list of taxon labels as strings from listofgenetreeschoices,listofgenetrees=dendropy TreeList
def dist_counter(l,listofgenetrees):
    M = listofgenetrees.taxon_set
    t1 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[1]+'),' +l[2]+ ',(' +l[3]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M) 
    t2 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[1]+'),' +l[3]+ ',(' +l[2]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t3 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[1]+'),' +l[4]+ ',(' +l[2]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t4 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[2]+'),' +l[1]+ ',(' +l[3]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t5 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[2]+'),' +l[3]+ ',(' +l[1]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t6 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[2]+'),' +l[4]+ ',(' +l[1]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t7 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[3]+'),' +l[1]+ ',(' +l[2]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t8 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[3]+'),' +l[2]+ ',(' +l[1]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t9 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[3]+'),' +l[4]+ ',(' +l[1]+ ',' +l[2]+ '));\n', schema = 'newick', taxon_set = M)
    t10 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[4]+'),' +l[1]+ ',(' +l[2]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t11 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[4]+'),' +l[2]+ ',(' +l[1]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t12 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[4]+'),' +l[3]+ ',(' +l[1]+ ',' +l[2]+ '));\n', schema = 'newick', taxon_set = M)
    t13 = dendropy.Tree.get_from_string('[&U] ((' +l[1]+ ',' +l[2]+'),' +l[0]+ ',(' +l[3]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t14 = dendropy.Tree.get_from_string('[&U] ((' +l[1]+ ',' +l[3]+'),' +l[0]+ ',(' +l[2]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t15 = dendropy.Tree.get_from_string('[&U] ((' +l[1]+ ',' +l[4]+'),' +l[0]+ ',(' +l[2]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    counter1 = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15]
    counter2 = [0 for i in range(15)]
    return [counter1, counter2]

# I should merge this next with dist_counter after I test them both
def get_dist(l,treelist):
    #treelist = dendropy.TreeList(inputtreelist)
    dist = dist_counter(l, treelist)
    #pp.pprint(dist)
    #viztrees = [dist[0][k].as_string("newick") for k in range(15)]
    #print viztrees
    garbage = []
    incompletes = []
    labels = [n.label for n in treelist.taxon_set]
    for i in range(len(treelist)):
        smalltree = copy.deepcopy(treelist[i])
        if len(list(smalltree.leaf_nodes())) == len(labels):
            smalltree.retain_taxa_with_labels(l)
            smalltree.update_splits()
            trace = [dendropy.treecalc.symmetric_difference(smalltree,dist[0][j]) for j in range(15)]
            if trace.count(0) > 1:
                garbage.append(['error: too many zeroes', smalltree.as_string("newick"), smalltree,  treelist[i], trace])
            elif 0 in trace:
                dist[1][trace.index(0)] = dist[1][trace.index(0)] + 1
            else:
                garbage.append([smalltree.as_newick_string(), treelist[i].as_newick_string(), min(trace)])
        else:
            incompletes.append(str(i)+'th gene tree incomplete')
    return [dist, garbage, incompletes]

#S is a rooted species tree, l is a list of five taxa labels like in get_dist, etc.Picking the edge and the quintet are unresolved
def basic_score_quintet(S,l, treelist):
    T = dendropy.Tree(S)
    T.retain_taxa_with_labels(l)
    T.ladderize(ascending=False)
    rooted_shape_str = T.as_newick_string()
    rooted_shape = ''
    y = ['(', ')', ',']
    for i in range(len(rooted_shape_str)):
        if rooted_shape_str[i] in y:
            rooted_shape = rooted_shape + rooted_shape_str[i]
    shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    U = get_dist(l,treelist)
    [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15] = U[0][1] 
    #print U[0][1]
    scorefuncs= [inv51, inv52, inv53]
    if rooted_shape == shapes[0]:
        score = inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    elif rooted_shape == shapes[1]:
        score = inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    elif rooted_shape == shapes[2]:
        score = inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        #score_function = scorefuncs[shapes.index(rooted_shape)]
        #score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    else:
        print 'error: quintet topology ' + T.as_newick_string() + '  not in list'
        score = 0
    return score

#score_quintet takes a rooted quintet tree Q as input-is ladderized output of get_rooted_quintet
#def score_quintet(Q,l, treelist):
    #rooted_shape_str = Q.as_newick_string()
    #rooted_shape = ''
    #y = ['(', ')', ',']
    #for i in range(len(rooted_shape_str)):
        #if rooted_shape_str[i] in y:
            #rooted_shape = rooted_shape + rooted_shape_str[i]
    #shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    #U = get_dist(l,treelist)
    #[u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15] = U[0][1]
    #scorefuncs= [inv51, inv52, inv53]
    #if rooted_shape in shapes:
        #score_function = scorefuncs[shapes.index(rooted_shape)]
        #score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    #else:
        #print 'error: quintet topology ' + Q.as_newick_string() + '  not in list'
        #score = 0
    #return score

#S is unrooted Astral tree, l is list of five taxa labels from S, i is the index of edge we want to root at in set 2n-3
def get_rooted_quintet(S,l,i):
    T = dendropy.Tree(S)
    edgelist = [e for e in T.postorder_edge_iter()]
    root_edge = edgelist[i]
    T.reroot_at_edge(root_edge)
    T.retain_taxa_with_labels(l)
    T.ladderize(ascending=False)
    return T

#This function scores an edge based on all five-element subsets of the leaf set of the species tree.  This will be replaced or 
#modified to score an edge based on 2n-3 five-elements subsets where n is the numbers of species.  i indexes in the edges in 
#find _best_edge_by_total_quintet_score(S,treelist):
def total_quintet_score(S,i,treelist):
    #T needs to be S rooted at edge i!!!!
    T = dendropy.Tree(S) 
    edgelist = [e for e in T.postorder_edge_iter()]
    #root_edge = edgelist[i]
    T.reroot_at_edge(edgelist[i])
    L = [n.taxon.label for n in T.leaf_nodes()]
    #print L
    Q1 = list(itertools.combinations(L,5))
    #print Q1 
    Q = [list(Q1[k]) for k in range(len(Q1))]
    #print Q
    #Qtreeslist =[get_rooted_quintet(T,l,i) for l in Q]
    scorelist = []
    print 'OOOGA BOOGA'
    #for j in range(len(Qtreeslist)):
    for j in range(len(Q)):
        H = basic_score_quintet(T, Q[j],treelist)
        print H
        scorelist.append(basic_score_quintet(T, Q[j],treelist))
    totalscore = sum(scorelist)
    return totalscore
    

def find_best_edge_by_total_quintet_score(S,treelist):
    T = dendropy.Tree(S) 
    edgelist = [e for e in T.postorder_edge_iter()]
    numedges = len(edgelist)
    Scores = [total_quintet_score(S,j,treelist) for j in range(numedges)] 
    best_edge = Scores.index(min(Scores))
    T.reroot_at_edge(edgelist[best_edge])
    return T


class Quartet:
    def __init__(self, qt):
#        ((self.a1, self.a2), (self.b1, self.b2)) = self.canonicalize(((a1, a2), (b1, b2)))
        if type(qt) == str:
            a, b = qt.split('),(') 
            a1, a2 = a.replace('(', '').replace(')', '').replace(';', '').split(',')
            b1, b2 = b.replace('(', '').replace(')', '').replace(';', '').split(',')
            self.data = self.canonicalize(((int(a1), int(a2)), (int(b1), int(b2))))
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

#quartetsfile is a string name of file, file must be in working directory
class QuartetsInfo:
    def __init__(self,quartetsfile):
        self.filename = quartetsfile
        self.quartet_dict_ = {}
        self.quartet_labels_dict_ = {}
        #with open(quartetsfile) as f:
            #for line in f:
                #(key,val) = line.split()
                #self.quartet_dict[key] = int(val)
    def quartet_dict(self):
        if len(self.quartet_dict_):
            return self.quartet_dict_
        d = {}
        with open(self.filename) as f:
            for line in f:
                (key,val) = line.split()
                d[Quartet(key)] = float(val)
        self.quartet_dict_ = d
        return d

#L is a list of four labels that are strings-alphabetical order required
    def get_freqs(self,L):
        s = self.quartet_dict()
        #q = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));']
        #u = [0,0,0]
        #for i in range(3):
        #    if q[i] in s:
        #        u[i] = s[q[i]]
        
        q0 = Quartet(((L[0], L[1]), (L[2], L[3])))
        q1 = Quartet(((L[0], L[2]), (L[1], L[3])))
        q2 = Quartet(((L[0], L[3]), (L[1], L[2])))

        u = [0, 0, 0]
        if q0 in s:
            u[0] = s[q0]
        else:
            print "ERROR", q0
        if q1 in s:
            u[1] = s[q1]
        else:
            print "ERROR", q1
        if q2 in s:
            u[2] = s[q2]
        else:
            print "ERROR", q2

        return u

        # u = [0,0,0]
        # q0 = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[1]+'),('+L[3]+','+L[2]+'));','(('+L[1]+','+L[0]+'),('+L[2]+','+L[3]+'));','(('+L[1]+','+L[0]+'),('+L[3]+','+L[2]+'));', '(('+L[2]+','+L[3]+'),('+L[0]+','+L[1]+'));','(('+L[2]+','+L[3]+'),('+L[1]+','+L[0]+'));','(('+L[3]+','+L[2]+'),('+L[0]+','+L[1]+'));','(('+L[3]+','+L[2]+'),('+L[1]+','+L[0]+'));']
        # q1 = ['(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[3]+','+L[1]+'));','(('+L[2]+','+L[0]+'),('+L[1]+','+L[3]+'));','(('+L[2]+','+L[0]+'),('+L[3]+','+L[1]+'));','(('+L[1]+','+L[3]+'),('+L[0]+','+L[2]+'));','(('+L[1]+','+L[3]+'),('+L[2]+','+L[0]+'));','(('+L[3]+','+L[1]+'),('+L[0]+','+L[2]+'));','(('+L[3]+','+L[1]+'),('+L[2]+','+L[0]+'));']
        # q2 = ['(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));', '(('+L[0]+','+L[3]+'),('+L[2]+','+L[1]+'));','(('+L[3]+','+L[0]+'),('+L[1]+','+L[2]+'));', '(('+L[3]+','+L[0]+'),('+L[2]+','+L[1]+'));','(('+L[1]+','+L[2]+'),('+L[0]+','+L[3]+'));','(('+L[1]+','+L[2]+'),('+L[3]+','+L[0]+'));','(('+L[2]+','+L[1]+'),('+L[0]+','+L[3]+'));','(('+L[2]+','+L[1]+'),('+L[3]+','+L[0]+'));']
        # for i in range(8):
        #     if q0[i] in s:
        #         u[0] = s[q0[i]]
        # for i in range(8):
        #     if q1[i] in s:
        #         u[1] = s[q1[i]]
        # for i in range(8):
        #     if q2[i] in s:
        #         u[2] = s[q2[i]]
        # return u
        
#L same as get_freqs, i,j are indices of L specified as cherry of interest, func is string with name of penalty function 
    def quartet_score(self,L,i,j,func):
        a = [i,j]
        a.sort()
        f1, f2, f3 = self.get_freqs(L)
        if a in [[0,1], [2,3]]:
            u1, u2, u3 = f1, f2, f3
        elif a in [[0,2], [1,3]]:
            u1, u2, u3 = f2, f1, f3
        elif a in [[0,3], [1,2]]:
            u1, u2, u3 = f3, f1, f2
        else:
            print "problem with quartet score input indices"
            print L, i, j

        #return u1
#        if func == 'inv3':
#            score = inv3(u1, u2, u3)
        return u1
        a12 = min(u1-u2,0)
        a13 = min(u1-u3,0)
        score = abs(u2-u3) + (-1)*a12 + (-1)*a13
        return score

    def quartet_labels_dict(self):
        if len(self.quartet_labels_dict_):
            return self.quartet_labels_dict_
        g = self.quartet_dict()
        h = {}
        gkeys = g.keys()
        # for i in range(len(gkeys)):

        #     keycopy = gkeys[i]
        #     labellist = keycopy.split(',')
        #     for j in range(4):
        #         labellist[j] = labellist[j].translate(None,string.punctuation)        
        #     h[gkeys[i]] = labellist

        for i in gkeys:
            h[i] = (i.data[0][0], i.data[0][1], i.data[1][0], i.data[1][1])

        self.quartet_labels_dict_ = h
        return h

#D1,D2  are two taxa labels as strings in alphabetical order
    def score_double(self,D1,D2):
        g = self.quartet_dict()
        h = self.quartet_labels_dict()
        #print g
        #print h
        score = 0
        for k in g:
            if '(' + D1 + ',' + D2 + ')' in k:
                #print k
                newscore = self.quartet_score(h[k], h[k].index(D1), h[k].index(D2),'inv3')
                #print newscore
                score = score + newscore
            if '(' + D2 + ',' + D1 + ')' in k:
                #print k
                newscore = self.quartet_score(h[k], h[k].index(D1), h[k].index(D2),'inv3')
                #print newscore
                score = score + newscore
        return score

    #def penalty_score(self, subset1, subset2):
        
    def score_tree(self,treequartetsfile):
        treelist = []
        with open(treequartetsfile) as f:
            for line in f:
                (keygen,val) = line.split()
                labellist = keygen.split(',')
                for j in range(4):
                    labellist[j] = labellist[j].translate(None,string.punctuation)        
                treelist.append(labellist)
        #print treelist
        score = 0
        for k in range(len(treelist)):
                newscore = self.quartet_score(treelist[k],0,1,'inv3')
                #print newscore
                score = score + newscore
        return score

class QuartetsInfoInequalitiesOnly:
    def __init__(self,quartetsfile):
        self.filename = quartetsfile
        #self.quartet_dict = {}
        #with open(quartetsfile) as f:
            #for line in f:
                #(key,val) = line.split()
                #self.quartet_dict[key] = int(val)
    def quartet_dict(self):
        d = {}
        with open(self.filename) as f:
            for line in f:
                (key,val) = line.split()
                d[key] = int(val)
        return d

#L is a list of four labels that are strings
    def get_freqs(self,L):
        s = self.quartet_dict()
        #q = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));']
        #u = [0,0,0]
        #for i in range(3):
        #    if q[i] in s:
        #        u[i] = s[q[i]]
        u = [0,0,0]
        q0 = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[1]+'),('+L[3]+','+L[2]+'));','(('+L[1]+','+L[0]+'),('+L[2]+','+L[3]+'));','(('+L[1]+','+L[0]+'),('+L[3]+','+L[2]+'));', '(('+L[2]+','+L[3]+'),('+L[0]+','+L[1]+'));','(('+L[2]+','+L[3]+'),('+L[1]+','+L[0]+'));','(('+L[3]+','+L[2]+'),('+L[0]+','+L[1]+'));','(('+L[3]+','+L[2]+'),('+L[1]+','+L[0]+'));']
        q1 = ['(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[3]+','+L[1]+'));','(('+L[2]+','+L[0]+'),('+L[1]+','+L[3]+'));','(('+L[2]+','+L[0]+'),('+L[3]+','+L[1]+'));','(('+L[1]+','+L[3]+'),('+L[0]+','+L[2]+'));','(('+L[1]+','+L[3]+'),('+L[2]+','+L[0]+'));','(('+L[3]+','+L[1]+'),('+L[0]+','+L[2]+'));','(('+L[3]+','+L[1]+'),('+L[2]+','+L[0]+'));']
        q2 = ['(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));', '(('+L[0]+','+L[3]+'),('+L[2]+','+L[1]+'));','(('+L[3]+','+L[0]+'),('+L[1]+','+L[2]+'));', '(('+L[3]+','+L[0]+'),('+L[2]+','+L[1]+'));','(('+L[1]+','+L[2]+'),('+L[0]+','+L[3]+'));','(('+L[1]+','+L[2]+'),('+L[3]+','+L[0]+'));','(('+L[2]+','+L[1]+'),('+L[0]+','+L[3]+'));','(('+L[2]+','+L[1]+'),('+L[3]+','+L[0]+'));']
        for i in range(8):
            if q0[i] in s:
                u[0] = s[q0[i]]
        for i in range(8):
            if q1[i] in s:
                u[1] = s[q1[i]]
        for i in range(8):
            if q2[i] in s:
                u[2] = s[q2[i]]
        return u
        
        
#L same as get_freqs, i,j are indices of L specified as cherry of interest
    def quartet_score(self,L,i,j):
        a = [i,j]
        a.sort()
        f1, f2, f3 = self.get_freqs(L)
        if a in [[0,1], [2,3]]:
            u1, u2, u3 = f1, f2, f3
        elif a in [[0,2], [1,3]]:
            u1, u2, u3 = f2, f1, f3
        elif a in [[0,3], [1,2]]:
            u1, u2, u3 = f3, f1, f2
        else:
            print "problem with quartet score input indices"
            print L, i, j
        a12 = min(u1-u2,0)
        a13 = min(u1-u3,0)
        score = (-1)*a12 + (-1)*a13
        return score

    def quartet_labels_dict(self):
        g = self.quartet_dict()
        h = {}
        gkeys = g.keys()
        for i in range(len(gkeys)):
            keycopy = gkeys[i]
            labellist = keycopy.split(',')
            for j in range(4):
                labellist[j] = labellist[j].translate(None,string.punctuation)        
            h[gkeys[i]] = labellist
        return h

#D1,D2  are two taxa labels as strings in alphabetical order
    def score_double(self,D1,D2):
        g = self.quartet_dict()
        h = self.quartet_labels_dict()
        #print g
        #print h
        score = 0
        for k in g:
            if '(' + D1 + ',' + D2 + ')' in k:
                #print k
                newscore = self.quartet_score(h[k], h[k].index(D1), h[k].index(D2))
                #print newscore
                score = score + newscore
        return score

    #def penalty_score(self, subset1, subset2):
        
    def score_tree(self,treequartetsfile):
        treelist = []
        with open(treequartetsfile) as f:
            for line in f:
                (keygen,val) = line.split()
                labellist = keygen.split(',')
                for j in range(4):
                    labellist[j] = labellist[j].translate(None,string.punctuation)        
                treelist.append(labellist)
        #print treelist
        score = 0
        for k in range(len(treelist)):
                newscore = self.quartet_score(treelist[k],0,1)
                #print newscore
                score = score + newscore
        return score


class SubsetPenalties:
    def __init__(self,labels,setlist,quartetsfile):
        self.quartetsinfo = QuartetsInfo(quartetsfile)
        self.setlist = setlist
        self.matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        self.matrix = self.matrix_maker.matrix()
        self.labels = self.matrix_maker.labels
        self.SM_ = None
    #def pairs(self):
        #clist = copy.copy(self.setlist)
        #clabels = copy.copy(self.labels)
        #p = []
        #while len(clist) > 0:
            #one = clist.pop(0)
            #two = list(set(clabels) - set(one))
            #if two in clist:
                #p.append([one, two])
                #clist.remove(two)
            #else:
                #print "missing pair from input setlist"
        #return p
        
        
    def add_quartets(self,set1,set2):
        stuntlabels = []
        for i in range(len(self.labels)):
            stuntlabels.append(self.labels[i])
        compl = set1 + set2
        #print stuntlabels, compl
        for s in compl:
            stuntlabels.remove(s)
            #print stuntlabels
        #for t in set2:
            #stuntlabels.remove(t)
            #print stuntlabels
        A = itertools.combinations(stuntlabels,2)
        AA = list(A)
        SminusApairs = [list(AA[j]) for j in range(len(AA))]
        A1A2 = []
        for k in range(len(set1)):
            for l in range(len(set2)):
                A1A2.append([set1[k],set2[l]])
        AQ = []
        for m in range(len(A1A2)):
            for n in range(len(SminusApairs)):
                tmp = A1A2[m] + SminusApairs[n]
                tmp.sort()
                AQ.append([tmp, tmp.index(A1A2[m][0]), tmp.index(A1A2[m][1])])
        #print AQ
        return AQ

    def subtract_quartets(self,set1,set2):
        a1 = itertools.combinations(set1,2)
        a2 = itertools.combinations(set2,2)
        aa1 = list(a1)
        aa2 = list(a2)
        A1 = [list(aa1[i]) for i in range(len(aa1))]
        A2 = [list(aa2[j]) for j in range(len(aa2))]
        SQ = []
        for m in range(len(A1)):
            for n in range(len(A2)):
                tmp = A1[m] + A2[n]
                tmp.sort()
                SQ.append([tmp,tmp.index(A1[m][0]), tmp.index(A1[m][1])])
        #print SQ
        return SQ

#need to enter score1, score2 for set1, set2 now: later can make this flexible? sets of size 1 and 2 are 0 and fixed. 
    def penalty_score(self,set1,set2,score1,score2):
#        g = self.quartetsinfo.quartet_dict()
 #       h = self.quartetsinfo.quartet_labels_dict()
        #labs = self.labels
        #print labs, set1, set2
        AQ = self.add_quartets(set1,set2)
        #print AQ
        SQ = self.subtract_quartets(set1,set2)
        #print SQ
        score = score1 + score2
        for i in range(len(AQ)):
            score = score+self.quartetsinfo.quartet_score(AQ[i][0], AQ[i][1], AQ[i][2], None)
            #print score
        for j in range(len(SQ)): 
            score = score-self.quartetsinfo.quartet_score(SQ[j][0], SQ[j][1], SQ[j][2], None)
        return score

#problem is finding min NONZERO entry here. 
    #def scored_matrix(self):
        #SM = np.copy(self.matrix)
        #for j in range(len(self.setlist)):
            #SM[j,j] = np.inf
        #for i in range(len(self.setlist)):
            #if len(self.setlist[i]) == 1:
                #SM[i,i] = 0
            #if len(self.setlist[i]) == 2:
                #m = self.setlist.index([self.setlist[i][0]])
                #n = self.setlist.index([self.setlist[i][1]])
                #SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                #SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            #if len(self.setlist[i]) > 2:
                #clabels = copy.copy(self.setlist[i])
                #clist = copy.copy(self.setlist)
                #for s in clist:
                    #if SM[i,self.setlist.index(s)] == np.inf:
                        #clist.remove(s)
                #while len(clist) > 0:
                        #one = clist.pop(0)
                        #two = list(set(clabels) - set(one))
                        #one.sort()
                        #two.sort()
                        #if two in clist:
                            #clist.remove(two)
                            #onedex = self.setlist.index(one)
                            #twodex = self.setlist.index(two)
                            #SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                            #SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
                            #SM[i,onedex] = np.inf
                            #SM[i,twodex] = np.inf
        #print SM
        #return SM

    def scored_matrix(self):
        if self.SM_ is not None:
            return self.SM_
        self.SM_ = np.zeros((len(self.setlist),len(self.setlist)))
        SM = self.SM_
        SM.fill(np.inf)
        for i in range(len(self.setlist)):
            if len(self.setlist[i]) == 1:
                SM[i,i] = 0
            if len(self.setlist[i]) == 2:
                m = self.setlist.index([self.setlist[i][0]])
                n = self.setlist.index([self.setlist[i][1]])
                SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            if len(self.setlist[i]) > 2:
                clabels = copy.copy(self.setlist[i])
                clist = copy.copy(self.setlist)
                while len(clist) > 0:
                        one = clist.pop(0)
                        onedex = self.setlist.index(one)
                        if self.matrix[i,onedex] == 1:
                            two = list(set(clabels) - set(one))
                            one.sort()
                            two.sort()
                            if two in clist:
                                clist.remove(two)
                                twodex = self.setlist.index(two)
                                SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                                SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
        return SM

    #def clades(self):
        #SM = self.scored_matrix()
        #L = []
        #temprows =[self.setlist.index(self.setlist[-1]), self.setlist.index(self.setlist[-1])]
        #while temprow > 0:
                #onedex = SM[temprow].argmin()
                #one = self.setlist[onedex]
                #two = list(set(self.setlist[temprow])-set(one))
                #two.sort()
                #twodex = self.setlist.index(two)



    def clades(self):
        SM = self.scored_matrix()
        #L = []
        L = [self.setlist[-1]]
        tempsets = [self.setlist[-1]]
        while len(tempsets) > 0:
            tmp = tempsets.pop(0)
            if len(tmp) == 1:
                L.append(tmp)
            else:
                rowdex = self.setlist.index(tmp)
                onedex = SM[rowdex].argmin()
                one = self.setlist[onedex]
                two = list(set(self.setlist[rowdex]) - set(one))
                two.sort()
                L.append(one)
                L.append(two)
                if len(one) > 1:
                    tempsets.append(one)
                if len(two) > 1:
                    tempsets.append(two)
        return L

    #def tree(self):
        #L = self.clades()
        #start = L.pop(0)
        #s = ''
        #for i in range(len(start)-1):
            #s = s + start[i]+ ','
        #s = s + start[-1]
        #s = '(' + s + ')'
        #while len(L) > 0:
            #tmp = L.pop(0)
            #if len(tmp)> 1:
                #spot = s.index(tmp[0])
                #t = '('
                #for j in range(len(tmp)-1):
                    #t = t + tmp[j] +','
                #t = t + tmp[-1] + ')'
                #for k in range(len(tmp)):
                    #s.replace(tmp[k] + ',','')
                    #s.replace(tmp[k],'')
                #s = s[:spot] + t + s[spot:]
        #return s

    def tree(self):
        L = self.clades()
        start = L.pop(0)
        s = ''
        for i in range(len(start)):
            s = s + start[i]+ ','
        #s = s + start[-1]
        s = '(' + s + ')'
        while len(L) > 0:
            tmp = L.pop(0)
            if len(tmp) > 1:
                first = tmp[0]
                #find the first taxon, to use first letter as string index place
                spot = s.index(first[0])
                t = '('
                for j in range(len(tmp)):
                    t = t + tmp[j] + ','
                t = t + '),'
                for k in range(len(tmp)):
                    s = s.replace(tmp[k] + ',', '')
                    s = s.replace(tmp[k], '')
                s = s[:spot] + t + s[spot:] 
        s = s.replace(',)',')')
        s = s + ';'
        return s

class SubsetPenaltiesInequalitiesOnly:
    def __init__(self,labels,setlist,quartetsfile):
        self.quartetsinfo = QuartetsInfoInequalitiesOnly(quartetsfile)
        self.setlist = setlist
        self.matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        self.matrix = self.matrix_maker.matrix()
        self.labels = self.matrix_maker.labels

    #def pairs(self):
        #clist = copy.copy(self.setlist)
        #clabels = copy.copy(self.labels)
        #p = []
        #while len(clist) > 0:
            #one = clist.pop(0)
            #two = list(set(clabels) - set(one))
            #if two in clist:
                #p.append([one, two])
                #clist.remove(two)
            #else:
                #print "missing pair from input setlist"
        #return p
        
        
    def add_quartets(self,set1,set2):
        stuntlabels = []
        for i in range(len(self.labels)):
            stuntlabels.append(self.labels[i])
        compl = set1 + set2
        #print stuntlabels, compl
        for s in compl:
            stuntlabels.remove(s)
            #print stuntlabels
        #for t in set2:
            #stuntlabels.remove(t)
            #print stuntlabels
        A = itertools.combinations(stuntlabels,2)
        AA = list(A)
        SminusApairs = [list(AA[j]) for j in range(len(AA))]
        A1A2 = []
        for k in range(len(set1)):
            for l in range(len(set2)):
                A1A2.append([set1[k],set2[l]])
        AQ = []
        for m in range(len(A1A2)):
            for n in range(len(SminusApairs)):
                tmp = A1A2[m] + SminusApairs[n]
                tmp.sort()
                AQ.append([tmp, tmp.index(A1A2[m][0]), tmp.index(A1A2[m][1])])
        #print AQ
        return AQ

    def subtract_quartets(self,set1,set2):
        a1 = itertools.combinations(set1,2)
        a2 = itertools.combinations(set2,2)
        aa1 = list(a1)
        aa2 = list(a2)
        A1 = [list(aa1[i]) for i in range(len(aa1))]
        A2 = [list(aa2[j]) for j in range(len(aa2))]
        SQ = []
        for m in range(len(A1)):
            for n in range(len(A2)):
                tmp = A1[m] + A2[n]
                tmp.sort()
                SQ.append([tmp,tmp.index(A1[m][0]), tmp.index(A1[m][1])])
        #print SQ
        return SQ

#need to enter score1, score2 for set1, set2 now: later can make this flexible? sets of size 1 and 2 are 0 and fixed. 
    def penalty_score(self,set1,set2,score1,score2):
        g = self.quartetsinfo.quartet_dict()
        h = self.quartetsinfo.quartet_labels_dict()
        #labs = self.labels
        #print labs, set1, set2
        AQ = self.add_quartets(set1,set2)
        #print AQ
        SQ = self.subtract_quartets(set1,set2)
        #print SQ
        score = score1 + score2
        for i in range(len(AQ)):
            score = score+self.quartetsinfo.quartet_score(AQ[i][0], AQ[i][1], AQ[i][2])
            #print score
        for j in range(len(SQ)): 
            score = score-self.quartetsinfo.quartet_score(SQ[j][0], SQ[j][1], SQ[j][2])
        return score

#problem is finding min NONZERO entry here. 
    #def scored_matrix(self):
        #SM = np.copy(self.matrix)
        #for j in range(len(self.setlist)):
            #SM[j,j] = np.inf
        #for i in range(len(self.setlist)):
            #if len(self.setlist[i]) == 1:
                #SM[i,i] = 0
            #if len(self.setlist[i]) == 2:
                #m = self.setlist.index([self.setlist[i][0]])
                #n = self.setlist.index([self.setlist[i][1]])
                #SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                #SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            #if len(self.setlist[i]) > 2:
                #clabels = copy.copy(self.setlist[i])
                #clist = copy.copy(self.setlist)
                #for s in clist:
                    #if SM[i,self.setlist.index(s)] == np.inf:
                        #clist.remove(s)
                #while len(clist) > 0:
                        #one = clist.pop(0)
                        #two = list(set(clabels) - set(one))
                        #one.sort()
                        #two.sort()
                        #if two in clist:
                            #clist.remove(two)
                            #onedex = self.setlist.index(one)
                            #twodex = self.setlist.index(two)
                            #SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                            #SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
                            #SM[i,onedex] = np.inf
                            #SM[i,twodex] = np.inf
        #print SM
        #return SM

    def scored_matrix(self):
        SM = np.zeros((len(self.setlist),len(self.setlist)))
        SM.fill(np.inf)
        for i in range(len(self.setlist)):
            if len(self.setlist[i]) == 1:
                SM[i,i] = 0
            if len(self.setlist[i]) == 2:
                m = self.setlist.index([self.setlist[i][0]])
                n = self.setlist.index([self.setlist[i][1]])
                SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            if len(self.setlist[i]) > 2:
                clabels = copy.copy(self.setlist[i])
                clist = copy.copy(self.setlist)
                while len(clist) > 0:
                        one = clist.pop(0)
                        onedex = self.setlist.index(one)
                        if self.matrix[i,onedex] == 1:
                            two = list(set(clabels) - set(one))
                            one.sort()
                            two.sort()
                            if two in clist:
                                clist.remove(two)
                                twodex = self.setlist.index(two)
                                SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                                SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
        return SM

    #def clades(self):
        #SM = self.scored_matrix()
        #L = []
        #temprows =[self.setlist.index(self.setlist[-1]), self.setlist.index(self.setlist[-1])]
        #while temprow > 0:
                #onedex = SM[temprow].argmin()
                #one = self.setlist[onedex]
                #two = list(set(self.setlist[temprow])-set(one))
                #two.sort()
                #twodex = self.setlist.index(two)



    def clades(self):
        SM = self.scored_matrix()
        #L = []
        L = [self.setlist[-1]]
        tempsets = [self.setlist[-1]]
        while len(tempsets) > 0:
            tmp = tempsets.pop(0)
            if len(tmp) == 1:
                L.append(tmp)
            else:
                rowdex = self.setlist.index(tmp)
                onedex = SM[rowdex].argmin()
                one = self.setlist[onedex]
                two = list(set(self.setlist[rowdex]) - set(one))
                two.sort()
                L.append(one)
                L.append(two)
                if len(one) > 1:
                    tempsets.append(one)
                if len(two) > 1:
                    tempsets.append(two)
        return L

    #def tree(self):
        #L = self.clades()
        #start = L.pop(0)
        #s = ''
        #for i in range(len(start)-1):
            #s = s + start[i]+ ','
        #s = s + start[-1]
        #s = '(' + s + ')'
        #while len(L) > 0:
            #tmp = L.pop(0)
            #if len(tmp)> 1:
                #spot = s.index(tmp[0])
                #t = '('
                #for j in range(len(tmp)-1):
                    #t = t + tmp[j] +','
                #t = t + tmp[-1] + ')'
                #for k in range(len(tmp)):
                    #s.replace(tmp[k] + ',','')
                    #s.replace(tmp[k],'')
                #s = s[:spot] + t + s[spot:]
        #return s

    def tree(self):
        L = self.clades()
        start = L.pop(0)
        s = ''
        for i in range(len(start)):
            s = s + start[i]+ ','
        #s = s + start[-1]
        s = '(' + s + ')'
        while len(L) > 0:
            tmp = L.pop(0)
            if len(tmp) > 1:
                first = tmp[0]
                #find the first taxon, to use first letter as string index place
                spot = s.index(first[0])
                t = '('
                for j in range(len(tmp)):
                    t = t + tmp[j] + ','
                t = t + '),'
                for k in range(len(tmp)):
                    s = s.replace(tmp[k] + ',', '')
                    s = s.replace(tmp[k], '')
                s = s[:spot] + t + s[spot:] 
        s = s.replace(',)',')')
        s = s + ';'
        return s

#dendropy_clades_tree takes list of Taxa and list of clades, returns a tree
def dendropy_clades_tree(Taxa,cladelist):
    trees = ''
    for i in range(len(cladelist)):
        #print i
        compl = list(set(Taxa)-set(cladelist[i]))
        #print compl
        if len(compl) > 0:
            half1 = '(('
            for k in range(len(cladelist[i])-1):
                half1 = half1 + str(cladelist[i][k]) + ','
            half1 = half1 + str(cladelist[i][-1]) +'),'
            half2 = '('
            for j in range(len(compl)-1):
                half2 = half2 + str(compl[j]) + ','
            half2  = half2 + str(compl[-1]) + '));'
            tree = half1+half2
            trees = trees+tree
    tmplist = dendropy.TreeList.get_from_string(trees,'newick')
    finaltree = tmplist.consensus(min_freq=0)
    return finaltree

        


#powerset makes the list of all subsets of a label set
def powerset(labellist):
    P = []
    for j in range(1,len(labellist)+1):
        tmplist = list(itertools.combinations(labellist,j))
        tmplist = [list(tmplist[i]) for i in range(len(tmplist))]
        for k in range(len(tmplist)):
            tmplist[k].sort()
        P = P + tmplist
    return P

class QuintetsMaker:
    def __init__(self,labels):
        self.labels = labels

#here l is list of five labels
    def quintets(self,l):
        t1 = '((' +l[0]+ ',' +l[1]+'),' +l[2]+ ',(' +l[3]+ ',' +l[4]+ '));\n'
        t2 = '((' +l[0]+ ',' +l[1]+'),' +l[3]+ ',(' +l[2]+ ',' +l[4]+ '));\n'
        t3 = '((' +l[0]+ ',' +l[1]+'),' +l[4]+ ',(' +l[2]+ ',' +l[3]+ '));\n'
        t4 = '((' +l[0]+ ',' +l[2]+'),' +l[1]+ ',(' +l[3]+ ',' +l[4]+ '));\n'
        t5 = '((' +l[0]+ ',' +l[2]+'),' +l[3]+ ',(' +l[1]+ ',' +l[4]+ '));\n'
        t6 = '((' +l[0]+ ',' +l[2]+'),' +l[4]+ ',(' +l[1]+ ',' +l[3]+ '));\n'
        t7 = '((' +l[0]+ ',' +l[3]+'),' +l[1]+ ',(' +l[2]+ ',' +l[4]+ '));\n'
        t8 = '((' +l[0]+ ',' +l[3]+'),' +l[2]+ ',(' +l[1]+ ',' +l[4]+ '));\n'
        t9 = '((' +l[0]+ ',' +l[3]+'),' +l[4]+ ',(' +l[1]+ ',' +l[2]+ '));\n'
        t10 = '((' +l[0]+ ',' +l[4]+'),' +l[1]+ ',(' +l[2]+ ',' +l[3]+ '));\n'
        t11 = '((' +l[0]+ ',' +l[4]+'),' +l[2]+ ',(' +l[1]+ ',' +l[3]+ '));\n'
        t12 = '((' +l[0]+ ',' +l[4]+'),' +l[3]+ ',(' +l[1]+ ',' +l[2]+ '));\n'
        t13 = '((' +l[1]+ ',' +l[2]+'),' +l[0]+ ',(' +l[3]+ ',' +l[4]+ '));\n'
        t14 = '((' +l[1]+ ',' +l[3]+'),' +l[0]+ ',(' +l[2]+ ',' +l[4]+ '));\n'
        t15 = '((' +l[1]+ ',' +l[4]+'),' +l[0]+ ',(' +l[2]+ ',' +l[3]+ '));\n'
        q = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15]
        return q
    
    def fives(self):
        s = list(itertools.combinations(self.labels,5))
        return s
    
    def allquintets(self):
        f = self.fives()
        All = []
        for j in range(len(f)):
            All = All + self.quintets(f[j])
        return All
#give filename as a string
    def makefile(self,filename):
        f = open(filename,'w')
        for k in range(len(self.allquintets())):
            f.write(self.allquintets()[k])
        f.close()

class FileClades:
    def __init__(self,filename):
        self.filename = filename

    def clades(self):
        C = []
        f = open(self.filename, 'r')
        for line in f:
            if '{' in line:
                C.append(line)
        f.close()
        for i in range(len(C)):
            C[i] = C[i].replace('{','')
        for j in range(len(C)):
            C[j] = C[j].replace('}','')
        for k in range(len(C)):
            C[k] = C[k].replace('\n','')
        for n in range(len(C)):
            C[n] = C[n].replace(' ','')
        D = []
        for a in range(len(C)):
            if len(C[a]) > 0:
                c = [int (i) for i in C[a].split(',')]
                c.sort()
                D.append(c)
        D.reverse()
        Taxa = []
        for i in range(len(D)):
            if len(D[i]) == 1:
                Taxa = Taxa + D[i]
        Taxa.sort()
        #Taxa.sort(key=float)
        if (Taxa in D) == False:
            D.append(Taxa)
        return D



#labels here is a list of just five labels
#def quintets(labels):
    

#this build tree takes an output from a clades function from SubsetPenalties.clades
    

#This class was not working well, replace with score_matrix function in class SubsetPenalties 
#matrix input called m  here is an instance of matrixmaker.MatrixMaker(labels,setlist)
#class ScoredMatrix:
    #def __init__(self,labels,setlist,quartetsfile):
        #self.start_matrix = matrixmaker.MatrixMaker(labels,setlist)
        #self.quartetsfile = quartetsfile
        #self.labels = self.start_matrix.labels
        #self.setlist = self.start_matrix.setlist
        #self.quartetsinfo = QuartetsInfo(quartetsfile)
        #self.matrixy = self.start_matrix.matrix()
        #self.quartet_dict = {}

    #def scored_matrix(self):
        #K  = self.matrixy
        #L = self.setlist
        #for i in range(len(L)): 
            #if len(L[i]) == 2:
               #for j in range(i-1):
                    #if K[i, j] == 1:
                        #K[i,j] = self.quartetsinfo.score_double(L[i][0], L[i][1])
        #return K            
        
