import dendropy
import itertools
import collections
import pprint
from dendropy import treecalc
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

#inv52 is for the caterpillar tree 5 leaves (((a,b),c),d),e)
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
    viztrees = [dist[0][k].as_string("newick") for k in range(15)]
    print viztrees
    garbage = []
    for i in range(len(treelist)):
        smalltree = dendropy.Tree(treelist[i])
        smalltree.retain_taxa_with_labels(l)
        smalltree.update_splits()
        trace = [dendropy.treecalc.symmetric_difference(smalltree,dist[0][j]) for j in range(15)]
        if trace.count(0) > 1:
            garbage.append(['error: too many zeroes', smalltree.as_string("newick"), smalltree,  treelist[i], trace])
        elif 0 in trace:
            dist[1][trace.index(0)] = dist[1][trace.index(0)] + 1
        else:
            garbage.append([smalltree.as_newick_string(), treelist[i].as_newick_string(), min(trace)])
    return [dist, garbage]

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
    scorefuncs= [inv51, inv52, inv53]
    if rooted_shape in shapes:
        score_function = scorefuncs[shapes.index(rooted_shape)]
        score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    else:
        print 'error: quintet topology not in list'
    return score

#score_quintet takes a rooted quintet tree Q as input-is ladderized output of get_rooted_quintet
def score_quintet(Q,l, treelist):
    rooted_shape_str = Q.as_newick_string()
    rooted_shape = ''
    y = ['(', ')', ',']
    for i in range(len(rooted_shape_str)):
        if rooted_shape_str[i] in y:
            rooted_shape = rooted_shape + rooted_shape_str[i]
    shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    U = get_dist(l,treelist)
    [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15] = U[0][1]
    scorefuncs= [inv51, inv52, inv53]
    if rooted_shape in shapes:
        score_function = scorefuncs[shapes.index(rooted_shape)]
        score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    else:
        print 'error: quintet topology not in list'
    return score

#S is unrooted Astral tree, l is list of five taxa labels from S, i is the index of edge we want to root at in set 2n-3
def get_rooted_quintet(S,l,i):
    T = dendropy.Tree(S)
    edgelist = [e for e in T.postorder_edge_iter()]
    root_edge = edgelist[i]
    T.reroot_at_edge(root_edge)
    T.retain_taxa_with_labels(l)
    T.ladderize(ascending=False)
    return T

#need to fix score_quintet function 'pipleline'
def total_quintet_score(S,i,treelist):
    T = dendropy.Tree(S)  
    L = [n.get_node_str() for n in T.leaf_iter()]
    Q1 = list(itertools.combinations(L,5))  
    Q = [list(Q1[k]) for k in range(len(Q1))]
    Qtreeslist =[get_rooted_quintet(T,l,i) for l in Q]
    scorelist = []
    for j in range(len(Qtreeslist)):
        scorelist.append(score_quintet(Qtreeslist[j], Q[j],treelist))
    totalscore = sum(scorelist)
    return totalscore
    

def find_best_edge_by_total_quintet_score(S,treelist):
    T = dendropy.Tree(S) 
    Scores = [total_quintet_score(S,j,treelist) for j in range(len((S.leaf_nodes())) - 3)] 
    best_edge = Scores.index(min(Scores))
    edgelist = [e for e in T.postorder_edge_iter()]
    T.reroot_at_edge(edgelist[best_edge])
    return T


    
    
