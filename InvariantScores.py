import dendropy
#from dendropy import *
#nuclear option to not type dendropy. I think this is frowned upon

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
    t1 = dendropy.Tree.get_from_string('((' +l[0]+ ',' +l[1]+'),' +l[2]+ ',(' +l[3]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M) 
    t2 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[1]+'),' +l[3]+ ',(' +l[2]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t3 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[1]+'),' +l[4]+ ',(' +l[2]+ ',' +l[3]+ '))', schema = 'newick', taxon_set = M)
    t4 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[2]+'),' +l[1]+ ',(' +l[3]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t5 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[2]+'),' +l[3]+ ',(' +l[1]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t6 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[2]+'),' +l[4]+ ',(' +l[1]+ ',' +l[3]+ '))', schema = 'newick', taxon_set = M)
    t7 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[3]+'),' +l[1]+ ',(' +l[2]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t8 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[3]+'),' +l[2]+ ',(' +l[1]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t9 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[3]+'),' +l[4]+ ',(' +l[1]+ ',' +l[2]+ '))', schema = 'newick', taxon_set = M)
    t10 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[4]+'),' +l[1]+ ',(' +l[2]+ ',' +l[3]+ '))', schema = 'newick', taxon_set = M)
    t11 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[4]+'),' +l[2]+ ',(' +l[1]+ ',' +l[3]+ '))', schema = 'newick', taxon_set = M)
    t12 = dendropy.tree.get_from_string('((' +l[0]+ ',' +l[4]+'),' +l[3]+ ',(' +l[1]+ ',' +l[2]+ '))', schema = 'newick', taxon_set = M)
    t13 = dendropy.tree.get_from_string('((' +l[1]+ ',' +l[2]+'),' +l[0]+ ',(' +l[3]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t14 = dendropy.tree.get_from_string('((' +l[1]+ ',' +l[3]+'),' +l[0]+ ',(' +l[2]+ ',' +l[4]+ '))', schema = 'newick', taxon_set = M)
    t15 = dendropy.tree.get_from_string('((' +l[1]+ ',' +l[4]+'),' +l[0]+ ',(' +l[2]+ ',' +l[3]+ '))', schema = 'newick', taxon_set = M)
    counter = {t1:0, t2:0, t3:0, t4:0, t5:0, t6:0,t7:0,t8:0, t9:0, t10:0, t11:0, t12:0, t13:0, t14:0, t15:0}
    return counter

# I should merge this next with dist_counter after I test them both
def get_dist(l,treelist):
    dist = dist_counter(l, treelist)
    for i in range(len(treelist)):
        smalltree = treelist[i].retain_taxa_with_labels(l)
        trace = [smalltree.symmetric_difference(dist.keys()[j]) for j in range(15)]
        if trace.count(0) > 1:
            garbage.append(['error: too many zeroes', smalltree, treelist[i]])
        elif 0 in trace:
            dist[trace.index(0)] = dist[trace.index(0)] + 1
        else:
            garbage.append([smalltree, treelist[i]])
    return [dist, garbage]

#S is a rooted species tree, l is a list of five taxa labels like in get_dist, etc.Picking the edge and the quintet are unresolved
def score_quintet(S,l, treelist):
    quintet_tree = S.retain_taxa_with_labels(l)
    quintet_tree.ladderize(ascending=False)
    rooted_shape_str = quintet_tree.as_newick_string()
    rooted_shape = ''
    y = ['(', ')', ',']
    for i in range(len(rooted_shape_str)):
        if rooted_shape_str[i] in y:
            rooted_shape = rooted_shape + rooted_shape_str[i]
    shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15] = get_dist(l, treelist)[0].values() 
    scorefuncs= [inv51, inv52, inv53]
    if rooted_shape in shapes:
        score_fucntion = scorefuncs[shapes.index(rooted_shape)]
        score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    else:
        print 'error: quintet topology not in list'
    return score

