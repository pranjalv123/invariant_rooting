import unittest
import InvariantScores
import dendropy
import matrixmaker
import pprint
import numpy as np
import copy 

class TestInvariantScores3(unittest.TestCase):
    def test_inv3(self):
        u1, u2, u3 = 2, 1, 1 #when u1 > u2 = u3 everything is right
        score = InvariantScores.inv3(u1,u2,u3)
        self.assertEqual(score, 0)

    def test_inv3_u2_ne_u3(self):
        u1, u2, u3 = 3,2,1
        score = InvariantScores.inv3(u1,u2,u3)
        self.assertEqual(score,1)

    def test_inv3_u1_lt_u2(self):
        u1, u2, u3 = 2,3,1
        score = InvariantScores.inv3(u1,u2,u3)
        self.assertEqual(score,3)

    def test_inv3_u1_lt_u3(self):
        u1, u2, u3 = 2,1,3
        score = InvariantScores.inv3(u1,u2,u3)
        self.assertEqual(score,3)

class TestInvariantScores5(unittest.TestCase):
    def test_inv51(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 43, 2, 1, 0, 3, 11, 9, 3, 3, 8, 10, 50, 23, 14, 7
        score = InvariantScores.inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,10)

    def test_inv51_2(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 30, 2, 1, 10, 3, 12, 9, 3, 0, 4, 4, 13, 23, 7, 7
        score = InvariantScores.inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,7)


    def test_inv52(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 43, 2, 1, 0, 3, 11, 9, 3, 3, 8, 10, 50, 23, 14, 7
        score = InvariantScores.inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,21)

    def test_inv52_2(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 30, 2, 1, 10, 3, 12, 9, 3, 0, 4, 4, 13, 23, 7, 7
        score = InvariantScores.inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,19)

    def test_inv53(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 43, 2, 1, 0, 3, 11, 9, 3, 3, 8, 10, 50, 23, 14, 7
        score = InvariantScores.inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,4)

    def test_inv53_2(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 30, 2, 1, 10, 3, 12, 9, 3, 0, 4, 4, 13, 23, 7, 7
        score = InvariantScores.inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,1)




'''
    def test_inv51(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 5,4,4,3,2,2,1,1,2,1,1,2,3,1,1
        score = InvariantScores.inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,0)

    def test_inv52(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 52, 42, 292, 16, 7, 7, 36, 15, 30, 185, 14, 34, 45, 46, 91
        score = InvariantScores.inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,709)

    def test_inv53(self):
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 52, 42, 292, 16, 7, 7, 36, 15, 30, 185, 14, 34, 45, 46, 91
        score = InvariantScores.inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,760)

'''

class TestUnrootedDistributions(unittest.TestCase):
    def test_dist_counter(self):
        l = ["A", "B", "C", "D", "E"]
        tree1 = dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        tree2 = dendropy.Tree.get_from_string('(((A,B),C),D,(E,(F,G)))', 'newick', taxon_set = tree1.taxon_set)
        tree3 = dendropy.Tree.get_from_string('(((A,B),C),(D,E),(F,G))', 'newick', taxon_set = tree1.taxon_set)
        treelist = dendropy.TreeList([tree1, tree2, tree3])
        counter = InvariantScores.dist_counter(l,treelist)
        h = []
        for i in range(15):
            h.append(counter[0][i].as_newick_string())
        self.assertEqual(h,['((A,B),C,(D,E))', '((A,B),D,(C,E))', '((A,B),E,(C,D))', '((A,C),B,(D,E))', '((A,C),D,(B,E))', '((A,C),E,(B,D))', '((A,D),B,(C,E))', '((A,D),C,(B,E))', '((A,D),E,(B,C))', '((A,E),B,(C,D))', '((A,E),C,(B,D))', '((A,E),D,(B,C))', '((B,C),A,(D,E))', '((B,D),A,(C,E))', '((B,E),A,(C,D))'])

    def test_get_dist(self):
        l = ["A", "B", "C", "D", "E"]
        tree1 = dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        tree2 = dendropy.Tree.get_from_string('(((A,B),C),D,(E,(F,G)))', 'newick', taxon_set = tree1.taxon_set)
        tree3 = dendropy.Tree.get_from_string('(((A,B),C),(D,E),(F,G))', 'newick', taxon_set = tree1.taxon_set)
        treelist = dendropy.TreeList([tree1, tree2, tree3])
        X = InvariantScores.get_dist(l,treelist)
        u = [X[0][0][i].as_newick_string() for i in range(15)]
        self.assertEqual(X[0][1], [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(X[1], [])
        self.assertEqual(u, ['((A,B),C,(D,E))', '((A,B),D,(C,E))', '((A,B),E,(C,D))', '((A,C),B,(D,E))', '((A,C),D,(B,E))', '((A,C),E,(B,D))', '((A,D),B,(C,E))', '((A,D),C,(B,E))', '((A,D),E,(B,C))', '((A,E),B,(C,D))', '((A,E),C,(B,D))', '((A,E),D,(B,C))', '((B,C),A,(D,E))', '((B,D),A,(C,E))', '((B,E),A,(C,D))'])


    def test_basic_score_quintet(self):
        tree1 = dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        tree2 = dendropy.Tree.get_from_string('(((A,B),C),D,(E,(F,G)))', 'newick', taxon_set = tree1.taxon_set)
        tree3 = dendropy.Tree.get_from_string('(((A,B),C),(D,E),(F,G))', 'newick', taxon_set = tree1.taxon_set)
        treelist = dendropy.TreeList([tree1, tree2, tree3])
        S = dendropy.Tree.get_from_string("[&R] ((((A,B), C), D), ((E,F), (G,H)));", 'newick', taxon_set = treelist.taxon_set)
        l = ["A", "B", "C", "D", "E"]
        W = [InvariantScores.get_rooted_quintet(S,l,i).as_newick_string() for i in range(13)]
        V = [InvariantScores.get_rooted_quintet(S,l,i).as_newick_string() for i in range(14)]
        score = InvariantScores.basic_score_quintet(S,l,treelist)
        self.assertEqual(score,0) 
        self.assertEqual(W,['((((D,E),C),B),A)', '((((D,E),C),A),B)', '(((D,E),C),(A,B))', '(((A,B),(D,E)),C)', '(((A,B),C),(D,E))', '((((A,B),C),E),D)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)'])
        self.assertEqual(V,['((((D,E),C),B),A)', '((((D,E),C),A),B)', '(((D,E),C),(A,B))', '(((A,B),(D,E)),C)', '(((A,B),C),(D,E))', '((((A,B),C),E),D)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)', '((((A,B),C),D),E)'])
    
#def test_inv52(self):
        #u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 

class TestQuartetStuff(unittest.TestCase):
    def testQuartetDict(self):
        info = InvariantScores.QuartetsInfo("output2.txt")
        d = info.quartet_dict()
        self.assertEqual(d['((a,b),(c,d));'], 2)
        self.assertEqual(d['((a,b),(c,e));'], 2)
        self.assertEqual(d['((a,b),(d,e));'], 3)

    def test2QuartetDict(self):
        info = InvariantScores.QuartetsInfo("output2.txt")
        d = info.quartet_dict()
        e = info.get_freqs(['a', 'b', 'c', 'd'])
        e1 = info.get_freqs(['a', 'b', 'c', 'e'])
        e2 = info.get_freqs(['a', 'b', 'd', 'e'])
        e3 = info.get_freqs(['a', 'c', 'd', 'e'])
        e4 = info.get_freqs(['b', 'c', 'd', 'e'])
        self.assertEqual(e, [2, 1, 0])
        self.assertEqual(e1, [2, 1, 0])
        self.assertEqual(e2, [3, 0, 0])
        self.assertEqual(e3, [2, 0, 1])
        self.assertEqual(e4, [2, 0, 1])

    def testQuartetScores(self):
        info = InvariantScores.QuartetsInfo("output2.txt")
        d = info.quartet_dict()
        e = info.get_freqs(['a', 'b', 'c', 'd'])
        f = info.get_freqs(['b', 'c', 'd', 'e'])
        escore = info.quartet_score(['a', 'b', 'c', 'd'],0,1)
        escore2 = info.quartet_score(['a', 'b', 'c', 'd'],2,3)
        escore3 = info.quartet_score(['a', 'b', 'c', 'd'],0,3)
        fscore = info.quartet_score(['b', 'c', 'd', 'e'],0,1)
        fscore2 = info.quartet_score(['b', 'c', 'd', 'e'],0,2)
        gscore = info.quartet_score(['a','b','c','e'],0,1)
        acdescore = info.quartet_score(['a','c','d','e'],0,1)
        aecdscore = info.quartet_score(['a','c','d','e'],0,3)
        hscore = info.quartet_score(['a','b','d','e'],0,1)
        becdscore = info.quartet_score(['b', 'c', 'd', 'e'],0,3)
        self.assertEqual(escore,1)
        self.assertEqual(fscore,1)
        self.assertEqual(fscore2,4) 
        self.assertEqual(escore2,1) 
        self.assertEqual(escore3,4)
        self.assertEqual(gscore,1) 
        self.assertEqual(hscore,0)
        self.assertEqual(acdescore,1)
        self.assertEqual(becdscore,3)
        self.assertEqual(aecdscore,3)

    def testQuartetScores2(self):
        info = InvariantScores.QuartetsInfo("R1.10gene.tre.quartets")
        score1 = info.quartet_score(['S1', 'S10', 'S2', 'S3'],0,1)
        score2 = info.quartet_score(['S2', 'S10', 'S1', 'S3'],0,1)
        self.assertEqual(score1,4)
        self.assertEqual(score2,2)
    
    def testQuartetScores3(self):
        info = InvariantScores.QuartetsInfo("11treesquartets")
        score1 = info.quartet_score(['S1','S2','S4','S11'],0,3)
        self.assertEqual(score1,3)
        score2 = info.quartet_score(['S2','S5','S1','S11'],2,3)
        self.assertEqual(score2,1)
        score3 = info.quartet_score(['S2','S5','S1','S11'],0,1) 
        self.assertEqual(score3,1)
        score4 = info.quartet_score(['S1','S2','S5','S11'],0,3)
        self.assertEqual(score4,1)
        score5 = info.quartet_score(['S1','S11','S4','S5'],0,1)
        self.assertEqual(score5,3) 
        score6 = info.quartet_score(['S1','S4','S2','S5'],2,3)
        self.assertEqual(score6,0)
        score7 = info.quartet_score(['S11','S4','S2','S5'],0,1)
        self.assertEqual(score7,1)


class TestQuartetLabelsDict(unittest.TestCase):
    def testlabelsdict(self):
        info = InvariantScores.QuartetsInfo("output2.txt")
        d = info.quartet_dict()
        g = info.quartet_labels_dict()
        self.assertEqual(g['((a,b),(c,d));'], ['a','b','c','d'])
        self.assertEqual(g['((a,b),(c,e));'], ['a','b','c','e'])
        self.assertEqual(g['((a,b),(d,e));'], ['a','b','d','e'])

        
class TestPenaltyFunctions(unittest.TestCase):
    def testPenaltyDoubleton(self):
        info = InvariantScores.QuartetsInfo("output2.txt")
        #print info.quartet_dict()
        #print info.quartet_labels_dict()
        fab = info.score_double('a', 'b')
        fac = info.score_double('a', 'c')
        fad = info.score_double('a', 'd')
        fae = info.score_double('a', 'e')
        fbc = info.score_double('b', 'c')
        fcd = info.score_double('c', 'd')
        fce = info.score_double('c', 'e')
        self.assertEqual(fab,2)
        self.assertEqual(fac,7)
        self.assertEqual(fad,0)
        self.assertEqual(fae,3)
        self.assertEqual(fbc,1)
        self.assertEqual(fcd,7)
        self.assertEqual(fce,1)
    
    def testPenaltyDoubleton2(self):
        info = InvariantScores.QuartetsInfo("11treesquartets")
        #print info.quartet_dict()
        #print info.quartet_labels_dict()
        f = info.score_double('S11','S1')
        self.assertEqual(f,7)
        g = info.score_double('S1','S11')
        self.assertEqual(g,7)
        h = info.score_double('S1','S2')
        self.assertEqual(h,0) 
 
    def testTreeScore(self):
        info = InvariantScores.QuartetsInfo("output2.txt") 
        treescore = info.score_tree("lonetreequartets")
        treescore2 = info.score_tree("lonetree2quartets")
        self.assertEqual(treescore,4)
        self.assertEqual(treescore2,10)

    def testTreeScore2(self):
        info = InvariantScores.QuartetsInfo("R1.10gene.tre.quartets")

    def testTreeScore3(self):
        info = InvariantScores.QuartetsInfo('11treesquartets')
        treescore = info.score_tree('11trees-1-quartets')
        self.assertEqual(treescore,8)


class TestMatrixScoring(unittest.TestCase):
    def testMatrixMaker(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        M = matrix_maker.matrix()
        self.assertEqual(M.shape, (20,20))

    def testMatrixScoring1(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        M = matrixmaker.MatrixMaker(labels,setlist).matrix()
        SM = A.scored_matrix()
        #print M
        pp = pprint.PrettyPrinter(indent=4,width=300)
        #pp.pprint(SM)
        self.assertEqual(SM[5,0],2)
        self.assertEqual(SM[10,2],7)
        self.assertEqual(SM[12,5],4)
    
    def testMatrixScoring2(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['b', 'd'],['b','e'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        M = matrixmaker.MatrixMaker(labels,setlist).matrix()
        SM = A.scored_matrix()
        #print M
        pp = pprint.PrettyPrinter(indent=4,width=300)
        #pp.pprint(SM)
        self.assertEqual(SM[5,0],2)
        self.assertEqual(SM[12,2],7)
        self.assertEqual(SM[14,5],4)

    def testMatrixScoring3(self):
        labels  = ['a','b','c','d']
        setlist = [['a'],['b'],['c'],['d'],['a','b'], ['a', 'c'], ['a','d'], ['b','c'],['b', 'd'], ['c', 'd'], ['a','b', 'c'], ['a','b', 'd'], ['a','c','d'], ['b','c','d'],['a', 'b', 'c', 'd']];
        A = InvariantScores.SubsetPenalties(labels,setlist,'output3')
        SM = A.scored_matrix()
        #print SM
        self.assertEqual(SM[3,3],0)
        self.assertEqual(SM[4,0],8)
        self.assertEqual(SM[7,2],2)
        self.assertEqual(SM[13,1],8)
        self.assertEqual(SM[12,5],6)
        self.assertEqual(SM[10,3],np.inf)
        self.assertEqual(SM[13,4],np.inf)
        
class TestTracebacks(unittest.TestCase):
    def testTraceback1(self):
        labels  = ['a','b','c','d']
        setlist = [['a'],['b'],['c'],['d'],['a','b'], ['a', 'c'], ['a','d'], ['b','c'],['b', 'd'], ['c', 'd'], ['a','b', 'c'], ['a','b', 'd'], ['a','c','d'], ['b','c','d'],['a', 'b', 'c', 'd']];
        A = InvariantScores.SubsetPenalties(labels,setlist,'output3')
        SM = A.scored_matrix()
        cladelist = A.clades()
        #print cladelist
        #L = [['a'],['b'],['c'],['d'],['c', 'd'],['b','c','d'], ['a', 'b', 'c', 'd']]
        L = [['a','b','c','d'],['a'],['b','c','d'],['d'],['b','c'], ['b'],['c']]
        self.assertEqual(L,cladelist)

class TestTreesFromClades(unittest.TestCase):
    def testBuildTree1(self):
        labels  = ['a','b','c','d']
        setlist = [['a'],['b'],['c'],['d'],['a','b'], ['a', 'c'], ['a','d'], ['b','c'],['b', 'd'], ['c', 'd'], ['a','b', 'c'], ['a','b', 'd'], ['a','c','d'], ['b','c','d'],['a', 'b', 'c', 'd']];
        A = InvariantScores.SubsetPenalties(labels,setlist,'output3')
        SM = A.scored_matrix()
        cladelist = A.clades()
        T = A.tree()
        #print T
        self.assertEqual(T, '(a,((b,c),d));')

    def testBuildTree2(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['b', 'd'],['b','e'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        T = A.tree()
        #print T
 
    #def testScoreDoubles(self):
        #labels = ['a','b','c','d','e']
        #setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e']]
        #matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        #M = matrix_maker.matrix()
        #info = InvariantScores.QuartetsInfo("output2.txt")
        #self.assertEqual(M.shape, (12,12))
        #SM = InvariantScores.ScoredMatrix(labels,setlist,"output2.txt").scored_matrix()
        #SM = ScoredClass.scored_matrix()
        #self.assertEqual(SM[5,0], 2)
        #self.assertEqual(SM[5.1], 2)
        #self.assertEqual(SM[6,0], 7)
        #self.assertEqual(SM[6,2], 7)
        #self.assertEqual(SM[8,0], 0)
        #self.assertEqual(SM[8,4], 0)
        #self.assertEqual(SM[11,2], 1)
        #self.assertEqual(SM[11,4], 1)

class TestSubsetPenaltyScores(unittest.TestCase):
    def testPenalty1(self):        
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        scoring = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        S1 = scoring.penalty_score(['a'],['b'],0,0)
        S2 = scoring.penalty_score(['a', 'b'], ['c'],2,0)
        S3 = scoring.penalty_score(['a', 'b'], ['c','d'],2,7)
        S4 = scoring.penalty_score(['c'],['d'],0,0)
        self.assertEqual(S1,2)
        self.assertEqual(S2,4)
        self.assertEqual(S3,8)
        self.assertEqual(S4,7)
        

    def testQuartetSets(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        scoring = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        AQ1 = scoring.add_quartets(['a', 'b'], ['c'])
        SQ1 = scoring.subtract_quartets(['a', 'b'], ['c'])
        SQ2 = scoring.subtract_quartets(['a', 'b'], ['c','d'])
        self.assertEqual(AQ1, [[['a','c','d','e'],0,1], [['b','c','d','e'],0,1]])
        self.assertEqual(SQ1, [])
        self.assertEqual(SQ2, [[['a','b','c','d'],0,1]])

#matrix and scored_matrix should NOT be the same.....
class TestClassCompositions(unittest.TestCase):
    def testSubsetPenaltiesinstance(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        M = matrixmaker.MatrixMaker(labels,setlist)
        Q = InvariantScores.QuartetsInfo('output2.txt')
        self.assertEqual(A.matrix[0,4], M.matrix()[0,4])
        self.assertEqual(A.matrix.all(), M.matrix().all())
        self.assertEqual(Q.quartet_dict(), A.quartetsinfo.quartet_dict())
        self.assertEqual(A.labels, M.labels)
        #self.assertEqual(A.matrix.all(), A.scored_matrix().all())
        #self.assertEqual(A.scored_matrix()[5,0],1) 

    def testClassStructure(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        M = matrixmaker.MatrixMaker(labels,setlist)
        #B = A.pairs()
        #print B
        #self.assertTrue([['a', 'b'], ['c', 'd']] in B)
        #self.assertTrue([['a', 'b'], ['c']] in B)
        self.assertEqual(setlist, A.setlist) 
         

class TestPowerSet(unittest.TestCase):
    def testPowerSet1(self):
        labellist = ['ana', 'ann', 'bob','sue']
        B = [['ana'], ['ann'], ['bob'], ['sue'], ['ana', 'ann'], ['ana', 'bob'], ['ana', 'sue'], ['ann', 'bob'], ['ann', 'sue'], ['bob', 'sue'], ['ana', 'ann', 'bob'], ['ana', 'ann', 'sue'], ['ana', 'bob', 'sue'], ['ann', 'bob', 'sue'], ['ana', 'ann', 'bob', 'sue']]
        C = InvariantScores.powerset(labellist)
        #print C
        self.assertEqual(B,C)

class TestQuintetTopologies(unittest.TestCase):
    def testQuintetsForOneSet(self):
        labels = ['a','b','c','d','e']
        Q = InvariantScores.QuintetsMaker(labels)
        q = Q.quintets(('a','b','c','d','e'))
        l = len(q)
        self.assertEqual(15,l)
        self.assertEqual('((a,b),c,(d,e));\n',q[0])
        self.assertEqual('((b,e),a,(c,d));\n',q[14])
        self.assertEqual('((a,e),b,(c,d));\n',q[9])

    
    def testFives(self):
        labels = ['a','b','c','d','e','f']
        Q = InvariantScores.QuintetsMaker(labels)
        s = Q.fives()
        self.assertEqual(len(s),6)
        self.assertEqual(('a','b','c','d','e'),s[0])
        self.assertEqual(('b','c','d','e','f'),s[5])
        self.assertEqual(('a', 'b', 'd', 'e', 'f'),s[3])



    def testAllQuintets1(self):
        labels = ['a','b','c','d','e','f']
        Q = InvariantScores.QuintetsMaker(labels)
        q = Q.allquintets()
        s = Q.fives()
        self.assertEqual(len(s),6)
        self.assertEqual(len(q),90)
        self.assertEqual(labels, Q.labels)
        #print q

    def testMakeFile(self):
        labels = ['a','b','c','d','e','f']
        Q = InvariantScores.QuintetsMaker(labels)
        Q.makefile('testfile')
        g = open('testfile','r')
        t =[]
        for k in range(20):
            t.append(g.readline())
        self.assertEqual(t[0],'((a,b),c,(d,e));\n')
        self.assertEqual(t[3],'((a,c),b,(d,e));\n')
        self.assertEqual(t[15],'((a,b),c,(d,f));\n')
        self.assertEqual(t[19],'((a,c),d,(b,f));\n')



class TestCladesFromFile(unittest.TestCase):
    def testClades1(self):
        C = InvariantScores.FileClades('fakeclusters')
        self.assertEqual(C.filename,'fakeclusters')
        #print C.clades()
        #self.assertEqual(C.clades()[-1], ['a','b','c','d'])
        self.assertEqual(C.clades()[0], ['d'])
        self.assertEqual(len(C.clades()),12)
        self.assertEqual(C.clades()[4],['c','d'])
 
    def testClades2(self):
        E = InvariantScores.FileClades('songclusters')
        E1 = InvariantScores.FileClades('songclusters1')
        E2 = InvariantScores.FileClades('songclusters2')
        E3 = InvariantScores.FileClades('songclusters3')
        E4 = InvariantScores.FileClades('songclusters4')
        C = E.clades()
        C1 = E1.clades()
        C2 = E2.clades()
        C3 = E3.clades()
        C4 = E4.clades()
        self.assertEqual(len(C1),1933)
        self.assertEqual(len(C2),1933)
        self.assertEqual(len(C3),1933)
        self.assertEqual(len(C4),1933)
        #print len(C
        self.assertEqual(len(C),1933)
        D = [C[i] for i in range(10)]
        #print D
        F = [C[-i] for i in range(1,11)]
        #print F
        self.assertEqual(C,C1)
        self.assertEqual(C,C2)
        self.assertEqual(C,C3)
        self.assertEqual(C,C4)
    
    def testClades3(self):
        F = InvariantScores.FileClades('mixedupclusters')
        D = F.clades()
        #print len(D)
        #print D
        self.assertEqual(len(D),12)
        self.assertEqual(D[-1],['a', 'b', 'c', 'd'])
        self.assertEqual(D[2],['a'])
        self.assertEqual(D[7],['a','d'])

    def testClades4(self):
        F = InvariantScores.FileClades('Ruth_weakILS_astralclades')
        D = F.clades()
        #print len(D)
        #print D
        self.assertEqual(len(D),555)
        W = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
        self.assertEqual(W,D[-1])
        self.assertTrue(['9'] in D)

class TestSubsetPenaltiesInequalitiesOnly(unittest.TestCase):
    def testNewClassPropertiesInequalitiesOnly(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        B = InvariantScores.SubsetPenaltiesInequalitiesOnly(labels, setlist, 'output2.txt')
        M = matrixmaker.MatrixMaker(labels,setlist)
        Q = InvariantScores.QuartetsInfo('output2.txt')
        P = InvariantScores.QuartetsInfoInequalitiesOnly('output2.txt')
        self.assertEqual(B.matrix[0,4], M.matrix()[0,4])
        self.assertEqual(B.matrix.all(), M.matrix().all())
        self.assertEqual(Q.quartet_dict(), B.quartetsinfo.quartet_dict())
        self.assertEqual(A.labels, B.labels)
        self.assertEqual(P.quartet_dict(), Q.quartet_dict())
        #self.assertEqual(B.scored_matrix()[-1,-2], A.scored_matrix()[-1,-2])
        #print A.scored_matrix()
        #print B.scored_matrix()


class TestTreesFromCladesInequalitiesOnly(unittest.TestCase):
    def testBuildTree1(self):
        labels  = ['a','b','c','d']
        setlist = [['a'],['b'],['c'],['d'],['a','b'], ['a', 'c'], ['a','d'], ['b','c'],['b', 'd'], ['c', 'd'], ['a','b', 'c'], ['a','b', 'd'], ['a','c','d'], ['b','c','d'],['a', 'b', 'c', 'd']];
        A = InvariantScores.SubsetPenaltiesInequalitiesOnly(labels,setlist,'output3')
        SM = A.scored_matrix()
        #print SM
        cladelist = A.clades()
        T = A.tree()
        #print T
        #self.assertEqual(T, '(a,((b,c),d));')

    def testBuildTree2(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['b', 'd'],['b','e'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenaltiesInequalitiesOnly(labels, setlist, 'output2.txt')
        SM = A.scored_matrix()
        #print SM
        T = A.tree()
        #print T
 
class TestDendropyTrees(unittest.TestCase):
    def testBuildDendropyTree1(self):
        Taxa = ['1','2','3','4']
        cladelist = [['1','2','3'],['2','3'],['1'],['2'],['3']]
        T = InvariantScores.dendropy_clades_tree(Taxa,cladelist)
        print T
        S = dendropy.Tree.get_from_string('(1,4,(2,3));','newick',taxon_set = T.taxon_set)
        d = S.symmetric_difference(T)
        self.assertEqual(d,0)

    def testBuildDendropyTree2(self):
        Taxa = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
        cladelist =[['1', '10', '11', '2', '3', '4', '5', '6', '7', '8'],['9'], ['5', '6', '7', '8'], ['1', '10', '11', '2', '3', '4'], ['7', '8'], ['5', '6'], ['10'], ['1', '11', '2', '3', '4'], ['7'], ['8'], ['5'], ['6'], ['11'], ['1', '2', '3', '4'], ['4'], ['1', '2', '3'], ['3'], ['1', '2'], ['1'], ['2']]
        T = InvariantScores.dendropy_clades_tree(Taxa,cladelist)
        print T
        S = dendropy.Tree.get_from_string('((((((1,2),3),4),11),10),((5,6),(7,8)),9);','newick',taxon_set = T.taxon_set)
        print S
        S1 = dendropy.Tree.get_from_string('((((((1,2),4),3),11),10),((5,6),(7,8)),9);','newick',taxon_set = T.taxon_set)
        print S1
        d = S.symmetric_difference(T)
        d1 = S1.symmetric_difference(T)
        self.assertEqual(d,0)
        self.assertEqual(d1,2)
        
if __name__ == '__main__':
     unittest.main(
)


############### Unit test Kajori ############### 
# START Unit testing for checking the following :-
#  1) Check if split_hash_bitmask_changes with deep copying -PASSED
#  2) Check if edge.oid with copying -FAILED
#  3)checking the taxons associated with the splits_hash_bitmask do not change by copying -PASSED
#  4) the splits_hash_bitmask associated with an edge donot change with copying
class TestSplitHashBitMask(unittest.TestCase):
    def test_split_hash_bitmask(self):
        S_old = dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        S_old.deroot()
        S_old.encode_splits()
        S_old.update_splits()
        DS_old = S_old.split_edges
        #print 'TestSplitHashBitMask'
        #print 'OLD split_bitmasks are'
        #print DS_old.keys()

        S_new=copy.deepcopy(S_old)
        S_new.deroot()
        S_new.encode_splits()
        S_new.update_splits()
        DS_new = S_new.split_edges
        #print 'NEW split_bitmasks are'
        #print DS_new.keys()

        #checking assertions
        self.assertEqual(len(set(DS_new.keys())-set(DS_old.keys())),0) # this checks DS_old.keys has all the elements of DS_new.keys
        self.assertEqual(len(set(DS_old.keys())-set(DS_new.keys())),0) # this checks DS_new.keys has all the elements of DS_old.keys

        #convert the keys to a list 
        key_old=[ k for k in DS_old.keys()]
        key_new=[ k for k in DS_new.keys()]
        self.assertEqual(len(key_old),len(key_new))
        # checking the taxons associated with the splits_hash_bitmask do not change by copying
        for index in range(len(key_old)):
            #print index
            taxon_set_old=InvariantScores.taxon_with_split_edge(S_old,key_old[index])
            #print taxon_set_old
            taxon_set_new=InvariantScores.taxon_with_split_edge(S_new,key_new[index])
            #print taxon_set_new
            self.assertEqual(set(taxon_set_old),set(taxon_set_new))
            self.assertEqual(len(taxon_set_old),len(taxon_set_new))
            self.assertEqual(key_old[index],key_new[index])
            
        
        '''
        edgeset_old= [e for e in S_old.preorder_edge_iter()]
        edgeset_new= [e for e in S_new.preorder_edge_iter()]
        for index in range(len(edgeset_old)):
            assert edgeset_old[index].oid==edgeset_new[index].oid #this assertion fails
            print index,edgeset_old[index],edgeset_new[index]
        '''
# START Module testing for InvariantScores.taxon_with_split_edge 
# taxon_with_split_edge takes as input a split_hash_bit_mask and returns taxons on on side of the bitmask
class Test_Taxon_With_Split_Edge(unittest.TestCase):
    def test_split_edge(self):
        S= dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        S.encode_splits()
        S.update_splits()
        DS = S.split_edges

        ES = [e for e in S.postorder_edge_iter()]
        INVDS = {v: k for k, v in DS.items()}

        #dict that gives bitmasks in same order as edge list should be
        ESBITS = [INVDS[ES[i]] for i in range(len(ES))]
        #S.print_plot()

        output=[[]]
        for k in  ESBITS:
            output.append(InvariantScores.taxon_with_split_edge(S,k))
        true_output=[[], ['A'], ['B'], ['A', 'B'], ['C'], ['D'], ['E'], ['D', 'E'], ['C', 'D', 'E'], ['F'], ['G'], ['F', 'G'], ['A', 'B', 'C', 'F', 'G', 'D', 'E']]
        self.assertEqual(output,true_output)
# End Module testing for InvariantScores.taxon_with_split_edge
 
# START Module testing for InvariantScores.total_quintet_score_distance_kajori
# total_quintet_score_distance_kajori takes as input an edge e and returns taxons closest to e
class Test_Total_Quintet_Score_Distance(unittest.TestCase):
    
    def test_nearest_quintet_kajori(self):
        S= dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        ES = [e for e in S.postorder_edge_iter()]
        output=[[]]
        for i in range(len(ES)):
            taxon_1,taxon_2=InvariantScores.nearest_quintet_kajori(S,i)
            #print 'list(set(taxon_1) & set(taxon_2)',list(set(taxon_1) & set(taxon_2))
            self.assertEqual(list(set(taxon_1) & set(taxon_2)),[]) #intersection of taxon_1 and taxon_2 is empty
            self.assertEqual(set(taxon_1) | set(taxon_2),set(['A', 'B', 'C', 'F', 'G', 'D', 'E'])) #intersection of taxon_1 and taxon_2 is empty
        #self.assertEqual(output,true_output)
        
    def test_sort_list_distance_kajori(self):
        S= dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        ES = [e for e in S.postorder_edge_iter()]
        taxon_1,taxon_2=InvariantScores.nearest_quintet_kajori(S,10)
        self.assertEqual(set(taxon_2),set(['A', 'B', 'C', 'D', 'E'])) 
        S.reroot_at_edge(ES[10])
        quintet=InvariantScores.sort_list_distance_kajori(S,taxon_2)
        self.assertEqual(quintet,['A', 'B', 'C', 'D', 'E']) 
    
    def test_find_edge_associated_with_outlier_kajori(self):
        S= dendropy.Tree.get_from_string('(((A,B),(C,(D,E)),(F,G)),H)', 'newick')
        pos=InvariantScores.find_edge_associated_with_outlier_kajori(S,'H')
        #print 'test_find_edge_associated_with_outlier_kajori',pos
        taxon_1,taxon_2=InvariantScores.nearest_quintet_kajori(S,pos)
        if (taxon_2==['H'] or taxon_1==['H'] ):
            print taxon_1,taxon_2
            self.assertTrue(False)
    
    
    
    def test_remove_split_bitmask_kajori(self):
        S = dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        S.deroot()
        S.encode_splits()
        S.update_splits()
        DS = S.split_edges
        ES = [e for e in S.postorder_edge_iter()]
        INVDS = {v: k for k, v in DS.items()}
        #print ' test_remove_split_bitmask_kajori DS keys =',DS.keys()
        #dict that gives bitmasks in same order as edge list should be
        ESBITS = [INVDS[ES[i]] for i in range(len(ES))]
        #true_output=[[], ['A'], ['B'], ['A', 'B'], ['C'], ['D'], ['E'], ['D', 'E'], ['C', 'D', 'E'], ['F'], ['G'], ['F', 'G'], ['A', 'B', 'C', 'F', 'G', 'D', 'E']]
        
        for k in range(len(ES)-1):
            #print 'ESBITS[k]',ESBITS[k]
            taxon_1=InvariantScores.taxon_with_split_edge(S,ESBITS[k])
        
            S_new=InvariantScores.remove_split_bitmask_kajori(copy.deepcopy(S),ESBITS[k])
            taxon_2=[ n.taxon.label for n in S_new.leaf_nodes()  ]
            
            #print ' taxon_1 , taxon_2 ',taxon_1,taxon_2
            self.assertEqual(len(set(taxon_1) & set(taxon_2)),0) 
            self.assertEqual(set(taxon_1) | set(taxon_2),set(['A', 'B', 'C', 'F', 'G', 'D', 'E'])) 
    
    
    def test_total_quintet_score_distance_kajori(self):
        S= dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        tree1 = dendropy.Tree.get_from_string('((A,B),(C,(D,E)),(F,G))', 'newick')
        tree2 = dendropy.Tree.get_from_string('(((A,B),C),D,(E,(F,G)))', 'newick', taxon_set = tree1.taxon_set)
        tree3 = dendropy.Tree.get_from_string('(((A,B),C),(D,E),(F,G))', 'newick', taxon_set = tree1.taxon_set)
        treelist = dendropy.TreeList([tree1, tree2, tree3])
        
        ES = [e for e in S.postorder_edge_iter()]
       
        quintet=InvariantScores.total_quintet_score_distance_kajori(S,2,treelist)
        print quintet
        #self.assertEqual(quintet,['A', 'B', 'C', 'D', 'E']) 
    
    
   
# End Module testing for InvariantScores.taxon_with_split_edge
