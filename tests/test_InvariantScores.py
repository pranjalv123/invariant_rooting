import unittest
import InvariantScores
import dendropy
import matrixmaker

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
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15 = 5,4,4,3,2,2,1,1,2,1,1,2,3,1,1
        score = InvariantScores.inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        self.assertEqual(score,0)

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
        hscore = info.quartet_score(['a','b','d','e'],0,1)
        self.assertEqual(escore,1)
        self.assertEqual(fscore,1)
        self.assertEqual(fscore2,4) 
        self.assertEqual(escore2,1) 
        self.assertEqual(escore3,4)
        self.assertEqual(gscore,1) 
        self.assertEqual(hscore,0)
        self.assertEqual(acdescore,1)

 
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
        self.assertEqual(fae,0)
        self.assertEqual(fbc,1)
        self.assertEqual(fcd,1)
        self.assertEqual(fce,1)

    def testTreeScore(self):
        info = InvariantScores.QuartetsInfo("output2.txt") 
        treescore = info.score_tree("lonetreequartets")
        treescore2 = info.score_tree("lonetree2quartets")
        self.assertEqual(treescore,4)
        self.assertEqual(treescore2,10)


class TestMatrixScoring(unittest.TestCase):
    def testScoreBig(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        M = matrix_maker.matrix()
        self.assertEqual(M.shape, (20,20))


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
        self.assertEqual(S1,2)
        self.assertEqual(S2,4)

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

class TestClassCompositions(unittest.TestCase):
    def testSubsetPenaltiesinstance(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        A = InvariantScores.SubsetPenalties(labels, setlist, 'output2.txt')
        M = matrixmaker.MatrixMaker(labels,setlist)
        Q = InvariantScores.QuartetsInfo('output2.txt')
        self.assertEqual(A.matrix[0,4], M.matrix()[0,4])
        self.assertEqual(Q.quartet_dict(), A.quartetsinfo.quartet_dict())
        self.assertEqual(A.labels, M.labels)
    

if __name__ == '__main__':
     unittest.main(
)
