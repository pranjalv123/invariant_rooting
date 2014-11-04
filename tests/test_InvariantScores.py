import unittest
import InvariantScores
import dendropy

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

if __name__ == '__main__':
    unittest.main(
)
