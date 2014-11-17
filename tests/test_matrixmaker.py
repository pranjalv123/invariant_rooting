import unittest
import matrixmaker

class TestmatrixMaker(unittest.TestCase):
    def test_binlabels(self):
        labels = ['MOM', 'DAD', 'BOB']
        setlist = [['MOM', 'BOB'], ['DAD']]
        matrix_maker = matrixmaker.MatrixMaker(labels, setlist)
        bin1 = matrix_maker.binlabels["MOM"]
        bin2 = matrix_maker.binlabels['BOB']
        self.assertEqual(bin1, 2**0)
        self.assertEqual(bin2, 2**2)
       

    def test_setlistlabels(self):
        labels = ['MOM', 'DAD', 'BOB']
        setlist = [['MOM', 'BOB'], ['DAD']]
        matrix_maker = matrixmaker.MatrixMaker(labels, setlist)
        setbin1 = matrix_maker.binset(['MOM', 'BOB'])
        self.assertEqual(setbin1, 2**0 + 2**2)

    def test_containment(self):
        labels = ['MOM', 'DAD', 'BOB']
        setlist = [['MOM', 'BOB'], ['DAD'],['MOM', 'DAD', 'BOB']]
        matrix_maker = matrixmaker.MatrixMaker(labels, setlist)        
        setbin1 = matrix_maker.binset(['MOM', 'BOB'])
        setbin2 = matrix_maker.binset(['MOM', 'DAD', 'BOB'])
        setbin3 = matrix_maker.binset(['DAD'])
        test1 = matrix_maker.inclusion(['MOM', 'BOB'], ['MOM', 'DAD', 'BOB'])
        test2 = matrix_maker.inclusion(['DAD'], ['MOM', 'BOB'])        
        self.assertEqual(test1, 1)
        self.assertEqual(test2, 0)

    def test_make_matrix(self):
        labels = ['MOM', 'DAD', 'BOB']
        setlist = [['MOM', 'BOB'], ['DAD'],['MOM', 'DAD', 'BOB']]
        matrix_maker = matrixmaker.MatrixMaker(labels, setlist)        
        M = matrix_maker.matrix()
        self.assertEqual(M.shape, (3,3))
        self.assertEqual(M[2,0], 1)
        self.assertEqual(M[1,0], 0)
        self.assertEqual(M[2,1], 1)

    def test2_make_matrix(self):
        labels = ['a','b','c','d','e']
        setlist = [['a'],['b'],['c'],['d'],['e'], ['a','b'], ['a', 'c'], ['a','d'], ['a','e'], ['b','c'],['c', 'd'], ['c', 'e'], ['a','b', 'c'], ['a','b','e'],['b', 'c', 'd'], ['b','c','e'], ['a', 'b', 'c', 'd'],['a','b', 'c', 'e'], ['b', 'c', 'd', 'e'], ['a', 'b', 'c', 'd', 'e']]
        matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        M = matrix_maker.matrix()
        self.assertEqual(M.shape, (20,20))
        self.assertEqual(M[4,5], 0)
        self.assertEqual(M[5,0], 1)
        self.assertEqual(M[5,1], 1)
        self.assertEqual(M[16,5], 1)
        self.assertEqual(M[16,6], 1)
        self.assertEqual(M[16,7], 1)
        self.assertEqual(M[16,8], 0)
        self.assertEqual(M[16,9], 1)
        self.assertEqual(M[16,12],1)
        self.assertEqual(M[16,13],0)



if __name__ == '__main__':
    unittest.main()
