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

#if __name__ == '__main__':
    #unittest.main()
