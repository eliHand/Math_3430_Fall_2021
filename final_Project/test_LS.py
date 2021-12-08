import unittest
from LS import *

matrixNotes = [[2,2,1],[-2,1,2],[1,3,1]]
matrixTest = [[1,2],[3,4]]
matrixTest2 = [[4]]
LSTest = [[1,0,-1,0],[0,1,0,-1],[1,1,1,1]]
b = [0,1,3,4]



matrix1 = [[1,2,3],[4,5,6],[7,8,9]]
matrix2OLD = [[6,8,4],[3,5,1],[7,2,9]]
matrix2 = [[1,2,3],[4,5,6],[7,2,9]]
matrix3 = [[4,5,6],[3,5,8],[2,4,1]]
matrixCJ = [[1,2,3],[4,5,6]]
test_vector_01 = [1, 2, 4]
test_vector_02 = [3, 1, 2]
test_vector_03 = [6,complex(3,2),7]
n = 9

class TestLeastSquares(unittest.TestCase):
    def test_leastSquares(self):
        self.assertAlmostEqual(leastSquares(matrixNotes,test_vector_01),[3.333333333333333, 2.6666666666666665, 1.666666666666667])
        self.assertAlmostEqual(leastSquares(matrix2,test_vector_01),[4.543441112511214, -0.4364357804719846, -0.40824829046386335])

if __name__ == "__main__":
    unittest.main()
