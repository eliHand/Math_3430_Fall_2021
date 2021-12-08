import unittest
from LA import *


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

class TestVectorOperations(unittest.TestCase):

    def test_addVectors(self):
        self.assertEqual(add_vectors(test_vector_01,test_vector_02),[4,3,6])
        self.assertEqual(add_vectors(test_vector_01,test_vector_03),[7,5+2j,11])
    
    def test_scalarVecotr(self):
        self.assertEqual(scalar_vector(test_vector_01,n),[9,18,36])
        self.assertEqual(scalar_vector(test_vector_02,n),[27,9,18])

    def test_scalarMatrixMult(self):
        self.assertEqual(scalar_matrix_mult(matrix1,n),[[9,18,27],[36,45,54],[63,72,81]])
        self.assertEqual(scalar_matrix_mult(matrix3,n),[[36,45,54],[27,45,72],[18,36,9]])

    def test_addMatrix(self):
        self.assertEqual(add_matrix(matrix1,matrix2),[[2, 4, 6], [8, 10, 12], [14, 10, 18]])
        self.assertEqual(add_matrix(matrix1,matrix3),[[5,7,9],[7,10,14],[9,12,10]])

    def test_matrixVectorMult(self):
        self.assertEqual(matrix_vector_mult(matrix1,test_vector_01),[37,44,51])
        self.assertEqual(matrix_vector_mult(matrix1,test_vector_02),[21,27,33])

    def test_matrixMatrixMult(self):
        self.assertEqual(matrix_matrix_mult(matrix1,matrix2),[[30, 36, 42], [66, 81, 96], [78, 96, 114]])
        self.assertEqual(matrix_matrix_mult(matrix1,matrix3),[[66, 81, 96], [79, 95, 111], [25, 32, 39]])

    def test_absoluteValue(self):
        self.assertEqual(absoluteValue(-3),3)
        self.assertAlmostEqual(absoluteValue(complex(3,2)),13**(1/2))

    def test_pNorm(self):
        self.assertAlmostEqual(pNorm(test_vector_01,3),73**(1/3))
        self.assertAlmostEqual(pNorm(test_vector_01),21**(1/2))

    def test_infNorm(self):
        self.assertEqual(infNorm(test_vector_01),4)
        self.assertEqual(infNorm(test_vector_02),3)

    def test_pNormOrInf(self):
        self.assertEqual(pNormOrInf(test_vector_01,infinite=True),4)
        self.assertAlmostEqual(pNormOrInf(test_vector_01),21**(1/2))

    def test_innerProduct(self):
        self.assertEqual(innerProduct(test_vector_01,test_vector_02),13)
        self.assertEqual(innerProduct(test_vector_01,test_vector_03),complex(40,4))
