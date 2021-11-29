import unittest
from QR import *
from LS import *
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

class TestLeastSquares(unittest.TestCase):
    def test_leastSquares(self):
        self.assertAlmostEqual(leastSquares(matrixNotes,test_vector_01),[3.333333333333333, 2.6666666666666665, 1.666666666666667])
        self.assertAlmostEqual(leastSquares(matrix2,test_vector_01),[4.543441112511214, -0.4364357804719846, -0.40824829046386335])

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


class TestQR(unittest.TestCase):
    
    def test_projection(self):
        self.assertAlmostEqual(projection(test_vector_01,test_vector_02),[13/21,26/21,52/21])
        self.assertAlmostEqual(projection(test_vector_01,test_vector_03),[(1.9047619047619047+0.19047619047619047j), (3.8095238095238093+0.38095238095238093j), (7.619047619047619+0.7619047619047619j)])

    def test_createEmptyMatrix(self):
        self.assertEqual(createEmptyMatrix(matrix1),[[0,0,0],[0,0,0],[0,0,0]])
        self.assertEqual(createEmptyMatrix(matrix2),[[0,0,0],[0,0,0],[0,0,0]])

    def test_transpose(self):
        self.assertEqual(transpose(matrix1),[[1,4,7],[2,5,8],[3,6,9]])
        self.assertEqual(transpose(matrix2),[[1, 4, 7], [2, 5, 2], [3, 6, 9]])

    def test_orthoNormalizeVectors(self):
        Q = orthoNormalizeVectors(matrix2)
        SolutionQ = [[14**(.5)/14,14**(.5)/7,3*14**(.5)/14],\
                [4*21**(.5)/21,21**(.5)/21,-2*21**(.5)/21],\
                [6**(.5)/6,-1*6**(.5)/3,6**(.5)/6]]
        for vector in range(len(Q)):
            for item in range(len(Q[vector])):
                self.assertAlmostEqual(Q[vector][item],SolutionQ[vector][item])
        Q = orthoNormalizeVectors(matrixNotes)
        SolutionQ = [[0.6666666666666666, 0.6666666666666666, 0.3333333333333333], [-0.6666666666666666, 0.3333333333333333, 0.6666666666666666], [-0.33333333333333337, 0.6666666666666667, -0.6666666666666666]]
        for vector in range(len(Q)):
            for item in range(len(Q[vector])):
                self.assertAlmostEqual(Q[vector][item],SolutionQ[vector][item])

    def test_GramSchmidtStable(self):
        QR = GramSchmidtStable(matrix2)
        SolutionQ = [[14**(.5)/14,14**(.5)/7,3*14**(.5)/14],\
                [4*21**(.5)/21,21**(.5)/21,-2*21**(.5)/21],\
                [6**(.5)/6,-1*6**(.5)/3,6**(.5)/6]]
        SolutionR = [[14**(.5),0,0],[(16*14**(.5))/7,(3*21**(.5))/7,0],[(19*14**(.5))/7,(4*21**(.5))/7,2*6**(.5)]]
        for vector in range(len(QR[0])):
            for item in range(len(QR[0][vector])):
                self.assertAlmostEqual(QR[0][vector][item],SolutionQ[vector][item])
        for vector in range(len(QR[1])):
            for item in range(len(QR[1][vector])):
                self.assertAlmostEqual(QR[1][vector][item],SolutionR[vector][item])
        QR = GramSchmidtStable(matrixNotes)
        SolutionQ = [[0.6666666666666666, 0.6666666666666666, 0.3333333333333333], [-0.6666666666666666, 0.3333333333333333, 0.6666666666666666], [-0.33333333333333337, 0.6666666666666667, -0.6666666666666666]]
        SolutionR = [[3.0, 0.0, 1.1102230246251565e-16], [0.0, 3.0, 2.220446049250313e-16], [3.0, 1.0, 0.9999999999999999]]
        for vector in range(len(QR[0])):
            for item in range(len(QR[0][vector])):
                self.assertAlmostEqual(QR[0][vector][item],SolutionQ[vector][item])
        for vector in range(len(QR[1])):
            for item in range(len(QR[1][vector])):
                self.assertAlmostEqual(QR[1][vector][item],SolutionR[vector][item])
        

    def test_copyMatrix(self):
        self.assertEqual(matrix1,copyMatrix(matrix1))
        self.assertEqual(matrix2,copyMatrix(matrix2))

    def test_makeColumn(self):
        self.assertEqual(makeColumn(test_vector_01,0),[1,0,0])
        self.assertEqual(makeColumn(test_vector_02,0),[1,0,0])

    def test_reconstructMatrix(self):
        self.assertEqual(reconstructMatrix(matrixTest,matrix1),[[1,0,0],[0,1,2],[0,3,4]])
        self.assertEqual(reconstructMatrix(matrixTest2,matrix1),[[1, 0, 0], [0, 1, 0], [0, 0, 4]])

    def test_makeIdentMatrix(self):
        self.assertEqual(makeIdentMatrix(matrix1),[[1,0,0],[0,1,0],[0,0,1]])
        self.assertEqual(makeIdentMatrix(matrix2),[[1,0,0],[0,1,0],[0,0,1]])

    def test_makeSubMatrix(self):
        self.assertEqual(makeSubMatrix(matrix1),[[5,6],[8,9]])
        self.assertEqual(makeSubMatrix(matrix2),[[5, 6], [2, 9]])

    def test_houseHolderCalc(self):
        self.assertAlmostEqual(houseHolderCalc(matrixNotes,0),[[0.6666666666666667, 0.6666666666666666, 0.3333333333333333],\
             [0.6666666666666666, -0.33333333333333326, -0.6666666666666666],\
             [0.3333333333333333, -0.6666666666666666, 0.6666666666666667]])
        self.assertAlmostEqual(houseHolderCalc(matrix2,0),[[0.2672612419124244, 0.5345224838248488, 0.8017837257372732], [0.5345224838248488, 0.6100734640269463, -0.5848898039595805], [0.8017837257372732, -0.5848898039595805, 0.12266529406062932]])

    def test_houseHolderDriver(self):
        QR = houseHolderDriver(matrixNotes)
        Q = QR[1]
        R = QR[0]
        self.assertAlmostEqual(R,[[3.0000000000000004, -1.1102230246251564e-16, 1.1102230246251565e-16],\
             [-2.220446049250313e-16, 3.0, 0.0], [3.0000000000000004, 0.9999999999999997, -1.0]])
        self.assertAlmostEqual(Q,[[0.6666666666666667, 0.6666666666666666, 0.3333333333333333], [-0.6666666666666666, 0.3333333333333332, 0.6666666666666666], [0.33333333333333337, -0.6666666666666666, 0.6666666666666667]])
        QR = houseHolderDriver(matrix2)
        Q = QR[1]
        R = QR[0]
        self.assertAlmostEqual(R,[[3.7416573867739413, 1.1517126582049866e-16, -1.8984042273872353e-16], [8.55235974119758, 1.963961012123932, -2.220446049250313e-16], [10.155927192672127, 2.618614682831909, -4.898979485566354]])
        self.assertAlmostEqual(Q,[[0.2672612419124244, 0.5345224838248488, 0.8017837257372732], [0.8728715609439697, 0.2182178902359923, -0.4364357804719847], [-0.4082482904638628, 0.8164965809277258, -0.408248290463863]])

if __name__ == "__main__":
    unittest.main()