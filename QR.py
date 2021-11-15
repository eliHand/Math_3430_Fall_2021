from LA import *
from typing import List
import unittest




def projection(vectorU:Vector,vectorV:Vector)->Vector:
    '''
    Projects the vector V orthonally onto the line spanned by
    vectorU.

    Args:
    vectorU: The spanning vector, stored as a List
    vectorV: The vector that is being projected. Stored as a List.
    
    Returns:
    Returns a vector stored as a list that is the result of the projection of V onto U.    
    '''
    numerator = innerProduct(vectorU,vectorV)
    denominator = innerProduct(vectorU,vectorU)
    quotient = numerator/denominator
    return scalar_vector(vectorU,quotient)

def addInverse(vectorA:Vector):
    '''
    Converts vectorA in place with each entry equal to
    its additive inverse.

    Args:
    vectorA: A vector stored as a list. Will be mutated so each entry
    is the additive inverse.

    Returns:
    None
    '''
    for i in range(len(vectorA)):
        vectorA[i] = -vectorA[i]

def createEmptyMatrix(m:Matrix)->Matrix:
    '''
    Creates a matrix the same size as m, but filled with 0s. Does not harm the input matrix.
    Created as a helper function for transpose.

    Args:
    m: A matrix stored as a list of lists.

    Returns:
    Matrix the same size as m filled with 0.
    '''
    returnMatrix = []
    for index in range(len(m)):
        returnMatrix.append(list())
        for j in range(len(m[index])):
            returnMatrix[index].append(0)
    return returnMatrix

def transpose(m:Matrix)->Matrix:
    '''
    Creates a new matrix that is simply the matrix m transposed.

    Args:
    m: A matrix stored as a list of lists

    Returns:
    Matrix that is the matrix m transposed.
    '''
    returnMatrix = createEmptyMatrix(m)
    for i in range(len(m)):
        for j in range(len(m)):
            returnMatrix[i][j] = m[j][i]
    return returnMatrix



def GramSchmidtStable(m:Matrix)->List:
    '''
    Implements the stable version of the Gram-Schmidt Algorithm for calculating the
    Q-R Decomposition of a matrix. 

    Input: A matrix m which is a list of lists.

    Output: Outputs a list containing two items. First Entry is Q in the Q-R Decomposition.
    The second is R, an upper triangular matrix from the Decomposition.
    '''
    V = []
    Q = []
    for i in m:
        V.append(i)
    for j in range(len(V)):
        tempLen = pNorm(V[j])
        tempVect = []
        for entry in range(len(V)):
            tempVect.append(V[j][entry]/tempLen)
        Q.append(tempVect)
        for k in range(j+1,len(V)):
            temp = projection(Q[j],V[k])
            addInverse(temp)
            V[k] = add_vectors(V[k],temp)
    # At this point, Q is now a set of orthonormal vectors.
    transposeOrthog = transpose(Q)
    R = matrix_matrix_mult(transposeOrthog,m)
    returnList = [Q,R]
    return returnList

def orthoNormalizeVectors(m:Matrix)->Matrix:
    '''
    This function will be supplied a list of vectors (a matrix), then will perform
    the Gram-Schmidt process, and will return Q from that process, which will be an
    orthonormal list of vectors which share the same span as those from the input matrix.

    Input:
        m: A matrix of vectors stored as a list of lists.

    Returns:
        Returns a Matrix containing vectors on the same span as the matrix m but
        orthonormalized.
    '''
    holder = GramSchmidtStable(m)
    return holder[0]

def copyMatrix(m:Matrix)->Matrix:
    temp = createEmptyMatrix(m)
    for i in range(len(m)):
        for j in range(len(m[i])):
            temp[i][j] = m[i][j]
    return temp



def makeColumn(v:Vector,n:int)->Vector:
    temp = []
    for i in v:
        temp.append(0)
    temp[n] = 1
    return temp

def reconstructMatrix(matrixToReconstruct:Matrix,BaseMatrix:Matrix)->Matrix:
    '''
    1 2     1 0 0    1 0 0 0
    3 4 ->  0 1 2 -> 0 1 0 0
            0 3 4    0 0 1 2
                     0 0 3 4
    '''
    copyMatrixToReconstruct = copyMatrix(matrixToReconstruct)
    while len(BaseMatrix) > len(copyMatrixToReconstruct):
        for i in copyMatrixToReconstruct:
            i.insert(0,0)
        tempCol = makeColumn(copyMatrixToReconstruct[0],0)
        copyMatrixToReconstruct.insert(0,tempCol)
    return copyMatrixToReconstruct
    

def makeIdentMatrix(m:Matrix)->Matrix:
    tempMatrix = createEmptyMatrix(m)
    for i in range(len(m)):
        for j in range(len(m[i])):
            if (i==j):
                tempMatrix[i][j] = 1
            else:
                tempMatrix[i][j] = 0
    return tempMatrix

def makeSubMatrix(m:Matrix)->Matrix:
    '''
    1 4 7    5 8
    2 5 8 -> 6 9 -> 9
    3 6 9
    '''
    temp = copyMatrix(m)
    temp.pop(0)
    for i in temp:
        i.pop(0)
    return temp

def houseHolderCalc(m:Matrix, n:int)->Matrix:
    x = m[n] # First column of m
    normX = pNorm(x)
    inverseX = scalar_vector(x,-1)
    e = makeColumn(x,n)
    v = add_vectors((scalar_vector(e,normX)),inverseX)
    # Now we calculate F
    numerator = []
    for entry in v:
        column = scalar_vector(v,entry)
        numerator.append(column)
    deno = 2/innerProduct(v,v)
    temp = scalar_matrix_mult(numerator,deno)
    ident = makeIdentMatrix(m)
    temp = scalar_matrix_mult(temp,-1)
    F = add_matrix(temp,ident)
    return F

def houseHolderDriver(m:Matrix)->Matrix:
    index = 2
    tempMatrix = []
    firstIteration = matrix_matrix_mult(houseHolderCalc(m,0),m)
    while index < len(firstIteration):
        tempMatrix = makeSubMatrix(firstIteration)
        tempMatrix = houseHolderCalc(tempMatrix,0)
        reconstructed = reconstructMatrix(tempMatrix,m)
        firstIteration = matrix_matrix_mult(reconstructed,firstIteration)
        index += 1 
    return firstIteration


matrixNotes = [[2,2,1],[-2,1,2],[1,3,1]]
matrixTest = [[1,2],[3,4]]




matrix1 = [[1,2,3],[4,5,6],[7,8,9]]
matrix2OLD = [[6,8,4],[3,5,1],[7,2,9]]
matrix2 = [[1,2,3],[4,5,6],[7,2,9]]
matrix3 = [[4,5,6],[3,5,8],[2,4,1]]
matrixCJ = [[1,2,3],[4,5,6]]
test_vector_01 = [1, 2, 4]
test_vector_02 = [3, 1, 2]
test_vector_03 = [6,complex(3,2),7]
n = 9

class TestQR(unittest.TestCase):
    
    def test_projection(self):
        self.assertAlmostEqual(projection(test_vector_01,test_vector_02),[13/21,26/21,52/21])

    def test_addInverse(self):
        addInverse(test_vector_01)
        self.assertEqual(test_vector_01,[-1,-2,-4])

    def test_createEmptyMatrix(self):
        self.assertEqual(createEmptyMatrix(matrix1),[[0,0,0],[0,0,0],[0,0,0]])

    def test_transpose(self):
        self.assertEqual(transpose(matrix1),[[1,4,7],[2,5,8],[3,6,9]])

    def test_orthoNormalizeVectors(self):
        Q = orthoNormalizeVectors(matrix2)
        SolutionQ = [[14**(.5)/14,14**(.5)/7,3*14**(.5)/14],\
                [4*21**(.5)/21,21**(.5)/21,-2*21**(.5)/21],\
                [6**(.5)/6,-1*6**(.5)/3,6**(.5)/6]]
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



if __name__ == '__main__':
    unittest.main()
    
