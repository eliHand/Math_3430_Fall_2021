from typing import List
import unittest

Vector = List[float]
Matrix = List[Vector]

def absoluteValue(scalar: float) -> float:
    """
    Calculates the absolute value of a scalar. Also works on complex numbers.
    Uses Python's builtin abs() function which gives the absolute value for normal
    floats, and returns the magnitude in the case of complex numbers.

    Args: 
    scalar: A floating point or complex number

    Returns:
    Floating point number representing the absolute value of the scalar, or
    the magnitude of a complex number.

    """ 
    return abs(scalar)

def pNorm(v1 : Vector, scalar:float=2) -> float:
    """
    Calculates the pNorm of a vector. Defaults to p = 2 if scalar value is not provided.
    Uses the formula: pNorm = (Sum(abs(entryOfVector)^scalar))^1/scalar.

    Args:
    v1: A Vector stored as a list.
    scalar: A floating point number. Defaults to 2 if none is provided.

    Returns:
    Floating point number representing the pNorm of v1.
    
    """
    total = 0
    for entry in v1:
        total += (abs(entry))**scalar
    return total**(1/scalar)


def infNorm(v1:Vector)->float:
    """
    Calculates the infinite Norm of v1. If we consider each entry of a vector to be
    its absolute value, then this turns out to be the maximum of those entries.

    Args:
    v1: Vector stored as a list.

    Returns:
    A float representing the infinite Norm of v1.
    """
    return abs(max(v1,key=lambda x: abs(x)))


def pNormOrInf(v1: Vector,scalar:float=2,infinite:bool=False)->float:
    """
    Chooses between the pNorm of Infinite Norm using a boolean value.
    Defaults to False -> returns the pNorm.

    Args:
    v1: Vector stored as a list.
    scalar: Float defaults to 2 if none is supplied.
    infinite: Boolean to decide to calculate infinite norm. Defaults to False.

    Returns:
    Returns a floating point number representing the infinite norm if infinite==True.
    Else returns a floating point number representing the pNorm of v1 with p=scalar.
    
    """
    if infinite:
        return infNorm(v1)
    else:
        return pNorm(v1,scalar)

def innerProduct(v1:Vector,v2:Vector)->float:
    """
    Calculates the inner product of v1 and v2. Simply sums together the multiplication
    of each corresponding entry of v1 and v2.

    Args:
    v1: A vector stored as a List.
    v2: A vector stored as a List.
    Both v1 and v2 must have the same size.

    Returns:
    Returns a float which represents the inner product (dot product) of
    v1 and v2.
    
    """
    sum = 0
    for i in range(len(v1)):
        sum += v1[i]*v2[i]
    return sum


def add_vectors(vector_a: Vector,vector_b: Vector) -> Vector:
    """Adds the two input vectors.

    Creates a result vector stored as a list of 0's the same length as the input 
    then overwrites each element of the result vector with the corresponding
    element of the sum of the input vectors. Achieves this using a for loop over
    the indices of result. 

    Args:
        vector_a: A vector stored as a list.
        vector_b: A vector, the same length as vector_a, stored as a list.

    Returns:
       The sum of the input vectors stored as a list. 
    """ 
    result = [0 for element in vector_a]
    for index in range(len(result)):
        result[index] = vector_a[index] + vector_b[index]
    return result


def scalar_vector(vector_a : Vector,n: int) -> Vector:
    """
    Performs scalar vector multiplication.

    Creates a list which is populated by a for loop looping over vector a and appending the resulting
    multiplication between vector_a[i] and n.

    Args:
        vector_a: A vector stored as a list.
        n: A scalar. Acceptable to be either a float or an integer.
    
    Returns:
        Returns the resulting vector stored as a list.
    """
    returnList = []
    for i in vector_a:
        returnList.append(i*n)
    return returnList

def scalar_matrix_mult(matrix_a: Matrix,n: int) -> Matrix:
    """
    Performs scalar matrix multiplication.

    Creates a return matrix as a list, whose elements are vectors (also stored as lists)
    These are populated by performing scalar vector multiplication on each column of matrix_a.

    Args:
        matrix_a: A matrix stored as a list of lists where each is a column vector.
        n: A scalar. Acceptable to be either a float or an integer.

    Returns:
        The return matrix, a list of lists. 
    """
    returnMatrix = []
    for i in matrix_a:
        returnMatrix.append(scalar_vector(i,n))
    return returnMatrix

def add_matrix(matrix_a: Matrix,matrix_b: Matrix) -> Matrix:
    """
    Adds the two input matricies.

    Creates a return matrix which is the populated by the columns of matrix_a
    and matrix_b added together using the add_vectors function.

    Args:
        matrix_a: A matrix, stored as a list of lists.
        matrix_b: A matrix, stored as a list of lists.

    Returns:
        A matrix stored as a list of lists of matrix_a and matrix_b added together.
    """
    returnMatrix = []
    for index in range(len(matrix_a)):
        tempRow = add_vectors(matrix_a[index],matrix_b[index])
        returnMatrix.append(tempRow)
    return returnMatrix

def matrix_vector_mult(matrix_a: Matrix,vector_a: Vector) -> Vector:
    """
    Performs matrix vector multiplication.

    Creates a temporary list. This is used to store the scalar vector multiplications
    of each column in matrix_a and each entry in vector_a. Then each of these entries are added
    together. This sum is then returned.

    Args: 
        matrix_a: A matrix stored as a list of lists.
        vector_a: A vector stored as a list.
    
    Returns:
        The resulting matrix vector multiplication as a vector stored as a list.
    """
    tempList = []
    for i in range(len(vector_a)):
        tempEntry = vector_a[i]
        tempCol = matrix_a[i]
        colAfter = scalar_vector(tempCol,tempEntry)
        tempList.append(colAfter)
    for i in range(len(tempList)):
        if ((i+1) == len(tempList)):
            break
        tempList[i+1] = add_vectors(tempList[i],tempList[i+1])
    return tempList[-1]

def matrix_matrix_mult(matrix_a: Matrix,matrix_b: Matrix) -> Matrix:
    """
    Performs matrix matrix multiplication.

    Creates a temporary list which is populated by the matrix vector multiplication
    of matrix_a and each column vector of matrix_b.

    Args:
        matrix_a: A matrix stored as a list of lists.
        matrix_b: A matrix stored as a list of lists.

    Returns:
        A matrix containing the matrix matrix multiplication.
    """
    tempList = []
    for i in range(len(matrix_b)):
        tempCol = matrix_b[i]
        solVector = matrix_vector_mult(matrix_a,tempCol)
        tempList.append(solVector)
    return tempList

#Test Inputs
matrix1 = [[1,2,3],[4,5,6],[7,8,9]]
matrix2 = [[6,8,4],[3,5,1],[7,2,9]]
matrix3 = [[4,5,6],[3,5,8],[2,4,1]]
test_vector_01 = [1, 2, 4]
test_vector_02 = [3, 1, 2]
test_vector_03 = [6,complex(3,2),7]
n = 9

class TestVectorOperations(unittest.TestCase):

    def test_addVectors(self):
        self.assertEqual(add_vectors(test_vector_01,test_vector_02),[4,3,6])
    
    def test_scalarVecotr(self):
        self.assertEqual(scalar_vector(test_vector_01,n),[9,18,36])
        self.assertEqual(scalar_vector(test_vector_02,n),[27,9,18])

    def test_scalarMatrixMult(self):
        self.assertEqual(scalar_matrix_mult(matrix1,n),[[9,18,27],[36,45,54],[63,72,81]])
        self.assertEqual(scalar_matrix_mult(matrix3,n),[[36,45,54],[27,45,72],[18,36,9]])

    def test_addMatrix(self):
        self.assertEqual(add_matrix(matrix1,matrix2),[[7,10,7],[7,10,7],[14,10,18]])
        self.assertEqual(add_matrix(matrix1,matrix3),[[5,7,9],[7,10,14],[9,12,10]])

    def test_matrixVectorMult(self):
        self.assertEqual(matrix_vector_mult(matrix1,test_vector_01),[37,44,51])
        self.assertEqual(matrix_vector_mult(matrix1,test_vector_02),[21,27,33])

    def test_matrixMatrixMult(self):
        self.assertEqual(matrix_matrix_mult(matrix1,matrix2),[[66, 84, 102], [30, 39, 48], [78, 96, 114]])
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



if __name__ == '__main__':
    unittest.main()
