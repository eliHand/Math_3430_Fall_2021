from typing import List
import unittest

Vector = List[float]
Matrix = List[Vector]

def absoluteValue(scalar: float) -> float:
    """
    Calculates the absolute value of a scalar. Also works on complex numbers.
    Uses the fact that the sqrt(real^2 + imag^2) == magnitude a.k.a the absolute value
    of complex numbers. Also uses the alternative definition of absolute value for real
    values sqrt(number^2). The positive solution is exactly the absolute value.

    Args: 
    scalar: A floating point or complex number

    Returns:
    Floating point number representing the absolute value of the scalar, or
    the magnitude of a complex number.

    """ 
    if (type(scalar) == type(complex())):
        return ((scalar.real)**2 + (scalar.imag)**2)**(1/2)
    else:
        return (scalar**2)**(1/2)

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
        total += (absoluteValue(entry))**scalar
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
    return absoluteValue(max(v1,key=lambda x: absoluteValue(x)))


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






if __name__ == '__main__':
    unittest.main()
