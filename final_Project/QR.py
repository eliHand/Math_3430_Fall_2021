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

def almostProjection(vectorU:Vector,vectorV:Vector)->float:
    '''
    Does the "projection formula" ((v.u1)/(u1.u1))*u1 without the final
    multiplication of Un.

    Args:
    vectorU: The spanning vector, stored as a List
    vectorV: The vector that is being projected. Stored as a List.
    
    Returns:
    Returns a float that is the result of the inner products ratioed together.
    Essentially the projection function without the final scalar_vector call.    
    '''
    numerator = innerProduct(vectorU,vectorV)
    denominator = innerProduct(vectorU,vectorU)
    quotient = numerator/denominator
    return quotient

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
            temp = scalar_vector(temp,-1)
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
    '''
    This function takes in a matrix m, and iterates entry
    by entry, copying it without reference to a new matrix
    which it then returns.

    Args:
        m: A matrix of vectors stored as list of lists

    Returns:
        Returns a matrix with the same entries as m without
        reference to them.
    '''
    temp = createEmptyMatrix(m)
    for i in range(len(m)):
        for j in range(len(m[i])):
            temp[i][j] = m[i][j]
    return temp



def makeColumn(v:Vector,n:int)->Vector:
    '''
    Creates a column vector of size v filled with 0's with n denoting
    the location of a 1. Then returns this vector.

    Args:
        v: A vector v stored as list
        n: An integer denoting where a 1 will be inserted into the new vector
    
    Returns
        Returns a vector composed of a 1 at position n and 0s everywhere else
    '''
    temp = []
    for i in v:
        temp.append(0)
    temp[n] = 1
    return temp

def reconstructMatrix(matrixToReconstruct:Matrix,BaseMatrix:Matrix)->Matrix:
    '''
    Reconstructs matrixToReconstruct to the same size as BaseMatrix
    with an "identity shell" as outlined below.

    1 2     1 0 0    1 0 0 0
    3 4 ->  0 1 2 -> 0 1 0 0
            0 3 4    0 0 1 2
                     0 0 3 4

    Args:
        matrixToReconstruct: A matrix stored as a list of lists
        BaseMatrix: A matrix stored as a list of lists

    Returns:
        Returns a matrix with an inner part as matrixToReconstruct and 
        the size BaseMatrix with an "identity shell" to make up the size
        difference.
    '''
    copyMatrixToReconstruct = copyMatrix(matrixToReconstruct)
    while len(BaseMatrix) > len(copyMatrixToReconstruct):
        for i in copyMatrixToReconstruct:
            i.insert(0,0)
        tempCol = makeColumn(copyMatrixToReconstruct[0],0)
        copyMatrixToReconstruct.insert(0,tempCol)
    return copyMatrixToReconstruct
    

def makeIdentMatrix(m:Matrix)->Matrix:
    '''
    Creates an identity matrix the same size as m.

    Args:
        m: A matrix stored as a list of lists

    Returns:
        Returns an identity matrix of the same size a m.
    '''
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
    Creates a submatrix out of m that slices off the first row
    and column of m. As outlined as below.

    1 4 7    5 8
    2 5 8 -> 6 9 -> 9
    3 6 9

    Args:
        m: A matrix stored as a list of lists.

    Returns:
        Returns the sliced matrix.
    '''
    temp = copyMatrix(m)
    temp.pop(0)
    for i in temp:
        i.pop(0)
    return temp

def houseHolderCalc(m:Matrix, n:int)->Matrix:
    '''
    Calculates the householder transformation matrix for m, Matrix.
    Is created by the formula outlined in class.

    Args:
        m: A matrix stored as a list of lists.
        n: An integer n (Should be set to 0 always)

    Returns:
        Returns the householder transformation matrix for m.
    '''
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

    


def houseHolderDriver(m:Matrix)->List:
    '''
    Drives the houseHolderCalc function. Basically slices and constructs the matricies
    to form each transformation matrix and multiplies it with m for the next iteration.
    Continues until the slices make a 2x2 matrix (the smallest that is useful)
    Args:
        m: A matrix, a list of lists
    
    Returns:
        Returns a list with the first entry being the transformed matrix from m. Will be upper-triangular.
        The second entry will be the Q matrix, containing the orthogonal vectors that combine with R to form 
        the input matrix.
    '''
    index = 2
    tempMatrix = []
    Q = houseHolderCalc(m,0)
    firstIteration = matrix_matrix_mult(Q,m)
    returnList = [0,0]
    while index < len(firstIteration):
        tempMatrix = makeSubMatrix(firstIteration)
        tempMatrix = houseHolderCalc(tempMatrix,0)
        reconstructed = reconstructMatrix(tempMatrix,m)
        Q = matrix_matrix_mult(Q,reconstructed)
        firstIteration = matrix_matrix_mult(reconstructed,firstIteration)
        index += 1
    returnList[0] = firstIteration
    returnList[1] = Q
    return returnList






