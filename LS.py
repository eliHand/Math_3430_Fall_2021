from QR import houseHolderDriver, almostProjection
from QR import Matrix,Vector

def leastSquares(inputMatrix:Matrix,inputVector:Vector)->Vector:
    '''
    Calculates the least squares solution. Uses the QR factorization of the
    householder function.

    Args:
        inputMatrix: A matrix stored as a list of lists which you want to find the least squares solution for
        inputVector: A vector stored as a list which corresponds with the inputMatrix which you want to find the least
        Squares solution for.
    
    Returns:
        Will return a Vector that is the least squares solution.
    '''
    leastSquaresVec = []
    # First Step we need to do is obtain the orthogonal matrix for the input Matrix
    orthog = houseHolderDriver(inputMatrix)[1]
    for vector in orthog:
        leastSquaresVec.append(almostProjection(vector,inputVector))
    return leastSquaresVec

