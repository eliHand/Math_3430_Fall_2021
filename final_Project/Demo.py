from LA import *
from QR import *
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



if __name__ == "__main__":
    print("My name is Elijah Hand, and this is the Linear Algebra library that I have written for Math 3430.")
    print("The following will be a demonstration of the various functions.\n\
    First off we have the LA.py file. This file performs various Linear Algebra Operations.")
    

    print("This is the add_vectors function. It takes in two vectors and returns their sum.")
    print("Take for example the two vectors [1,2,4] and [3,1,2]. The function will return:")
    print(add_vectors(test_vector_01,test_vector_02))
    print("This is the scalar_vector function. It takes in a vector and scalar and returns vec*scalar")
    print("Take for example the vector [1,2,4] and scalar 9. The function will return:")
    print(scalar_vector(test_vector_01,n))
    print("This is the scalar_matrix_mult function. It takes in a matrix and scalar and returns mat*scalar.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,8,9]] and scalar 9. The function will return:")
    print(scalar_matrix_mult(matrix1,n))
    print("This is the add_matrix function. It takes in two matricies and returns their sum.")
    print("Take for example the two matricies [[1,2,3],[4,5,6],[7,8,9]] and [[1,2,3],[4,5,6],[7,2,9]]. The function will return:")
    print(add_matrix(matrix1,matrix2))
    print("This is the matrix_vector_mult function. It takes in a matrix and vector and returns their product.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,8,9]] and vector [1,2,4] The function will return:")
    print(matrix_vector_mult(matrix1,test_vector_01))
    print("This is the matrix_matrix_mult function. It takes in two matricies and returns their product.")
    print("Take for example the two matricies [[1,2,3],[4,5,6],[7,8,9]] and [[1,2,3],[4,5,6],[7,2,9]]. The function will return:")
    print(matrix_matrix_mult(matrix1,matrix2))
    print("This is the pNorm function. It takes in a vector and scalar returns its scalar-norm (defaults to 2).")
    print("Take for example the vector [1,2,4] and scalar 2. The function will return:")
    print(pNorm(test_vector_01))
    print("This is the infNorm function. It takes in a vector and returns its infinite Norm.")
    print("Take for example the two vector [1,2,4]. The function will return:")
    print(infNorm(test_vector_01))
    print("This is the pNormOrInf function. It takes in a vector, scalar, and boolean and returns either its infinite or n-th norm.")
    print("Take for example the vector [1,2,4], scalar 2 and boolean False. The function will return:")
    print(pNormOrInf(test_vector_01))
    print("This is the innerProduct function. It takes in two vectors and returns their inner product.")
    print("Take for example the two vectors [1,2,4] and [3,1,2]. The function will return:")
    print(innerProduct(test_vector_01,test_vector_02))
    


    print("This is the projection function. It takes in two vectors and returns their projection.")
    print("Take for example the two vectors [1,2,4] and [3,1,2]. The function will return:")
    print(projection(test_vector_01,test_vector_02))
    print("This is the createEmptyMatrix function. It takes in a matrix and returns an empty matrix the same size.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,8,9]]. The function will return:")
    print(createEmptyMatrix(matrix1))
    print("This is the transpose function. It takes in a matrix and returns the transpose of the matrix.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,8,9]]. The function will return:")
    print(transpose(matrix1))
    print("This is the orthoNormalizeVectors function. It takes in a matrix and  Returns a Matrix containing vectors on the same span as the matrix m but orthonormalized.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,2,9]]. The function will return:")
    print(orthoNormalizeVectors(matrix2))
    print("This is the GramSchmidtStable function. It takes in a matrix and returns a list with the Q-R Decomposition. Q in the first entry and R in the second entry.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,2,9]]. The function will return:")
    print(GramSchmidtStable(matrix2))
    print("This is the copy function. It takes in a matrix and returns an identical matrix with the same entries.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,8,9]]. The function will return:")
    print(copyMatrix(matrix1))
    print("This is the makeColumn function. It takes in a vector and an integer returns a vector the same size with 0's and a 1 at the integer position.")
    print("Take for example the vector [1,2,4] and 0. The function will return:")
    print(makeColumn(test_vector_01,0))
    print("This is the reconstructMatrix function. It takes in two matricies and returns a reconstructed matrix the size of the second argument, with an identity shell around the first argument.")
    print("Take for example the two matricies [[1,2],[3,4]] and [[1,2,3],[4,5,6],[7,8,9]]. The function will return:")
    print(reconstructMatrix(matrixTest,matrix1))
    print("This is the makeIdentMatrix function. It takes in a matrix and returns an identity matrix the same size.")
    print("Take for example the matrix [[1,2,3],[4,5,6],[7,8,9]]. The function will return:")
    print(makeIdentMatrix(matrix1))
    print("This is the makeSubMatrix function. It takes in a matrix and returns the sub-matrix, which is taking off the outer entries.")
    print("Take for example the two vectors [[1,2,3],[4,5,6],[7,8,9]]. The function will return:")
    print(makeSubMatrix(matrix1))
    print("This is the houseHolderCalc function. It takes in a matrix and an integer and returns the householder transformation matrix.")
    print("Take for example the matrix [[2,2,1],[-2,1,2],[1,3,1]]. The function will return:")
    print(houseHolderCalc(matrixNotes,0))
    print("This is the houseHolderDriver function. It takes in a matrix and returns the QR decomposition using HouseHolder transformations.")
    print("Take for example the matrix [[2,2,1],[-2,1,2],[1,3,1]]. The function will return:")
    print(houseHolderDriver(matrixNotes))


    print("This is the leastSquares function. It takes in a matrix and vector and returns the least squares solution.")
    print("Take for example the matrix [[2,2,1],[-2,1,2],[1,3,1]] and vector [1,2,4]. The function will return:")
    print(leastSquares(matrixNotes,test_vector_01))
