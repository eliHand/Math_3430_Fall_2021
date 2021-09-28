"""
For this homework assignment we will take our work from HW01 and use it to
prepare a python script which will implement our algoirthms as python functions. 

For Problems #0-5 from HW01, Do the following.



1) Write your answer from HW01 in a comment.

2) Below the comment write a function which implements the algorithm from your
comment. If you find that you need to change your algorithm for your python
code, you must edit your answer in the comment. 

3) Test each of your functions on at least 2 inputs. 

4) Upload your .py file to a github repo named "Math_3430_Fall_2021"

This assignment is due by 11:59pm 09/27/2021. Do NOT upload an updated version to github
after that date. 
"""


#Example:

#Problem 00

"""
-The Three Questions

Q1: What do we have?

A1: Two vectors stored as lists. Denoted by the names vector_a and vector_b. 

Q2: What do we want?

A2: Their sum stored as a list.

Q3: How will we get there?

A3: We will create an empty list of the appropriate size and store the sums of
the corresponding components of vector_a and vector_b. 

-PsuedoCode

def add_vectors(vector_a,vector_b):

Initialize a result vector of 0's which is the same size as vector_a. Call this
vector result.

# Set each element of result to be equal to the desired sum.
for index in range(length(result)):
  result[index] = vector_a[index] + vector_b[index]

Return the desired result.
"""

def add_vectors(vector_a,vector_b):
  result = [0 for element in vector_a]
  for index in range(len(result)):
    result[index] = vector_a[index] + vector_b[index]
  return result

#End Example



#Problem 01
""""
Write an algorithm to implement scalar-vector multiplication.

Q1: What do we have?
A1: A vector (vector_a) stored as a list and a floating point/integer number (n)

Q2: What do we want?
A2: The scalar-vector multiplication stored as a list.

Q3: How do we get there?
A3: Create a new empty list. Take a scalar number n and then multiply it by every entry in vector_a,
appending each value in turn to the new list (which is the vector after the fact). Return this list.

def scalar_vector(vector_a,n):
	returnList = []
	for i in vector_a:
		returnList.append(i*n)
	return returnList
"""
def scalar_vector(vector_a,n):
	returnList = []
	for i in vector_a:
		returnList.append(i*n)
	return returnList

#Problem 02
"""
Write an algorithm to implement scalar-matrix multiplication.

Q1: What do we have?
A1: A matrix (matrix_a) stored as a 2-dimensional list and a scalar k.

Q2: What do we want?
A2: The scalar-matrix multiplication stored as a matrix

Q3: How do we get there?
A3: Take the first row of matrix_a as a vector. Apply Problem #1. Store this
new list as the first row of the new 2-dimensional list. Continue n times. 


def scalar_matrix_mult(matrix_a,n):
    returnMatrix = []
    for i in matrix_a:
        returnMatrix.append(scalar_vector(i,n))
    return returnMatrix
"""
def scalar_matrix_mult(matrix_a,n):
    returnMatrix = []
    for i in matrix_a:
        returnMatrix.append(scalar_vector(i,n))
    return returnMatrix

#Problem 03
"""
Write an algorithm to implement matrix addition.

Q1: What do we have?
A1: Two matricies (matrix_a, matrix_b) of the same size.

Q2: What do we want?
A2: A new matrix that is matrix_a + matrix_b elementwise.

Q3: How do we get there?
A3: Take row 1 from matrix_a and row 1 from matrix_b. Treat these as vectors and apply
Problem #0. This will be row 1 of the new matrix. Continue n times.

def add_matrix(matrix_a,matrix_b):
    returnMatrix = []
    for index in range(len(matrix_a)):
        tempRow = add_vectors(matrix_a[index],matrix_b[index])
        returnMatrix.append(tempRow)
    return returnMatrix
"""
def add_matrix(matrix_a,matrix_b):
    returnMatrix = []
    for index in range(len(matrix_a)):
        tempRow = add_vectors(matrix_a[index],matrix_b[index])
        returnMatrix.append(tempRow)
    return returnMatrix

#Problem 04
"""
Write an algorithm to implement matrix-vector multiplication using the linear
combination of columns method. You must use the algorithms from Problem #0 and
Problem #1.  

Q1: What do we have?
A1: A matrix (matrix_a) and a vector (vector_a).

Q2: What do we want?
A2: A new matrix that is the result of the matrix-vector multiplication

Q3: How do we get there?
A3: Take the first column of matrix_a, and the first entry of vector_a.
Then apply Problem #1. Store this temporarily in a list. Do the same
thing with the next column and entry. Continue n times.
Then apply Problem #0 to each pair of entries in the temporary list.
That sum will be the answer. Return that.

def column(matrix,col):
    return [row[col] for row in matrix]

def matrix_vector_mult(matrix_a,vector_a):
    tempList = []
    returnVector = vector_a
    for i in returnVector:
        returnVector = 0
    for i in range(len(vector_a)):
        tempCol = column(matrix_a,i)
        tempEntry = vector_a[i]
        colAfter = scalar_vector(tempCol,tempEntry)
        tempList.append(colAfter)
    for i in range(len(tempList)):
        if ((i+1) == len(tempList)):
            break
        tempList[i+1] = add_vectors(tempList[i],tempList[i+1])
    return tempList[-1]
"""
def column(matrix,col):
    return [row[col] for row in matrix]

def matrix_vector_mult(matrix_a,vector_a):
    tempList = []
    returnVector = vector_a
    for i in returnVector:
        returnVector = 0
    for i in range(len(vector_a)):
        tempCol = column(matrix_a,i)
        tempEntry = vector_a[i]
        colAfter = scalar_vector(tempCol,tempEntry)
        tempList.append(colAfter)
    for i in range(len(tempList)):
        if ((i+1) == len(tempList)):
            break
        tempList[i+1] = add_vectors(tempList[i],tempList[i+1])
    return tempList[-1]

#Problem 05
"""
Write an algorithm to implement matrix-matrix multiplication using your
algorithm from Problem #4.  

Q1: What do we have?
A1: Two matricies (matrix_a,matrix_b)

Q2: What do we want?
A2: A matrix that is the solution to the matrix multiplication of matrix_a and matrix_b.

Q3: How do we get there?
A3: Take the first column of matrix_b, this will serve as a vector. Apply Problem #4.
Store this value. This is the first column of the solution. Continue with the next column n times.
Then return the matrix with each entry being the column vector.

def matrix_matrix_mult(matrix_a,matrix_b):
    tempList = []
    for i in range(len(matrix_b)):
        tempCol = column(matrix_b,i)
        solVector = matrix_vector_mult(matrix_a,tempCol)
        tempList.append(solVector)
    return tempList
"""
def matrix_matrix_mult(matrix_a,matrix_b):
    tempList = []
    for i in range(len(matrix_b)):
        tempCol = column(matrix_b,i)
        solVector = matrix_vector_mult(matrix_a,tempCol)
        tempList.append(solVector)
    return tempList

#Test Inputs
matrix1 = [[1,2,3],[4,5,6],[7,8,9]]
matrix2 = [[6,8,4],[3,5,1],[7,2,9]]
matrix3 = [[4,5,6],[3,5,8],[2,4,1]]
test_vector_01 = [1, 2, 4]
test_vector_02 = [3, 1, 2]
n = 9

# add_vectors(test_vector_01,test_vector_02) should output [4,3,6]
print('Test Output for add_vectors: ' + str(add_vectors(test_vector_01,test_vector_02)))
print('Should have been [4,3,6]')

# scalar_vector(test_vector_01,n) should output [9,18,36]
print('Test Output for scalar_vector: ' + str(scalar_vector(test_vector_01,n)))
print('Should have been [9,18,36]')

# scalar_vector(test_vector_02,n) should output [27,9,18]
print('Test Output for scalar_vector: ' + str(scalar_vector(test_vector_02,n)))
print('Should have been [27,9,18]')

# scalar_matrix_mult(matrix1,n) should output [[9,18,27][36,45,54][63,72,81]]
print('Test Output for scalar_matrix_mult: ' + str(scalar_matrix_mult(matrix1,n)))
print('Should have been [[9,18,27][36,45,54][63,72,81]]')

# scalar_matrix_mult(matrix3,n) should output [[36,45,54][27,45,72][18,36,9]]
print('Test Output for scalar_matrix_mult: ' + str(scalar_matrix_mult(matrix3,n)))
print('Should have been [[36,45,54][27,45,72][18,36,9]]')

# add_matrix(matrix1,matrix2) should output [[7,10,7][7,10,7][14,10,18]]
print('Test Output for add_matrix: ' + str(add_matrix(matrix1,matrix2)))
print('should output [[7,10,7][7,10,7][14,10,18]]')

# add_matrix(matrix1,matrix3) should output [[5,7,9][7,10,14][9,12,10]]
print('Test Output for add_matrix: ' + str(add_matrix(matrix1,matrix3)))
print('should output [[5,7,9][7,10,14][9,12,10]]')

# matrix_vector_mult(matrix1,test_vector_01) should output [17,38,59]
print('Test Output for matrix_vector_mult: ' + str(matrix_vector_mult(matrix1,test_vector_01)))
print('Should have been [17,38,59]')

# matrix_vector_mult(matrix1,test_vector_02) should output [11,29,47]
print('Test Output for matrix_vector_mult: ' + str(matrix_vector_mult(matrix1,test_vector_02)))
print('Should have been [11,29,47]')

# matrix_matrix_mult(matrix1,matrix2) should output [[33,81,129][24,69,114][33,75,117]]
print('Test Output for matrix_matrix_mult: ' + str(matrix_matrix_mult(matrix1,matrix2)))
print('Should have been [[33,81,129][24,69,114][33,75,117]]')

# matrix_matrix_mult(matrix1,matrix3) should output [[16,43,70][27,69,111][25,70,115]]
print('Test Output for matrix_matrix_mult: ' + str(matrix_matrix_mult(matrix1,matrix3)))
print('Should have been [[16,43,70][27,69,111][25,70,115]]')