from math import *

#TODO:
#UNITTESTING!!!!!!!!!!!!!!!!
#Solve systems of equations - DONE ( I think, needs more testing )
#Gram-Schmidt Orthogonalization
#Finding image/kernel of matrix - DONE ( I think, needs more testing )
#Finding rank of matrix - DONE ( I think, needs more testing )
#Find inverse of matrix - DONE ( I think, needs more testing )
#Find eigenvalues / eigenvectors of 2x2 and 3x3 matrices (if they exist)
#Implement GUI
#GCD and LCM algorithms, because why not
#root finding algorithms for polynomials
#square root and more? not sure

"""Assumes matrix is well formed, and rows and columns have consistent length"""
#  Prints matrices in a nicely-formatted manner
def printmatrix(mat):
    if len(mat) == 0:
        print("[]")
        return

    dashmult = 4 #controls length of dashlines relative to matrix size
    h_len = len(mat[0])

    #print pre-matrix dashline
    print(" ", end = '')
    for i in range(h_len*dashmult + 4):
        print("-", end = '') #no newline
    print()

    #print matrix
    for i in range(len(mat)): #each row
        print(" |", end = '')
        for j in range(h_len): #each individual element
            print(f'{mat[i][j]:>4}', end = '')
        print("  | ")

    #print post-matrix dashline
    print(" ", end = '')
    for i in range(h_len*dashmult + 4):
        print("-", end = '') #no newline
    print()


"""Assumes matrix is 2x2"""
def det2x2(arr):
    if(len(arr) == 2 and len(arr[0]) == 2 and len(arr[1]) == 2):
        det = arr[0][0]*arr[1][1] - arr[0][1]*arr[1][0]
        return det
    else:
        raise ValueError('Array is not 2x2')

"""Assumes matrix is square"""
def detBareiss(arr):
    pivot = 1
    arrsize = len(arr)
    negmult = 1
    for k in range(arrsize-1):
        if pivot == 0:
            return 0
        for i in range(k+1, arrsize):
            for j in range(k+1, arrsize):
                #print("k, i, j: ", k, i, j)
                arr[i][j] = (arr[i][j]*arr[k][k] - arr[i][k]*arr[k][j])//pivot
        pivot = arr[k][k]
        p = k+2
        while (arr[k+1][k+1] == 0 and p<arrsize):
            #print('swapping rows')
            if (arr[p][k+1] != 0):
                temp = arr[k+1]
                arr[k+1] = arr[p]
                arr[p] = temp
                # print(arr[k+1])
                # print(arr[p])
                negmult *= -1
            else:
                p+=1 

    return arr[arrsize-1][arrsize-1]*negmult

"""Transforms a matrix into reduced row echelon form - needs testing"""
def rref(arr):
    divisor = 1
    arrsize = len(arr) #number of equations
    numvars = len(arr[0]) - 1
    poffset = 0
    for k in range(arrsize):
        # Below: Swapping rows, incrementing poffset if necessary (Happens when next divisor to-be is 0)
        while (k+poffset < numvars and arr[k][k+poffset] == 0): # not on last column, next divisor would be 0
            p = k + 1 #row below - to swap
            if k+poffset == numvars: #??? Should never reach here
                return "System is inconsistent and has no solution"
            while (arr[k][k+poffset] == 0 and p<arrsize):
                if (arr[p][k+poffset] != 0): # if row to swap has valid divisor, then swap
                    # printmatrix(arr)
                    # print('swapping rows')
                    temp = arr[k]
                    arr[k] = arr[p]
                    arr[p] = temp
                else:
                    p += 1
            if arr[k][k+poffset] == 0: # if we haven't found a suitable swap, increment the poffset to next column
                                           # should only happen when we have a column of 0s below divisor to-be
                poffset += 1

        # if we have a row of 0s?? Should only happen when we're gonna end
        if k+poffset == numvars:
            break # row of 0s, break and end
            # print("Poffset is", poffset)
            # print("K is", k)
            # if arr[k][k+poffset] != 0:
            #     break # break out of for loop
            #     return "System is inconsistent and has no solution"
            # else:
            #     poffset -= 1
            #     break
        
        for i in range(arrsize):
            if i != k: # don't modify row of current pivot
                for j in range(numvars+1):
                    if j != k+poffset: # don't modify column of current pivot
                                    #print(arr[i][k + poffset])
                        arr[i][j] = (arr[i][j]*arr[k][k+poffset] - arr[i][k+poffset]*arr[k][j])//divisor
                                    #print('arr[{i}][{j}] = {val}'.format(i=i, j=j, val = arr[i][j]))
        for l in range(arrsize): #set the column of current pivot to 0
            if l != k:
                arr[l][k+poffset] = 0

        divisor = arr[k][k+poffset]

"""Note: The following method assumes that the number of equations is equal to the number of variables"""
#Edit: it might not - still have to test
"""Also assumes that the matrix is well-formed - # of rows and columns is consistent"""
def linearSystem(arr):
    arrsize = len(arr) #number of equations
    numvars = len(arr[0]) - 1

    rref(arr)
    #matrix should now be in reduced row-echelon form
    finVector = []
    for x in range(arrsize):
        div = arr[x][x]
        ans = arr[x][-1]
        #print('ans: ', ans, 'div: ', div)
        p = 1
        while div == 0 and (x+p)<numvars-1:
            div = arr[x][x+p]
            p += 1
        if div == 0:
            if ans == 0:
                return "System has an infinite solution set"
            else:
                return "System is inconsistent and has no solution"
        else:
            finVector.append(ans/div)
    return finVector


"""Creates identity matrix given size"""
def identity_mat(size):
    if size == 0:
        return []
    arr = [[0 for i in range(size)] for j in range(size)]
    for i in range(size):
        arr[i][i] = 1
    return arr

"""Duplicates array"""
def dup_mat(arr):
    new_arr = [None] * len(arr)
    for i in range(len(arr)):
        new_arr[i] = arr[i].copy()
    return new_arr

"""Returns true if the matrix is square, false otherwise"""
def isSquare(arr):
    numRows = len(arr)
    for i in range(numRows):
        if len(arr[i]) != numRows:
            return False
    return True

"""Returns the inverse of a matrix, assumes it is square"""
def invert(arr):
    arrsize = len(arr)

    if not isSquare(arr):
        raise ValueError("Matrix is not square")
    determinant = detBareiss(dup_mat(arr)) #technically doing double work here, but it's clearer
    if determinant == 0:
        raise ValueError("Matrix is not invertible (determinant 0)")

    #append an identity matrix to the right
    identity = identity_mat(len(arr))
    for i in range(len(arr)):
        arr[i].extend(identity[i])
    rref(arr)

    #chop off unused left size, divide the whole thing by determinant
    for i in range(len(arr)):
        arr[i] = arr[i][-arrsize:]
        for j in range(arrsize):
            arr[i][j] = arr[i][j] / determinant
            #convert to int if possible, to make it look nice
            if arr[i][j].is_integer():
                arr[i][j] = int(arr[i][j])

    return arr

"""Find the rank of a matrix"""
def rank(arr):
    rref(arr)
    rank = 0
    for i in range(len(arr)):
        for j in range(len(arr[0])):
            if arr[i][j] != 0:
                rank += 1
                break
    return rank

"""Find the kernel of a matrix"""
def kernel(arr):
    if len(arr) == 0:
        return []
    rref(arr)
    kernel = []
    arrlen = len(arr)
    redun_cols = [None]*len(arr[0]) #stores redundant columns and their position
    # if arr[arrlen-1][len(arr[0])-1] != 0: #aka determinant is not 0, so kernel is empty
    #     return kernel
    # # we know the kernel is not 0??????????????
    offset = 0
    for i in range(arrlen):

        while i+offset < len(arr[0]) and arr[i][i+offset] == 0: #for every column with no leading 1
            kernel.append([0]*len(arr[0])) #initialize empty vector
            redun_cols[i+offset] = offset # mark this column as redundant
            kernel[offset][i+offset] = 1
            offset += 1

    # in case there aren't enough equations, we still need to find other redundant columns w/ no 0 below them
    # only kicks in for a situation such as [[1, 0, 5, 7]
    #                                        [0, 1, 2, 8]]
    for i in range(arrlen+offset, len(arr[0])):
        kernel.append([0]*len(arr[0])) #initialize empty vector
        redun_cols[i] = offset
        kernel[offset][i] = 1
        offset += 1

    offset = 0
    for i in range(arrlen):
        leading_col = None
        for j in range(len(arr[0])):
            if leading_col != None and redun_cols[j] != None: # this is a redundant column after the leading 1
                temp = -arr[i][j]/leader
                if temp.is_integer():
                    temp = int(temp)
                kernel[redun_cols[j]][leading_col] = temp
                # note - dividing by the leader presents a possibility of ending up with non-integers in
                # the kernel. There is a more elegant solution here that involves scaling temp and the leader,
                # but it would likely require implementing a LCM algorithm, which would require more computational
                # complexity. No plans on implementing this feature unless deemed necessary. 
            if leading_col == None and arr[i][j] != 0:
                leading_col = j
                leader = arr[i][j]

    if len(kernel) > 0:
        kernel = transpose(kernel)
    return kernel

"""Returns the basis of the image of the matrix"""
def image(arr):
    ker_arr = dup_mat(arr)
    rref(ker_arr)
    offset = 0
    non_redun_cols = []
    for i in range(len(arr)):
        while i+offset < len(ker_arr[0]) and ker_arr[i][i+offset] == 0:
            offset += 1
        if i+offset < len(ker_arr[0]):
            non_redun_cols.append(i+offset)
    #at this point we have the list of all the non-redundant columns
    #we just need to eliminate the redundant ones from the original matrix

    for i in range(len(arr)):
        new_row = []
        for j in non_redun_cols:
            new_row.append(arr[i][j])
        arr[i] = new_row
    return arr

"""Transposes a matrix"""
def transpose(arr):
    newarr = []
    if len(arr) == 0:
        return newarr
    for column in range(len(arr[0])):
        newarr.append([])
        for row in range(len(arr)):
            newarr[column].append(arr[row][column])
    return newarr

"""Assuems matrix columns and rows have consistent lengths"""
def multiply(arr1, arr2):
    if len(arr1[0]) != len(arr2):
        raise ValueError("Matrices cannot be multiplied (dimensions do not match")
    #initialize product array
    prod = [[0 for i in range(len(arr1))] for j in range(len(arr2[0]))]
    #multiply matrices
    for i in range(len(prod)):
        for j in range(len(prod[0])):
            for k in range(len(arr1[0])):
                prod[i][j] += arr1[i][k] * arr2[k][j]
    return prod


"""Checks if matrix is in reduced row echelon form"""
def confirm_rref(arr):
    poffset = 0
    for i in range(len(arr)):
        j = 0
        while j < len(arr[i]) and arr[i][j] == 0:
            j += 1
        if j < i+poffset:
            return False
        else:
            poffset = j-i-1 #-1 to account for i being incremented next loop
            #check for all 0s in column with a leading number
            for x in range(len(arr)):
                if x != i and j<len(arr[x]) and arr[x][j] != 0:
                    return False
    return True

"""Performs Gram-Schmidt orthogonalization on matrix"""
def GSorthog(arr):
    pass

def dot(v1, v2):
    if not len(v1) == len(v2):
        raise ValueError("Vector lengths do not match")
    else:
        sum = 0
        for i in v1:
            sum += v1[i] * v2[i]
        return sum

"""Finds the magnitude of a vector"""
def magnitude(vec):
    sum = 0
    for i in vec:
        sum += vec[i]*vec[i]
    return sqrt(vec)

"""Finds the greatest common divisor of two numbers using the Euclidean algorithm"""
def gcd(a, b):
    if a == 0:
        return b
    if b == 0:
        return a
    else:
        q = int(floor(a/b))
        r = int(a - b*q)
        print("q: ", q, " r: ", r)
        return gcd(b, r)

"""Scales a vector by a scalar and returns the result"""
def scale(vec, scalar):
    newvec = []
    for i in vec:
        newvec.append(vec[i]*scalar)
    return newvec


def main():
    # arr = [[1, -4, 1, 2], [-1, 4, 4, 1], [3, 3, 3, 4], [2, 5, 2, -1]]
    # arr2 = [[1, 1, 1], [2, -3, 3], [1, 2, 2]]
    # sysEq = [[1, 1, 1, 3], [2, -3, 3, 8], [1, 2, 2, 5]]
    # badSys = [[1, 1, 1, 3], [1, 1, 1, 4], [1, 2, 2, 5]]
    # mat1 = [[1, 2], [3, 4]]
    # mat5 = [[1, 2, 3, 4, 5], [0,0,0,0,0]]
    # arr3 = [[3, 1, 9], [2, 1, 7]]
    # print(detBareiss(arr2))
    # print(linearSystem(arr3))
    arr = [[1, 2, 3, 4]]
    done = False
    while (not done):
        print("Greatest common divisor, enter two integers: ")
        int1 = float(input("int1: "))
        if int1 == 66:
            done = True
        int2 = float(input("int2: "))
        print ("result: ", gcd(int1, int2))



if __name__ == "__main__":
    main()