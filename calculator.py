#TODO:
#UNITTESTING!!!!!!!!!!!!!!!!
#Solve systems of equations - DONE ( I think, needs more testing )
#Gram-Schmidt Orthogonalization
#Finding image/kernel of matrix
#Finding rank of matrix - DONE ( I think, needs more testing )
#Find inverse of matrix - DONE ( I think, needs more testing )
#Find eigenvalues / eigenvectors of 2x2 and 3x3 matrices (if they exist)
#Implement GUI

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
    determinant = detBareiss(dup_mat(arr))
    if determinant == 0:
        raise ValueError("Matrix is not invertible (determinant 0)")

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
    redun_cols = [None]*len(arr[0])
    # if arr[arrlen-1][len(arr[0])-1] != 0: #aka determinant is not 0, so kernel is empty
    #     return kernel
    # # we know the kernel is not 0??????????????
    offset = 0
    for i in range(arrlen):
        #BROKEN RN, NEED TO FIX
        if arr[i][i+offset] == 0: # we have a column with no leading 1
            kernel.append([0]*len(arr[0])) #initialize empty vector
            redun_cols[i+offset] = offset # mark this column as redundant
            # for k in range(arrlen):
            #     kernel[offset].append(-(arr[k][i+offset]))
            # kernel[offset][i] = 1
            offset += 1
    offset = 0
    # for i in range(arrlen):
    #     for j in range(len(arr[0])):
    #         if redun_cols[j] != None: # this is a redundant column
    #             kernel[redun_cols[j]][i] = - 
    if len(kernel) > 0:
        kernel = transpose(kernel)
    print(redun_cols)
    return kernel

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



if __name__ == "__main__":
    main()