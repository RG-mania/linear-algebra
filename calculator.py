#TODO:
#UNITTESTING!!!!!!!!!!!!!!!!
#Solve systems of equations
#Gram-Schmidt Orthogonalization
#Finding image/kernel of matrix
#Finding rank of matrix
#Find inverse of matrix
#Find eigenvalues / eigenvectors of 2x2 and 3x3 matrices (if they exist)

"""Assumes matrix is well formed, and rows and columns have consistent length"""
#  Prints matrices in a nicely-formatted manner
def printmatrix(mat):
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
        raise NameError('Array is not 2x2')

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

"""Note: The following method assumes that the number of equations is equal to the number of variables"""
#Edit: it might not - still have to test
"""Also assumes that the matrix is well-formed - # of rows and columns is consistent"""
def linearSystem(arr):
    divisor = 1
    arrsize = len(arr) #number of equations
    numvars = len(arr[0]) - 1
    poffset = 0
    for k in range(arrsize):
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
        # Below: Swapping rows, incrementing poffset if necessary (Happens when next divisor to-be is 0)
        # SWAP TO THE START OF FOR LOOP, INSTEAD OF END OF PREVIOUS (case where arr[0][0] == 0)
        while (k+1+poffset < numvars and arr[k+1][k+1+poffset] == 0): # not on last column, next divisor would be 0
            p = k + 2 #row below - to swap
            if k+poffset == numvars: #??? Should never reach here
                return "System is inconsistent and has no solution"
            while (arr[k+1][k+1+poffset] == 0 and p<arrsize):
                if (arr[p][k+1+poffset] != 0): # if row to swap has valid divisor, then swap
                    print('swapping rows')
                    temp = arr[k+1]
                    arr[k+1] = arr[p]
                    arr[p] = temp
                else:
                    p += 1
            if arr[k+1][k+1+poffset] == 0: # if we haven't found a suitable swap, increment the poffset to next column
                                           # should only happen when we have a column of 0s below divisor to-be
                poffset += 1

            # if we have a row of 0s?? Should only happen when we're gonna end
            if k+1+poffset == numvars:
                # print("Poffset is", poffset)
                # print("K is", k)
                if arr[k+1][k+1+poffset] != 0:
                    return "System is inconsistent and has no solution"
                else:
                    poffset -= 1


    #matrix should now be in row-echelon form
    finVector = []
    for x in range(arrsize):
        div = arr[x][x]
        ans = arr[x][-1]
        #print('ans: ', ans, 'div: ', div)
        p = 1
        while div == 0 and p<numvars-1:
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

def main():
    arr = [[1, -4, 1, 2], [-1, 4, 4, 1], [3, 3, 3, 4], [2, 5, 2, -1]]
    arr2 = [[1, 1, 1], [2, -3, 3], [1, 2, 2]]
    sysEq = [[1, 1, 1, 3], [2, -3, 3, 8], [1, 2, 2, 5]]
    badSys = [[1, 1, 1, 3], [1, 1, 1, 4], [1, 2, 2, 5]]
    mat1 = [[1, 2], [3, 4]]
    mat5 = [[1, 2, 3, 4, 5], [0,0,0,0,0]]
    arr3 = [[3, 1, 9], [2, 1, 7]]
    # print(detBareiss(arr2))
    print(linearSystem(arr3))

if __name__ == "__main__":
    main()