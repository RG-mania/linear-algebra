#TODO:
#Solve systems of equations
#Gram-Schmidt Orthogonalization
#Finding image/kernel of matrix
#Finding rank of matrix
#Find inverse of matrix
#Find eigenvalues / eigenvectors of 2x2 and 3x3 matrices (if they exist)

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
                arr[i][j] = (arr[i][j]*arr[k][k] - arr[i][k]*arr[k][j])/pivot
        pivot = arr[k][k]
        p = k+2
        while (arr[k+1][k+1] == 0 and p<arrsize):
            print('swapping rows')
            temp = arr[k+1]
            arr[k+1] = arr[p]
            arr[p] = temp
            print(arr[k+1])
            print(arr[p])
            negmult *= -1
            p+=1

    return arr[arrsize-1][arrsize-1]*negmult

"""Note: The following method assumes that the number of equations is equal to the number of variables"""
#Edit: it might not - still have to test
def sysEquationsBareiss(arr):
    divisor = 1
    arrsize = len(arr) #number of equations
    numvars = len(arr[0]) - 1
    poffset = 0
    for k in range(2):
        for i in range(arrsize):
            if i != k:
                for j in range(numvars+1):
                    if j != k+poffset:
                        #print(arr[i][k + poffset])
                        arr[i][j] = (arr[i][j]*arr[k][k+poffset] - arr[i][k+poffset]*arr[k][j])/divisor
                        #print('arr[{i}][{j}] = {val}'.format(i=i, j=j, val = arr[i][j]))
        for l in range(arrsize):
            if l != k:
                arr[l][k+poffset] = 0

        divisor = arr[k][k+poffset]
        while (k+1+poffset < numvars and arr[k+1][k+1+poffset] == 0):
            p = k + 2
            if k+poffset == numvars:
                return "System is inconsistent and has no solution"
            while (arr[k+1][k+1+poffset] == 0 and p<arrsize):
                if (arr[p][k+1+poffset] != 0):
                    print('swapping rows')
                    temp = arr[k+1]
                    arr[k+1] = arr[p]
                    arr[p] = temp
                else:
                    p += 1
            if arr[k+1][k+1+poffset] == 0:
                poffset += 1


        if k+1+poffset == numvars:
            if arr[k+1][k+1+poffset] != 0:
                return "System is inconsistent and has no solution"
            else:
                poffset -= 1

    #
    # for h in range(arrsize):
    #     print(arr[h][-4])

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
    # print(detBareiss(arr2))
    print(sysEquationsBareiss(badSys))

if __name__ == "__main__":
    main()