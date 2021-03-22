

def det2x2(arr):
    if(len(arr) == 2 and len(arr[0]) == 2 and len(arr[1]) == 2):
        det = arr[0][0]*arr[1][1] - arr[0][1]*arr[1][0]
        return det
    else:
        raise NameError('Array is not 2x2')

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

def main():
    arr = [[1, -4, 1, 2], [-1, 4, 4, 1], [3, 3, 3, 4], [2, 5, 2, -1]]
    arr2 = [[1, 1, 1], [2, -3, 3], [1, 2, 2]]
    print(detBareiss(arr))

if __name__ == "__main__":
    main()