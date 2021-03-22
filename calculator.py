

def det2x2(arr):
    if(len(arr) == 2 and len(arr[0]) == 2 and len(arr[1]) == 2):
        det = arr[0][0]*arr[1][1] - arr[0][1]*arr[1][0]
        return det
    else:
        raise NameError('Array is not 2x2')

def detBareiss(arr):
    pivot = 1
    arrsize = len(arr)
    for k in range(arrsize-1):
        for i in range(k+1, arrsize):
            for j in range(k+1, arrsize):
                #print("k, i, j: ", k, i, j)
                arr[i][j] = (arr[i][j]*arr[k][k] - arr[i][k]*arr[k][j])/pivot
        pivot = arr[k][k]
    return arr[arrsize-1][arrsize-1]

def main():
    arr = [[1, 1, 1], [2, -3, 3], [1, 2, 2]]
    print(detBareiss(arr))

if __name__ == "__main__":
    main()