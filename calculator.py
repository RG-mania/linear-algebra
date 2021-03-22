

def det2x2(arr):
    if(len(arr) == 2 and len(arr[0]) == 2 and len(arr[1]) == 2):
        det = arr[0][0]*arr[1][1] - arr[0][1]*arr[1][0]
        return det
    else:
        raise NameError('Array is not 2x2')



def main():
    arr = [[1, 2], [3, 4]]
    print(det2x2(arr))

if __name__ == "__main__":
    main()