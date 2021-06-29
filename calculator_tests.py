import unittest
from calculator import *
from math import *

#error messages
inf_sol_msg = "System has an infinite solution set"
no_sol_msg = "System is inconsistent and has no solution"

"""Assumes v1, v2 are vectors"""
def vectorsEqual(v1, v2, tol = 0.1):
    if len(v1) != len(v2):
        return False
    if tol < 0:
        raise ArithmeticError("Tolerance must be greater than 0")
    for i in range(len(v1)):
        if not isclose(v1[i], v2[i], rel_tol = tol):
            return False
    return True

class TestCalculator(unittest.TestCase):

    def test_2x2det_01(self):
        arr1 = [[1, 2], [3, 4]]
        arr2 = [[0, 3], [0, 6]]
        arr3 = [[2, 1], [8, 4]]
        arr4 = [[15, 2], [4, 3]]
        arr5 = [[2, 2], [4, 6]]
        self.assertEqual(det2x2(arr1), -2)
        self.assertEqual(det2x2(arr2), 0)
        self.assertEqual(det2x2(arr3), 0)
        self.assertEqual(det2x2(arr4), 37)
        self.assertEqual(det2x2(arr5), 4)

    def test_detBareiss_01(self):
        arr1 = [[1, 1, 1], [2, -3, 3], [1, 2, 2]]
        arr2 = [[0, 3], [0, 6]]
        arr3 = [[2, 1], [8, 4]]
        arr4 = [[15, 2], [4, 3]]
        arr5 = [[2, 2], [4, 6]]
        arr6 = [[1, 2], [3, 4]]
        self.assertEqual(detBareiss(arr1), -6)
        self.assertEqual(detBareiss(arr2), 0)
        self.assertEqual(detBareiss(arr3), 0)
        self.assertEqual(detBareiss(arr4), 37)
        self.assertEqual(detBareiss(arr5), 4)
        self.assertEqual(detBareiss(arr6), -2)

    def test_detBareiss_02(self):
        arr1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        arr2 = [[1, 2, 3], [4, 5, 6], [7, 8, 3]]
        self.assertEqual(detBareiss(arr1), 0)
        self.assertEqual(detBareiss(arr2), 18)

    def test_detBareiss_4x4(self):
        arr1 = [[1, 2, 3, 1], [4, 5, 6, 7], [7, 8, 3, 3], [3, 7, 8, 1]]
        arr2 = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 7]]
        arr3 = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 2]]
        self.assertEqual(detBareiss(arr1), 72)
        self.assertEqual(detBareiss(arr2), 0)
        self.assertEqual(detBareiss(arr3), 180)

    def test_detBareiss_5x5(self):
        arr1 = [[1, 2, 3, 4, 2],
                [5, 6, 7, 8, 1],
                [9, 1, 2, 3, 3],
                [4, 5, 6, 2, 5],
                [2, 7, 8, 9, 6]]

        arr2 = [[1, 2, 3, 4, 2],
                [0, 6, 7, 8, 1],
                [0, 0, 2, 3, 3],
                [0, 0, 0, 2, 5],
                [0, 0, 0, 0, 0]]

        arr3 = [[1, 2, 3, 0, 0],
                [5, 6, 7, 0, 0],
                [3, 1, 2, 0, 0],
                [4, 5, 6, 2, 5],
                [2, 7, 8, 3, 6]]

        arr4 = [[1, 2, 0, 4, 2],
                [5, 6, 0, 8, 1],
                [9, 1, 0, 3, 3],
                [4, 5, 0, 2, 5],
                [2, 7, 0, 9, 6]]

        arr5 = [[1, 2, 3, 4, 2],
                [5, 6, 7, 8, 1],
                [9, 1, 2, 3, 3],
                [0, 0, 0, 0, 0],
                [2, 7, 8, 9, 6]]

        self.assertEqual(detBareiss(arr1), 1005)
        self.assertEqual(detBareiss(arr2), 0)
        self.assertEqual(detBareiss(arr3), 36)
        self.assertEqual(detBareiss(arr4), 0)
        self.assertEqual(detBareiss(arr5), 0)


    #-------------------------------------------------------
    #Linear System tests
    
    def test_linear01_one_eq(self):
        #test with one vector where var is 1 - result should match input
        arr1 = [[1, 3]]
        self.assertTrue(vectorsEqual(linearSystem(arr1), [3]))
        arr2 = [[1, 4]]
        self.assertTrue(vectorsEqual(linearSystem(arr2), [4]))
        arr3 = [[1, 20]]
        self.assertTrue(vectorsEqual(linearSystem(arr3), [20]))

    def test_linear02(self):
        arr1 = [[3, 1, 9], [2, 1, 7]]
        a1solve = linearSystem(arr1)
        self.assertTrue(vectorsEqual(a1solve, [2.0, 3.0]))
        arr2 = [[1, 1, 1, 3], [2, -3, 3, 8], [1, 2, 2, 5]]
        a2solve = linearSystem(arr2)
        self.assertTrue(vectorsEqual(a2solve, [1, 0, 2]))

    def test_linear03(self):
        arr1 = [[1, 1, 1, 3], [1, 1, 1, 4], [1, 2, 2, 5]]
        self.assertEqual(linearSystem(arr1), no_sol_msg)
        arr2 = [[1, 1, 1, 3], [2, 2, 3, 8], [1, 2, 2, 5]]
        a2solve = linearSystem(arr2)
        self.assertTrue(vectorsEqual(a2solve, [1, 0, 2]))
        #breaks on infinite solutions
        arr3 = [[1, 1, 1, 3], [1, 1, 1, 3], [1, 2, 2, 5]]
        self.assertEqual(linearSystem(arr3), inf_sol_msg)




if __name__ == '__main__':
   unittest.main()