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

"""Assumes arr1, arr2 are 2d arrays representing matrices"""
def matricesEqual(arr1, arr2, tol = 0.1):
    if len(arr1) != len(arr2):
        return False
    for i in range(len(arr1)):
        if not vectorsEqual(arr1[i], arr2[i], tol):
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
        a1solve = linearSystem(arr1)
        self.assertEqual(a1solve, no_sol_msg)
        arr2 = [[1, 1, 1, 3], [2, 2, 3, 8], [1, 2, 2, 5]]
        a2solve = linearSystem(arr2)
        self.assertTrue(vectorsEqual(a2solve, [1, 0, 2]))
        arr3 = [[1, 1, 1, 3], [1, 1, 1, 3], [1, 2, 2, 5]]
        a3solve = linearSystem(arr3)
        self.assertEqual(a3solve, inf_sol_msg)

    def test_linear04_inf(self):
        arr1 = [[0, 0]]
        self.assertEqual(linearSystem(arr1), inf_sol_msg)
        arr2 = [[0, 1, 1, 3],
                [1, 1, 1, 4],
                [1, 2, 2, 5]]
        self.assertEqual(linearSystem(arr2), no_sol_msg)
        arr3 = [[0, 1, 1, 3],
                [0, 0, 2, 5],
                [1, 1, 1, 4]]
        a3solve = linearSystem(arr3)
        self.assertEqual(a3solve, [1, 0.5, 2.5])

    def test_rref01(self):
        arr1 = [[1, 1, 1, 3],
                [1, 1, 2, 4],
                [1, 1, 3, 5]]
        arr1sol = [[1, 1, 0, 2],
                   [0, 0, 1, 1],
                   [0, 0, 0, 0]]
        rref(arr1)
        self.assertTrue(matricesEqual(arr1, arr1sol))
        self.assertTrue(confirm_rref(arr1sol))
        self.assertEqual(rank(arr1), 2)

    def test_confirm_rref_01(self):
        arr1 = [[1, 0, 0],
                [0, 1, 0],
                [0, 0, 2]]
        arr2 = [[1, 2, 3],
                [2, 1, 3],
                [0, 0, 2]]
        arr3 = [[0, 2, 3],
                [0, 0, 0],
                [2, 0, 0]]
        arr4 = [[1, 0, 0],
                [0, 1, 0],
                [0, 0, 0]]
        arr5 = [[1, 2, 0],
                [0, 1, 0],
                [0, 0, 4]]
        arr6 = [[1, 2, -3],
                [0, 1, 3],
                [0, 0, 4]]
        arr7 = [[1, 0, 0],
                [0, 0, 0],
                [0, 0, 0]]
        arr8 = [[1, 0, 3, 0],
                [0, 2, 1, 0],
                [0, 0, 0, 1],
                [0, 0, 0, 0]]
        arr9 = [[1, 0, 0, 3],
                [0, 4, 0, 2],
                [0, 0, 3, 1],
                [0, 0, 0, 0]]
        self.assertTrue(confirm_rref(arr1))
        self.assertFalse(confirm_rref(arr2))
        self.assertFalse(confirm_rref(arr3))
        self.assertTrue(confirm_rref(arr4))
        self.assertFalse(confirm_rref(arr5))
        self.assertFalse(confirm_rref(arr6))
        self.assertTrue(confirm_rref(arr7))
        self.assertTrue(confirm_rref(arr8))
        self.assertTrue(confirm_rref(arr9))


    #----------------------------------------------------------------

    def test_create_identity01(self):
        i0 = []
        i1 = [[1]]
        i2 = [[1, 0],
              [0, 1]]
        i3 = [[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]]
        i4 = [[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]]
        i5 = [[1, 0, 0, 0, 0],
              [0, 1, 0, 0, 0],
              [0, 0, 1, 0, 0],
              [0, 0, 0, 1, 0],
              [0, 0, 0, 0, 1]]
        self.assertTrue(matricesEqual(i0, identity_mat(0)))
        self.assertTrue(matricesEqual(i1, identity_mat(1)))
        self.assertTrue(matricesEqual(i2, identity_mat(2)))
        self.assertTrue(matricesEqual(i3, identity_mat(3)))
        self.assertTrue(matricesEqual(i4, identity_mat(4)))
        self.assertTrue(matricesEqual(i5, identity_mat(5)))

    def test_dup_arr01(self):
        a1 = identity_mat(3)
        new_a1 = dup_mat(a1)
        self.assertTrue(matricesEqual(a1, new_a1))
        new_a1[1][0] = 2
        self.assertFalse(matricesEqual(a1, new_a1))

    #----------------------------------------------------------

    def test_invert_matrix01(self):
        arr1 = [[1, 2],
                [3, 4]]
        arr2 = [[2, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]
        arr3 = identity_mat(7)
        invert(arr1)
        invert(arr2)
        invert(arr3)
        arr1_inv = [[-2, 1],
                    [1.5, -0.5]]
        arr2_inv = [[0.5, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]]

        self.assertTrue(matricesEqual(arr1, arr1_inv))
        self.assertTrue(matricesEqual(arr2, arr2_inv))
        self.assertTrue(matricesEqual(arr3, identity_mat(7)))


    def test_multiply_01(self):
        arr1 = [[1, 2],
                [3, 4]]
        arr2 = [[2, 0],
                [0, 2]]
        prod12 = [[2, 4],
                  [6, 8]]
        self.assertTrue(matricesEqual(multiply(arr1, arr2), prod12))

    def test_rank01(self):
        arr1 = identity_mat(4)
        self.assertEqual(rank(arr1), 4)
        arr2 = identity_mat(1)
        self.assertEqual(rank(arr2), 1)

    def test_kernel01(self):
        nokernel = [[1, 0],
                    [0, 1]]
        self.assertTrue(matricesEqual(kernel(nokernel), []))
        arr01 = [[1, 0, 3],
                 [0, 1, 2],
                 [0, 0, 0]]
        ker01 = transpose([[-3, -2, 1]])
        self.assertTrue(matricesEqual(kernel(arr01), ker01))
        arr02 = [[1, 2, 3, 4],
                 [0, 1, 2, 3],
                 [0, 0, 0, 1]]
        ker02 = transpose([[1, -2, 1, 0]])
        self.assertTrue(matricesEqual(kernel(arr02), ker02))
        arr03 = [[1, 2, 2, -5, 6],
                 [-1, -2, -1, 1, -1],
                 [4, 8, 5, -8, 9],
                 [3, 6, 1, 5, -7]]
        ker03 = [[-2, -3, 4],
                 [1, 0, 0],
                 [0, 4, -5],
                 [0, 1, 0],
                 [0, 0, 1]]
        self.assertTrue(matricesEqual(kernel(arr03), ker03))

    def test_kernel02(self):
        arr1 = [[1, 2],
                [3, 4]]
        self.assertTrue(matricesEqual(kernel(arr1), []))
        arr2 = [[1, 2, 3]]
        ker2 = [[-2, -3],
                [1, 0],
                [0, 1]]
        self.assertTrue(matricesEqual(kernel(arr2), ker2))
        arr3 = [[0, 0],
                [0, 0]]
        ker3 = [[1, 0],
                [0, 1]]
        self.assertTrue(matricesEqual(kernel(arr3), ker3))
        arr4 = [[2, 3],
                [6, 9]]
        ker4 = transpose([[-1.5, 1]])
        self.assertTrue(matricesEqual(kernel(arr4), ker4))
        arr5 = [[1, 0, 5, 7],
                [0, 1, 2, 8]]
        ker5 = [[-5, -7],
                [-2, -8],
                [1,   0],
                [0,   1]]
        self.assertTrue(matricesEqual(kernel(arr5), ker5))

    def test_image01(self):
        arr1 = [[1, 2],
                [3, 4]]
        im1 = dup_mat(arr1)
        self.assertTrue(matricesEqual(image(arr1), im1))
        arr2 = [[1, 2, 3, 4],
                [0, 1, 2, 2]]
        im2 = [[1, 2],
               [0, 1]]
        self.assertTrue(matricesEqual(image(arr2), im2))
        arr3 = [[1, 2, 3]]
        im3 = [[1]]
        self.assertTrue(matricesEqual(image(arr3), im3))
        arr4 = [[1, 3, 2, 8],
                [0, 0, 3, 2]]
        im4 = [[1, 2],
               [0, 3]]
        self.assertTrue(matricesEqual(image(arr4), im4))

if __name__ == '__main__':
   unittest.main()