import unittest
from calculator import *

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


if __name__ == '__main__':
   unittest.main()