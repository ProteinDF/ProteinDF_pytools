#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from proteindf_tools.pdfmath import Math


class TestPdfMath(unittest.TestCase):
    def test_dbfact_zero(self):
        # Math.dbfact(0) returns 1.0 (0!! is mathematically 1.0).
        self.assertEqual(Math.dbfact(0), 1.0)

    def test_dbfact_known_values(self):
        # Known values for n!!:
        # 1!! = 1
        # 2!! = 2
        # 3!! = 3
        # 4!! = 8
        # 5!! = 15
        # 6!! = 48
        # 7!! = 105
        # 8!! = 384
        known_values = {
            1: 1.0,
            2: 2.0,
            3: 3.0,
            4: 8.0,
            5: 15.0,
            6: 48.0,
            7: 105.0,
            8: 384.0,
        }
        for n, expected in known_values.items():
            with self.subTest(n=n):
                res = Math.dbfact(n)
                self.assertEqual(res, expected)
                self.assertIsInstance(res, float)

    def test_dbfact_even_odd_boundaries(self):
        # Boundary for even numbers: 2 (loop while x > 2 does not execute)
        self.assertEqual(Math.dbfact(2), 2.0)
        # Boundary for odd numbers: 1 (loop while x > 1 does not execute)
        self.assertEqual(Math.dbfact(1), 1.0)

    def test_dbfact_negative(self):
        # (-1)!! is mathematically 1.0
        self.assertEqual(Math.dbfact(-1), 1.0)

    def test_dbfact_float(self):
        # Floating point inputs:
        # Math.dbfact casts input with int(x), truncating the fractional part towards zero.
        self.assertEqual(Math.dbfact(1.5), 1.0)  # int(1.5) == 1 -> returns 1.0
        self.assertEqual(Math.dbfact(2.5), 2.0)  # int(2.5) == 2 -> returns 2.0
        self.assertEqual(Math.dbfact(3.7), 3.0)  # int(3.7) == 3 -> returns 3.0
        self.assertEqual(Math.dbfact(-1.5), 1.0)  # int(-1.5) == -1 -> returns 1.0

    def test_dbfact_float_zero(self):
        # For inputs whose truncated integer part is 0 (0.5, -0.5), Math.dbfact returns 1.0.
        self.assertEqual(Math.dbfact(0.5), 1.0)
        self.assertEqual(Math.dbfact(-0.5), 1.0)


if __name__ == "__main__":
    unittest.main()
