#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import textwrap

from proteindf_tools.basis2 import Basis2


class TestBasis2(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_crete(self):
        bs2 = Basis2()

    def test_DZVP(self):
        expect = """
        O-DZVP.H
         2 0 0 0 0 0 0 0 0 0 0
         4
             5.099918e+01  9.660500e-03
             7.483218e+00  7.372890e-02
             1.777468e+00  2.958581e-01
             5.193295e-01  7.159053e-01
         1
             1.541100e-01  1.000000e+00
        """
        expect = textwrap.dedent(expect).strip()

        bs2 = Basis2()
        dzvp = bs2.get_basisset("O-DZVP.H")
        self.assertEqual(str(dzvp).strip(), expect)

    def test_DZVP2(self):
        expect = """
        O-DZVP2.C
         3 2 1 0 0 0 0 0 0 0 0
         7
             5.784157e+03  8.190000e-04
             8.693035e+02  6.293500e-03
             1.985116e+02  3.178120e-02
             5.642990e+01  1.172734e-01
             1.828546e+01  3.034763e-01
             6.448715e+00  4.535214e-01
             2.341860e+00  2.430591e-01
         2
             5.459533e+00 -7.780440e-02
             4.781968e-01  5.714947e-01
         1
             1.457301e-01  1.000000e+00
         5
             3.425856e+01  5.804300e-03
             7.863895e+00  4.064030e-02
             2.344519e+00  1.550219e-01
             7.961715e-01  3.531444e-01
             2.726804e-01  4.550062e-01
         1
             8.926050e-02  1.000000e+00
         1
             6.000000e-01  1.000000e+00
        """
        expect = textwrap.dedent(expect).strip()

        bs2 = Basis2()
        dzvp2 = bs2.get_basisset("O-DZVP2.C")
        self.assertEqual(str(dzvp2).strip(), expect)


if __name__ == "__main__":
    unittest.main()
