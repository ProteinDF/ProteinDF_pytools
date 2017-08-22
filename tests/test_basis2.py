#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from pdfpytools.basis2 import Basis2

class TestBasis2(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
        
    def test_crete(self):
        bs2 = Basis2()

    def test_DZVP(self):
        bs2 = Basis2()
        dzvp = bs2.get_basisset("O-DZVP.H")
        print(dzvp)

    def test_DZVP2(self):
        bs2 = Basis2()
        dzvp2 = bs2.get_basisset("O-DZVP2.C")
        print(dzvp2)
        
if __name__ == "__main__":
    unittest.main()
    
