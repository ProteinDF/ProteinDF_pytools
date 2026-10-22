#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from proteindf_tools import constants
from proteindf_tools.pdfarchive_object import PdfArchiveObject
from proteindf_tools.qmsim import QmSim


class TestConstants(unittest.TestCase):
    def test_au2ang_value(self):
        self.assertEqual(constants.AU2ANG, 0.52917721067)

    def test_ang2au_value(self):
        self.assertEqual(constants.ANG2AU, 1.0 / constants.AU2ANG)

    def test_conversion_product(self):
        self.assertAlmostEqual(constants.AU2ANG * constants.ANG2AU, 1.0)

    def test_known_conversion(self):
        # 1.0 Angstrom is approximately 1.8897261 Bohr (ANG2AU)
        self.assertAlmostEqual(1.0 * constants.ANG2AU, 1.8897261, places=6)
        self.assertAlmostEqual(1.8897261 * constants.AU2ANG, 1.0, places=6)


class TestQmSim(unittest.TestCase):
    def test_default_method(self):
        qmsim = QmSim()
        self.assertEqual(qmsim.method, "hf")

    def test_set_method_string(self):
        qmsim = QmSim()
        qmsim.method = "dft"
        self.assertEqual(qmsim.method, "dft")

    def test_set_method_non_string(self):
        qmsim = QmSim()
        qmsim.method = 123
        self.assertEqual(qmsim.method, "123")
        qmsim.method = 3.14
        self.assertEqual(qmsim.method, "3.14")

    def test_instances_independent(self):
        q1 = QmSim()
        q2 = QmSim()
        q1.method = "b3lyp"
        self.assertEqual(q1.method, "b3lyp")
        self.assertEqual(q2.method, "hf")


class TestPdfArchiveObject(unittest.TestCase):
    def test_init_without_args(self):
        obj = PdfArchiveObject()
        self.assertIsInstance(obj, PdfArchiveObject)

    def test_init_with_rhs(self):
        """__init__(self, rhs=None) accepts rhs argument but does not use it."""
        obj = PdfArchiveObject(rhs="dummy")
        self.assertIsInstance(obj, PdfArchiveObject)
        obj2 = PdfArchiveObject(rhs=123)
        self.assertIsInstance(obj2, PdfArchiveObject)


if __name__ == "__main__":
    unittest.main()
