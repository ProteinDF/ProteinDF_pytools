#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

import proteindf_bridge as bridge
from proteindf_tools.gauparam import GaussianParam


class TestGaussianParam(unittest.TestCase):
    def test_default_values(self):
        p = GaussianParam()
        self.assertEqual(p.title, "")
        self.assertEqual(p.jobtype, "sp")
        self.assertEqual(p.charge, 0)
        self.assertEqual(p.method, "hf")
        self.assertEqual(p.basisset, "6-31G")
        self.assertEqual(p.ext_basisset, "")

    def test_setters(self):
        p = GaussianParam()

        p.title = "Sample Molecule"
        self.assertEqual(p.title, "Sample Molecule")

        p.jobtype = "opt"
        self.assertEqual(p.jobtype, "opt")

        p.charge = -1
        self.assertEqual(p.charge, -1)
        p.charge = "2"
        self.assertEqual(p.charge, 2)

        p.method = "b3lyp"
        self.assertEqual(p.method, "b3lyp")

        p.basisset = "cc-pVDZ"
        self.assertEqual(p.basisset, "cc-pVDZ")

        p.ext_basisset = "DZ"
        self.assertEqual(p.ext_basisset, "DZ")

    def test_get_inputfile_contents_h2(self):
        p = GaussianParam()
        p.title = "H2 calculation"
        p.jobtype = "sp"
        p.charge = 0
        p.method = "rhf"
        p.basisset = "6-31G"

        mol = bridge.AtomGroup()
        h1 = bridge.Atom()
        h1.symbol = "H"
        h1.xyz = bridge.Position(0.0, 0.0, 0.0)
        mol.set_atom("1", h1)

        h2 = bridge.Atom()
        h2.symbol = "H"
        h2.xyz = bridge.Position(0.0, 0.0, 0.74)
        mol.set_atom("2", h2)

        p.molecule = mol

        contents = p.get_inputfile_contents()
        self.assertIn("#P rhf/6-31G sp", contents)
        self.assertIn("5D GFPrint GFInput NoSymm", contents)
        self.assertIn("H2 calculation", contents)
        self.assertIn("0 1", contents)
        self.assertIn("H 0.0 0.0 0.0", contents)
        self.assertIn("H 0.0 0.0 0.74", contents)

    def test_get_inputfile_contents_gen_basisset(self):
        p = GaussianParam()
        p.title = "H2 DZ test"
        p.method = "b3lyp"
        p.basisset = "gen"
        p.ext_basisset = "DZ"

        mol = bridge.AtomGroup()
        h1 = bridge.Atom()
        h1.symbol = "H"
        h1.xyz = bridge.Position(0.0, 0.0, 0.0)
        mol.set_atom("1", h1)

        p.molecule = mol

        contents = p.get_inputfile_contents()
        self.assertIn("#P b3lyp/gen sp", contents)
        self.assertIn("H 0\n S 3", contents)
        self.assertIn("****", contents)


if __name__ == "__main__":
    unittest.main()
