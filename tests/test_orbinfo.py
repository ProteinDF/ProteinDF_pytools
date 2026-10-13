#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

import proteindf_bridge as bridge
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.orbinfo import OrbInfo
from proteindf_tools.pdfparam import PdfParam


class TestOrbInfo(unittest.TestCase):
    def _create_h2_param(self):
        """Create minimal PdfParam representing an H2 molecule with STO-3G-like basis."""
        param = PdfParam()

        mol = bridge.AtomGroup()
        atom1 = bridge.Atom()
        atom1.symbol = "H"
        atom1.label = "1"
        mol.set_atom("1", atom1)

        atom2 = bridge.Atom()
        atom2.symbol = "H"
        atom2.label = "2"
        mol.set_atom("2", atom2)

        param.molecule = mol

        # BasisSet for H (1 s-type CGTO)
        bs = BasisSet()
        bs.name = "sto-3g-H"
        cgto = ContractedGTO("s", 1)
        cgto[0] = PrimitiveGTO(3.42525, 0.15433)
        bs.append(cgto)

        param.set_basisset("H@1", bs)
        param.set_basisset("H@2", bs)

        return param

    def _create_hetero_param(self):
        """Create minimal PdfParam with H (s) and C (s, p) to test multiple shell types."""
        param = PdfParam()

        mol = bridge.AtomGroup()
        atom_h = bridge.Atom()
        atom_h.symbol = "H"
        atom_h.label = "1"
        mol.set_atom("1", atom_h)

        atom_c = bridge.Atom()
        atom_c.symbol = "C"
        atom_c.label = "2"
        mol.set_atom("2", atom_c)

        param.molecule = mol

        # BasisSet for H: 1 s orbital
        bs_h = BasisSet()
        bs_h.name = "basis-H"
        cgto_h_s = ContractedGTO("s", 1)
        cgto_h_s[0] = PrimitiveGTO(3.42525, 0.15433)
        bs_h.append(cgto_h_s)

        # BasisSet for C: 1 s orbital + 1 p orbital (3 components: px, py, pz)
        bs_c = BasisSet()
        bs_c.name = "basis-C"
        cgto_c_s = ContractedGTO("s", 1)
        cgto_c_s[0] = PrimitiveGTO(71.6168, 0.15433)
        cgto_c_p = ContractedGTO("p", 1)
        cgto_c_p[0] = PrimitiveGTO(2.9412, 0.15592)
        bs_c.append(cgto_c_s)
        bs_c.append(cgto_c_p)

        param.set_basisset("H@1", bs_h)
        param.set_basisset("C@2", bs_c)

        return param

    def test_h2_orbinfo(self):
        param = self._create_h2_param()
        orb = OrbInfo(param)

        self.assertEqual(orb.get_num_of_orbitals(), 2)

        # Orbital 0: H@1, s shell, s basis
        self.assertEqual(orb.get_atom_id(0), 0)
        self.assertEqual(orb.get_shell_type(0), "s")
        self.assertEqual(orb.get_basis_type(0), "s")
        atom0 = orb.get_atom(0)
        self.assertIsNotNone(atom0)
        self.assertEqual(atom0.symbol, "H")

        # Orbital 1: H@2, s shell, s basis
        self.assertEqual(orb.get_atom_id(1), 1)
        self.assertEqual(orb.get_shell_type(1), "s")
        self.assertEqual(orb.get_basis_type(1), "s")
        atom1 = orb.get_atom(1)
        self.assertIsNotNone(atom1)
        self.assertEqual(atom1.symbol, "H")

    def test_hetero_orbinfo(self):
        param = self._create_hetero_param()
        orb = OrbInfo(param)

        # Total orbitals: H has 1 (s), C has 1 (s) + 3 (p) = 5
        self.assertEqual(orb.get_num_of_orbitals(), 5)

        # Orbital 0: H (s)
        self.assertEqual(orb.get_atom_id(0), 0)
        self.assertEqual(orb.get_shell_type(0), "s")
        self.assertEqual(orb.get_basis_type(0), "s")

        # Orbital 1: C (s)
        self.assertEqual(orb.get_atom_id(1), 1)
        self.assertEqual(orb.get_shell_type(1), "s")
        self.assertEqual(orb.get_basis_type(1), "s")

        # Orbitals 2-4: C (p) -> px, py, pz
        expected_p_basis = ["px", "py", "pz"]
        for idx, expected_basis in enumerate(expected_p_basis, start=2):
            with self.subTest(orb_index=idx):
                self.assertEqual(orb.get_atom_id(idx), 1)
                self.assertEqual(orb.get_shell_type(idx), "p")
                self.assertEqual(orb.get_basis_type(idx), expected_basis)

    def test_out_of_range(self):
        param = self._create_h2_param()
        orb = OrbInfo(param)

        self.assertIsNone(orb.get_atom_id(99))
        self.assertIsNone(orb.get_atom(99))
        self.assertIsNone(orb.get_shell_type(99))
        self.assertIsNone(orb.get_basis_type(99))

    def test_empty_orbinfo(self):
        orb = OrbInfo()
        self.assertEqual(orb.get_num_of_orbitals(), 0)
        self.assertIsNone(orb.get_atom_id(0))
        self.assertIsNone(orb.get_atom(0))
        self.assertIsNone(orb.get_shell_type(0))
        self.assertIsNone(orb.get_basis_type(0))


if __name__ == "__main__":
    unittest.main()
