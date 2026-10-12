#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import unittest

import proteindf_bridge as bridge
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_object import PdfParamObject, error_threshold


class TestPdfParamModel(unittest.TestCase):
    def _create_h2_param(self):
        """Create minimal PdfParam representing an H2 molecule with STO-3G-like basis."""
        param = PdfParam()

        mol = bridge.AtomGroup()
        atom1 = bridge.Atom()
        atom1.symbol = "H"
        atom1.label = "1"
        atom1.xyz = [0.0, 0.0, 0.0]
        mol.set_atom("1", atom1)

        atom2 = bridge.Atom()
        atom2.symbol = "H"
        atom2.label = "2"
        atom2.xyz = [0.0, 0.0, 0.74]
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

    def _create_sample_basisset(self, name="test-basis"):
        bs = BasisSet()
        bs.name = name
        cgto = ContractedGTO("s", 1)
        cgto[0] = PrimitiveGTO(1.2345, 0.6789)
        bs.append(cgto)
        return bs

    # --------------------------------------------------------------------------
    # 1. Molecule, atoms, electrons, AOs, MOs, dummy atoms, find_atom_index
    # --------------------------------------------------------------------------

    def test_molecule_property(self):
        param = PdfParamObject()
        self.assertIsInstance(param.molecule, bridge.AtomGroup)

        mol = bridge.AtomGroup()
        atom = bridge.Atom()
        atom.symbol = "H"
        atom.label = "1"
        mol.set_atom("1", atom)

        param.molecule = mol
        self.assertIs(param.molecule, mol)

        with self.assertRaises(AssertionError):
            param.molecule = "invalid_type"

    def test_num_of_atoms(self):
        param = PdfParamObject()
        self.assertIsNone(param.num_of_atoms)

        param.num_of_atoms = 10
        self.assertEqual(param.num_of_atoms, 10)
        self.assertIsInstance(param.num_of_atoms, int)

        param.num_of_atoms = "5"
        self.assertEqual(param.num_of_atoms, 5)

    def test_num_of_electrons_explicit(self):
        param = PdfParamObject()
        param.num_of_electrons = 42
        self.assertEqual(param.num_of_electrons, 42)

    def test_num_of_electrons_computed_from_molecule(self):
        param = self._create_h2_param()
        # H (atomic number 1) * 2 = 2.0
        self.assertEqual(param.num_of_electrons, 2.0)

    def test_num_of_electrons_with_dummy_atoms(self):
        param = self._create_h2_param()
        # Add dummy atom with charge 1.5
        dummy = bridge.Atom()
        dummy.symbol = "X"
        dummy.label = "dummy1"
        dummy.charge = 1.5
        param.molecule.set_atom("dummy1", dummy)

        # 2.0 (from 2 H atoms) + 1.5 (dummy atom charge) = 3.5
        self.assertEqual(param.num_of_electrons, 3.5)

    def test_get_charge_of_dummy_atoms_hierarchical(self):
        param = PdfParamObject()

        root = bridge.AtomGroup()
        x1 = bridge.Atom()
        x1.symbol = "X"
        x1.charge = 0.75
        root.set_atom("x1", x1)

        subgroup = bridge.AtomGroup()
        x2 = bridge.Atom()
        x2.symbol = "X"
        x2.charge = 1.25
        subgroup.set_atom("x2", x2)

        # Non-dummy atom should not contribute to dummy charge
        h = bridge.Atom()
        h.symbol = "H"
        h.charge = 0.5
        subgroup.set_atom("h1", h)

        root.set_group("sub1", subgroup)

        total_charge = param._get_charge_of_dummy_atoms(root)
        self.assertAlmostEqual(total_charge, 2.0)

    def test_num_of_AOs(self):
        param = PdfParamObject()
        self.assertIsNone(param.num_of_AOs)

        param.num_of_AOs = 15
        self.assertEqual(param.num_of_AOs, 15)
        self.assertIsInstance(param.num_of_AOs, int)

    def test_num_of_MOs(self):
        param = PdfParamObject()
        self.assertIsNone(param.num_of_MOs)

        param.num_of_MOs = 10
        self.assertEqual(param.num_of_MOs, 10)
        self.assertIsInstance(param.num_of_MOs, int)

        param.num_of_MOs = None
        self.assertEqual(param.num_of_MOs, 0)

    def test_find_atom_index(self):
        param = self._create_h2_param()

        # Retrieve atom 1
        a1 = param.molecule.get_atom("1")
        self.assertEqual(param.find_atom_index(a1), 0)

        # Retrieve atom 2
        a2 = param.molecule.get_atom("2")
        self.assertEqual(param.find_atom_index(a2), 1)

        # Atom not in molecule
        other_atom = bridge.Atom()
        other_atom.symbol = "O"
        other_atom.xyz = [0.0, 0.0, 1.5]
        self.assertEqual(param.find_atom_index(other_atom), -1)

    def test_find_atom_index_in_subgroup(self):
        param = PdfParamObject()
        mol = bridge.AtomGroup()

        sub = bridge.AtomGroup()
        atom_sub = bridge.Atom()
        atom_sub.symbol = "C"
        atom_sub.label = "C1"
        atom_sub.xyz = [1.0, 2.0, 3.0]
        sub.set_atom("C1", atom_sub)
        mol.set_group("sub", sub)

        param.molecule = mol
        self.assertEqual(param.find_atom_index(atom_sub), 0)

    # --------------------------------------------------------------------------
    # 2. BasisSet get/set (standard, j, xc, gridfree), names, atomlabels, basis2
    # --------------------------------------------------------------------------

    def test_basisset_object_roundtrip(self):
        param = PdfParamObject()
        bs_std = self._create_sample_basisset("std")
        bs_j = self._create_sample_basisset("j")
        bs_xc = self._create_sample_basisset("xc")
        bs_gf = self._create_sample_basisset("gf")

        param.set_basisset("H@1", bs_std)
        param.set_basisset_j("H@1", bs_j)
        param.set_basisset_xc("H@1", bs_xc)
        param.set_basisset_gridfree("H@1", bs_gf)

        self.assertEqual(param.get_basisset("H@1"), bs_std)
        self.assertEqual(param.get_basisset_j("H@1"), bs_j)
        self.assertEqual(param.get_basisset_xc("H@1"), bs_xc)
        self.assertEqual(param.get_basisset_gridfree("H@1"), bs_gf)

        # Unset atom label returns default empty BasisSet
        empty_bs = param.get_basisset("C@2")
        self.assertIsInstance(empty_bs, BasisSet)
        self.assertEqual(len(empty_bs), 0)

    def test_basisset_names(self):
        param = PdfParamObject()

        param.set_basisset_name("H@1", "sto-3g")
        param.set_basisset_j_name("H@1", "def2-universal-j")
        param.set_basisset_xc_name("H@1", "def2-universal-xc")
        param.set_basisset_gridfree_name("H@1", "gf-basis")

        self.assertEqual(param.get_basisset_name("H@1"), "sto-3g")
        self.assertEqual(param.get_basisset_j_name("H@1"), "def2-universal-j")
        self.assertEqual(param.get_basisset_xc_name("H@1"), "def2-universal-xc")
        self.assertEqual(param.get_basisset_gridfree_name("H@1"), "gf-basis")

        # Unset atom label returns empty string
        self.assertEqual(param.get_basisset_name("C@2"), "")
        self.assertEqual(param.get_basisset_j_name("C@2"), "")
        self.assertEqual(param.get_basisset_xc_name("C@2"), "")
        self.assertEqual(param.get_basisset_gridfree_name("C@2"), "")

    def test_basisset_atomlabels(self):
        param = PdfParamObject()
        bs = self._create_sample_basisset()

        param.set_basisset("H@1", bs)
        param.set_basisset_j("H@2", bs)
        param.set_basisset_xc("H@3", bs)
        param.set_basisset_gridfree("H@4", bs)

        self.assertIn("H@1", list(param.get_basisset_atomlabels()))
        self.assertIn("H@2", list(param.get_basisset_j_atomlabels()))
        self.assertIn("H@3", list(param.get_basisset_xc_atomlabels()))
        self.assertIn("H@4", list(param.get_basisset_gridfree_atomlabels()))

    def test_basis2_string_roundtrip(self):
        param = PdfParamObject()
        bs = self._create_sample_basisset("roundtrip-basis")

        param.set_basisset("H@1", bs)
        retrieved_bs = param.get_basisset("H@1")

        self.assertEqual(retrieved_bs.get_basis2(), bs.get_basis2())

    # --------------------------------------------------------------------------
    # 3. Total energy, TEs, gradient, gradient RMS, Mulliken population
    # --------------------------------------------------------------------------

    def test_total_energy_by_iteration_and_TEs(self):
        param = PdfParamObject()

        param.set_total_energy(1, -1.132)
        param.set_total_energy(2, -1.135)
        param.set_total_energy("3", -1.136)

        self.assertAlmostEqual(param.get_total_energy(1), -1.132)
        self.assertAlmostEqual(param.get_total_energy(2), -1.135)
        self.assertAlmostEqual(param.get_total_energy(3), -1.136)
        self.assertIsNone(param.get_total_energy(4))

        # TEs property
        tes = param.TEs
        self.assertEqual(tes[1], -1.132)
        self.assertEqual(tes[2], -1.135)
        self.assertEqual(tes[3], -1.136)

        param.TEs = {1: -2.0, 2: -2.5}
        self.assertAlmostEqual(param.get_total_energy(1), -2.0)
        self.assertAlmostEqual(param.get_total_energy(2), -2.5)

    def test_gradient_set_and_get(self):
        param = PdfParamObject()
        param.num_of_atoms = 2

        param.set_gradient(0, 0.01, 0.02, 0.03)
        param.set_gradient(1, -0.01, -0.02, -0.03)

        self.assertEqual(param.get_gradient(0), [0.01, 0.02, 0.03])
        self.assertEqual(param.get_gradient(1), [-0.01, -0.02, -0.03])

        with self.assertRaises(AssertionError):
            param.set_gradient(2, 0.0, 0.0, 0.0)

        with self.assertRaises(AssertionError):
            param.get_gradient(2)

    def test_gradient_rms_when_unset(self):
        param = PdfParamObject()
        param.num_of_atoms = 2
        self.assertIsNone(param.get_gradient_rms())

    @unittest.expectedFailure
    def test_gradient_rms_computed_expected_failure(self):
        """get_gradient_rms fails with NameError: name 'math' is not defined.

        Expected: sqrt((0.01^2 + 0.02^2 + 0.03^2 + 0.01^2 + 0.02^2 + 0.03^2) / (3 * 2))
        Actual: NameError: name 'math' is not defined in pdfparam_object.py:791
        """
        param = PdfParamObject()
        param.num_of_atoms = 2
        param.set_gradient(0, 0.01, 0.02, 0.03)
        param.set_gradient(1, -0.01, -0.02, -0.03)

        expected_rms = math.sqrt(
            (0.01**2 + 0.02**2 + 0.03**2 + (-0.01)**2 + (-0.02)**2 + (-0.03)**2) / (3 * 2)
        )
        self.assertAlmostEqual(param.get_gradient_rms(), expected_rms, places=7)

    def test_mulliken_atom_population(self):
        param = PdfParamObject()

        param.set_mulliken_atom_population(0, 0.987)
        param.set_mulliken_atom_population(1, 1.013)

        self.assertAlmostEqual(param.get_mulliken_atom_population(0), 0.987)
        self.assertAlmostEqual(param.get_mulliken_atom_population(1), 1.013)
        self.assertIsNone(param.get_mulliken_atom_population(2))

    # --------------------------------------------------------------------------
    # 4. Comparison methods: compare_info, compare_energy, compare_pop,
    #    compare_grad, _check
    # --------------------------------------------------------------------------

    def test_check_scalars(self):
        param = PdfParamObject()

        # float within / outside default error_threshold (1.0e-5)
        self.assertTrue(param._check(1.0, 1.0, "float equal"))
        self.assertTrue(param._check(1.0, 1.0 + error_threshold * 0.5, "float within threshold"))
        self.assertFalse(param._check(1.0, 1.0 + error_threshold * 1.5, "float outside threshold"))

        # custom threshold
        self.assertTrue(param._check(1.0, 1.05, "custom threshold pass", threshold=0.1))
        self.assertFalse(param._check(1.0, 1.05, "custom threshold fail", threshold=0.01))

        # int, bool, str equality
        self.assertTrue(param._check(10, 10, "int match"))
        self.assertFalse(param._check(10, 20, "int mismatch"))
        self.assertTrue(param._check(True, True, "bool match"))
        self.assertFalse(param._check(True, False, "bool mismatch"))
        self.assertTrue(param._check("abc", "abc", "str match"))
        self.assertFalse(param._check("abc", "xyz", "str mismatch"))

    def test_check_list_odd_length(self):
        param = PdfParamObject()
        # Odd-length list (3 elements): XOR logic happens to yield correct results
        # when all match (False ^ True ^ True ^ True -> True)
        self.assertTrue(param._check([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], "odd list match"))
        # 1 element differs: (False ^ False ^ True ^ True -> False)
        self.assertFalse(param._check([1.0, 2.0, 3.0], [1.0, 2.0, 9.0], "odd list 1 diff"))

    @unittest.expectedFailure
    def test_check_list_even_length_identical_expected_failure(self):
        """_check fails for identical even-length lists due to XOR accumulation bug.

        In _check, `answer = answer ^ self._check(...)` starts with answer = False.
        For a list of length 2 where both elements match:
          i=0: False ^ True = True
          i=1: True ^ True = False
        Expected: True
        Actual: False
        """
        param = PdfParamObject()
        self.assertTrue(param._check([1.0, 2.0], [1.0, 2.0], "even list identical"))

    @unittest.expectedFailure
    def test_check_list_odd_length_two_diffs_expected_failure(self):
        """_check returns True for list with 2 differences due to XOR accumulation bug.

        In _check with length 3 and 2 elements differing:
          i=0: False ^ True = True (matched)
          i=1: True ^ False = True (mismatched)
          i=2: True ^ False = True (mismatched)
        Expected: False (they do not match)
        Actual: True
        """
        param = PdfParamObject()
        self.assertFalse(param._check([1.0, 2.0, 3.0], [1.0, 9.0, 9.0], "odd list 2 diffs"))

    def test_compare_info(self):
        p1 = PdfParamObject()
        p1.num_of_AOs = 10
        p1.num_of_MOs = 5
        p1._set_scf_converged(True)

        p2 = PdfParamObject()
        p2.num_of_AOs = 10
        p2.num_of_MOs = 5
        p2._set_scf_converged(True)

        self.assertTrue(p1.compare_info(p2))

        # Different AOs
        p2.num_of_AOs = 12
        self.assertFalse(p1.compare_info(p2))
        p2.num_of_AOs = 10

        # Different MOs
        p2.num_of_MOs = 6
        self.assertFalse(p1.compare_info(p2))
        p2.num_of_MOs = 5

        # Different scf_converged
        p2._set_scf_converged(False)
        self.assertFalse(p1.compare_info(p2))

    def test_compare_energy(self):
        p1 = PdfParamObject()
        p1.iterations = 5
        p1.set_total_energy(5, -100.0)

        p2 = PdfParamObject()
        p2.iterations = 5
        p2.set_total_energy(5, -100.0)

        self.assertTrue(p1.compare_energy(p2))

        # Within error threshold
        p2.set_total_energy(5, -100.0 + error_threshold * 0.5)
        self.assertTrue(p1.compare_energy(p2))

        # Outside error threshold
        p2.set_total_energy(5, -100.0 + error_threshold * 2.0)
        self.assertFalse(p1.compare_energy(p2))

    @unittest.expectedFailure
    def test_compare_pop_expected_failure(self):
        """compare_pop fails with AttributeError: 'PdfParamObject' object has no attribute 'get_population'.

        compare_pop calls self.get_population('mulliken', itr1), but get_population
        is not defined on PdfParamObject.
        Expected: True for matching population data
        Actual: AttributeError: 'PdfParamObject' object has no attribute 'get_population'
        """
        p1 = PdfParamObject()
        p1.num_of_atoms = 2
        p1.iterations = 1

        p2 = PdfParamObject()
        p2.num_of_atoms = 2
        p2.iterations = 1

        self.assertTrue(p1.compare_pop(p2))

    @unittest.expectedFailure
    def test_compare_grad_expected_failure(self):
        """compare_grad fails because get_gradient_rms raises NameError (missing math import).

        Expected: True for matching gradients
        Actual: NameError: name 'math' is not defined in get_gradient_rms
        """
        p1 = PdfParamObject()
        p1.num_of_atoms = 1
        p1.set_gradient(0, 0.01, 0.02, 0.03)

        p2 = PdfParamObject()
        p2.num_of_atoms = 1
        p2.set_gradient(0, 0.01, 0.02, 0.03)

        self.assertTrue(p1.compare_grad(p2))


if __name__ == "__main__":
    unittest.main()
