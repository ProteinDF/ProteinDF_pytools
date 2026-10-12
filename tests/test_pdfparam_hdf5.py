#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

try:
    import h5py
    import numpy as np

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

import proteindf_bridge as bridge
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.matrix import Matrix, SymmetricMatrix
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5
from proteindf_tools.vector import Vector


@unittest.skipUnless(HAVE_H5PY, "h5py or numpy not available")
class TestPdfParamH5(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.h5_path = os.path.join(self.tmpdir.name, "test.h5")

    def tearDown(self):
        self.tmpdir.cleanup()

    # -------------------------------------------------------------------------
    # 1. save_model -> load_model roundtrip
    # -------------------------------------------------------------------------
    def test_save_load_model_control_parameters(self):
        """Verify that control parameters are preserved across save_model and load_model."""
        p = PdfParam_H5()
        p.open(self.h5_path)

        p.step_control = "scf"
        p.comment = "test comment"
        p.num_of_atoms = 10
        p.num_of_AOs = 20
        p.num_of_MOs = 20
        p.method = "rks"
        p.cut_value = 1.0e-10
        p.CDAM_tau = 1.0e-5
        p.CD_epsilon = 1.0e-4
        p.guess = "core"
        p.xc_functional = "b3lyp"
        p.iterations = 0
        p.max_iterations = 50
        p.orbital_independence_threshold = 1.0e-6
        p.orbital_independence_threshold_canonical = 1.0e-7
        p.orbital_independence_threshold_lowdin = 1.0e-8
        p.scf_acceleration = "damping"

        p.save_model()
        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)

        self.assertEqual(p2.step_control, "scf")
        self.assertEqual(p2.comment, "test comment")
        self.assertEqual(p2.num_of_atoms, 10)
        self.assertEqual(p2.num_of_AOs, 20)
        self.assertEqual(p2.num_of_MOs, 20)
        self.assertEqual(p2.method, "rks")
        self.assertAlmostEqual(p2.cut_value, 1.0e-10)
        self.assertAlmostEqual(p2.CDAM_tau, 1.0e-5)
        self.assertAlmostEqual(p2.CD_epsilon, 1.0e-4)
        self.assertEqual(p2.guess, "core")
        self.assertEqual(p2.xc_functional, "b3lyp")
        self.assertEqual(p2.iterations, 0)
        self.assertEqual(p2.max_iterations, 50)
        self.assertAlmostEqual(p2.orbital_independence_threshold, 1.0e-6)
        self.assertAlmostEqual(p2.orbital_independence_threshold_canonical, 1.0e-7)
        self.assertAlmostEqual(p2.orbital_independence_threshold_lowdin, 1.0e-8)
        self.assertEqual(p2.scf_acceleration, "damping")
        p2._close()

    def test_save_load_model_molecule(self):
        """Verify that molecule structure is preserved across save_model and load_model."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.iterations = 0
        p.num_of_atoms = 2
        p.num_of_AOs = 1
        p.num_of_MOs = 1

        mol = bridge.AtomGroup()
        atom_o = bridge.Atom()
        atom_o.symbol = "O"
        atom_o.name = "1"
        atom_o.xyz = bridge.Position(0.0, 0.5, 1.0)
        atom_o.charge = 8.0
        mol.set_atom("1", atom_o)

        atom_h = bridge.Atom()
        atom_h.symbol = "H"
        atom_h.name = "2"
        atom_h.xyz = bridge.Position(0.0, -0.5, 2.0)
        atom_h.charge = 1.0
        mol.set_atom("2", atom_h)
        p.molecule = mol

        p.save_model()
        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)
        self.assertEqual(p2.molecule.get_number_of_atoms(), 2)

        atom1 = p2.molecule.get_atom("1")
        self.assertIsNotNone(atom1)
        self.assertEqual(atom1.symbol, "O")
        self.assertEqual(atom1.name, "1")
        self.assertAlmostEqual(atom1.xyz.x, 0.0)
        self.assertAlmostEqual(atom1.xyz.y, 0.5)
        self.assertAlmostEqual(atom1.xyz.z, 1.0)
        self.assertAlmostEqual(atom1.charge, 8.0)

        atom2 = p2.molecule.get_atom("2")
        self.assertIsNotNone(atom2)
        self.assertEqual(atom2.symbol, "H")
        self.assertEqual(atom2.name, "2")
        self.assertAlmostEqual(atom2.xyz.x, 0.0)
        self.assertAlmostEqual(atom2.xyz.y, -0.5)
        self.assertAlmostEqual(atom2.xyz.z, 2.0)
        self.assertAlmostEqual(atom2.charge, 1.0)
        p2._close()

    def test_save_load_model_basisset(self):
        """Verify that basis sets are preserved across save_model and load_model."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.iterations = 0
        p.num_of_atoms = 2
        p.num_of_AOs = 1
        p.num_of_MOs = 1

        bs_o = BasisSet()
        bs_o.name = "sto-3g-O"
        cgto_o_s = ContractedGTO("s", 2)
        cgto_o_s.scale_factor = 1.1
        cgto_o_s[0] = PrimitiveGTO(130.70932, 0.154328)
        cgto_o_s[1] = PrimitiveGTO(23.80886, 0.535328)
        bs_o.append(cgto_o_s)

        cgto_o_p = ContractedGTO("p", 2)
        cgto_o_p.scale_factor = 1.2
        cgto_o_p[0] = PrimitiveGTO(5.03315, 0.444635)
        cgto_o_p[1] = PrimitiveGTO(1.16959, 0.612101)
        bs_o.append(cgto_o_p)
        p.set_basisset("O@1", bs_o)

        bs_h = BasisSet()
        bs_h.name = "sto-3g-H"
        cgto_h_s = ContractedGTO("s", 2)
        cgto_h_s.scale_factor = 0.95
        cgto_h_s[0] = PrimitiveGTO(3.42525, 0.15433)
        cgto_h_s[1] = PrimitiveGTO(0.62391, 0.53533)
        bs_h.append(cgto_h_s)
        p.set_basisset("H@2", bs_h)

        p.save_model()
        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)
        labels = p2.get_basisset_atomlabels()
        self.assertIn("O@1", labels)
        self.assertIn("H@2", labels)

        bs_o_loaded = p2.get_basisset("O@1")
        self.assertEqual(bs_o_loaded.name, "sto-3g-O")
        self.assertEqual(len(bs_o_loaded), 2)

        cgto0 = bs_o_loaded[0]
        self.assertEqual(cgto0.shell_type, "s")
        self.assertAlmostEqual(cgto0.scale_factor, 1.1)
        self.assertEqual(len(cgto0), 2)
        self.assertAlmostEqual(cgto0[0].exp, 130.70932)
        self.assertAlmostEqual(cgto0[0].coef, 0.154328)
        self.assertAlmostEqual(cgto0[1].exp, 23.80886)
        self.assertAlmostEqual(cgto0[1].coef, 0.535328)

        cgto1 = bs_o_loaded[1]
        self.assertEqual(cgto1.shell_type, "p")
        self.assertAlmostEqual(cgto1.scale_factor, 1.2)
        self.assertEqual(len(cgto1), 2)
        self.assertAlmostEqual(cgto1[0].exp, 5.03315)
        self.assertAlmostEqual(cgto1[0].coef, 0.444635)
        self.assertAlmostEqual(cgto1[1].exp, 1.16959)
        self.assertAlmostEqual(cgto1[1].coef, 0.612101)

        bs_h_loaded = p2.get_basisset("H@2")
        self.assertEqual(bs_h_loaded.name, "sto-3g-H")
        self.assertEqual(len(bs_h_loaded), 1)

        cgto_h = bs_h_loaded[0]
        self.assertEqual(cgto_h.shell_type, "s")
        self.assertAlmostEqual(cgto_h.scale_factor, 0.95)
        self.assertEqual(len(cgto_h), 2)
        self.assertAlmostEqual(cgto_h[0].exp, 3.42525)
        self.assertAlmostEqual(cgto_h[0].coef, 0.15433)
        self.assertAlmostEqual(cgto_h[1].exp, 0.62391)
        self.assertAlmostEqual(cgto_h[1].coef, 0.53533)
        p2._close()

    # -------------------------------------------------------------------------
    # 2. occ and energy_level roundtrip
    # -------------------------------------------------------------------------
    def test_occ_vector_and_homo_level_roundtrip(self):
        """Verify save_occ, get_occ_vector, and get_HOMO_level roundtrip."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.iterations = 0
        p.num_of_atoms = 1
        p.num_of_AOs = 5
        p.num_of_MOs = 5
        p.save_model()

        occ = Vector([2.0, 2.0, 1.0, 0.0, 0.0])
        p.save_occ("rks", occ)
        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)
        occ_loaded = p2.get_occ_vector("rks")
        self.assertIsInstance(occ_loaded, Vector)
        self.assertEqual(len(occ_loaded), 5)
        for i in range(len(occ)):
            self.assertAlmostEqual(occ[i], occ_loaded[i])

        # HOMO level is 0-indexed; last orbital with occ > 0 is index 2
        self.assertEqual(p2.get_HOMO_level("rks"), 2)
        p2._close()

    def test_energy_level_roundtrip(self):
        """Verify save_energy_level and get_energy_level roundtrip."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.iterations = 0
        p.num_of_atoms = 1
        p.num_of_AOs = 5
        p.num_of_MOs = 5
        p.save_model()

        ene = Vector([-10.5, -1.2, -0.5, 0.2, 0.8])
        p.save_energy_level("rks", 1, ene)
        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)
        ene_loaded = p2.get_energy_level("rks", 1)
        self.assertIsInstance(ene_loaded, Vector)
        self.assertEqual(len(ene_loaded), 5)
        for i in range(len(ene)):
            self.assertAlmostEqual(ene[i], ene_loaded[i])
        p2._close()

    # -------------------------------------------------------------------------
    # 3. matrix roundtrip
    # -------------------------------------------------------------------------
    def test_symmetric_matrices_roundtrip(self):
        """Verify save and get for symmetric matrices (S, h, h2, P, F)."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.iterations = 0
        p.num_of_atoms = 1
        p.num_of_AOs = 2
        p.num_of_MOs = 2
        p.save_model()

        arr_s = np.array([[1.0, 0.2], [0.2, 1.0]], dtype=float)
        mat_s = SymmetricMatrix(arr_s)
        p.save_s_matrix(mat_s)

        arr_h = np.array([[2.0, 0.5], [0.5, 2.0]], dtype=float)
        mat_h = SymmetricMatrix(arr_h)
        p.save_h_matrix(mat_h)

        arr_h2 = np.array([[3.0, 0.1], [0.1, 3.0]], dtype=float)
        mat_h2 = SymmetricMatrix(arr_h2)
        p.save_h2_matrix(mat_h2)

        arr_p = np.array([[0.4, 0.3], [0.3, 0.4]], dtype=float)
        mat_p = SymmetricMatrix(arr_p)
        p.save_density_matrix("rks", 1, mat_p)

        arr_f = np.array([[-1.5, -0.2], [-0.2, -1.5]], dtype=float)
        mat_f = SymmetricMatrix(arr_f)
        p.save_f_matrix("rks", 1, mat_f)

        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)

        for name, original, loaded in [
            ("s", mat_s, p2.get_s_matrix()),
            ("h", mat_h, p2.get_h_matrix()),
            ("h2", mat_h2, p2.get_h2_matrix()),
            ("density", mat_p, p2.get_density_matrix("rks", 1)),
            ("f", mat_f, p2.get_f_matrix("rks", 1)),
        ]:
            self.assertIsInstance(loaded, SymmetricMatrix, msg=f"{name} should be SymmetricMatrix")
            self.assertEqual(loaded.rows, original.rows)
            self.assertEqual(loaded.cols, original.cols)
            self.assertEqual(loaded.type, "SY")
            self.assertTrue(np.allclose(original.data, loaded.data), msg=f"{name} data mismatch")
        p2._close()

    def test_general_matrix_roundtrip(self):
        """Verify save and get for general matrix (C matrix)."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.iterations = 0
        p.num_of_atoms = 1
        p.num_of_AOs = 2
        p.num_of_MOs = 2
        p.save_model()

        arr_c = np.array([[0.8, 0.6], [-0.6, 0.8]], dtype=float)
        mat_c = Matrix(arr_c)
        p.save_c_matrix("rks", 1, mat_c)
        p._close()

        p2 = PdfParam_H5()
        p2.open(self.h5_path)
        loaded_c = p2.get_c_matrix("rks", 1)

        self.assertIsInstance(loaded_c, Matrix)
        self.assertNotIsInstance(loaded_c, SymmetricMatrix)
        self.assertEqual(loaded_c.rows, mat_c.rows)
        self.assertEqual(loaded_c.cols, mat_c.cols)
        self.assertEqual(loaded_c.type, "GE")
        self.assertTrue(np.allclose(mat_c.data, loaded_c.data))
        p2._close()

    # -------------------------------------------------------------------------
    # 4. Nonexistent files and unsaved items behavior
    # -------------------------------------------------------------------------
    def test_open_nonexistent_filepath(self):
        """Verify open() behavior for nonexistent paths."""
        p = PdfParam_H5()

        # Empty path raises ValueError
        with self.assertRaises(ValueError):
            p.open("")

        # Path in non-existent directory raises OSError
        nonexistent_dir_path = os.path.join(self.tmpdir.name, "no_such_dir", "file.h5")
        with self.assertRaises(OSError):
            p.open(nonexistent_dir_path)

        # Non-existent file in existing directory is created and opened in mode 'a'
        new_file_path = os.path.join(self.tmpdir.name, "new_file.h5")
        p.open(new_file_path)
        self.assertTrue(os.path.exists(new_file_path))
        self.assertEqual(p.filepath, new_file_path)
        p._close()

    def test_read_unsaved_items(self):
        """Verify behavior when reading unsaved items from an open HDF5 file."""
        p = PdfParam_H5()
        p.open(self.h5_path)

        # Unsaved vectors return None
        self.assertIsNone(p.get_occ_vector("rks"))
        self.assertIsNone(p.get_energy_level("rks", 1))
        self.assertIsNone(p.get_TEs())

        # get_HOMO_level raises TypeError when occ vector is not saved
        with self.assertRaises(TypeError):
            p.get_HOMO_level("rks")

        # Unsaved matrices raise AttributeError because dataset is None
        with self.assertRaises(AttributeError):
            p.get_s_matrix()
        with self.assertRaises(AttributeError):
            p.get_h_matrix()
        with self.assertRaises(AttributeError):
            p.get_c_matrix("rks", 1)
        with self.assertRaises(AttributeError):
            p.get_density_matrix("rks", 1)
        with self.assertRaises(AttributeError):
            p.get_f_matrix("rks", 1)
        p._close()

    def test_open_file_without_saved_model(self):
        """Verify that opening an existing HDF5 file without model attributes raises TypeError."""
        with h5py.File(self.h5_path, "w") as f:
            f.attrs["unrelated"] = 1

        p = PdfParam_H5()
        # _load_model attempts int(None) for num_of_atoms, raising TypeError
        with self.assertRaises(TypeError):
            p.open(self.h5_path)

    def test_open_file_with_iterations_but_unsaved_tes(self):
        """Verify that opening a file with iterations > 0 but unsaved TEs raises TypeError."""
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.num_of_atoms = 1
        p.num_of_AOs = 1
        p.num_of_MOs = 1
        p.iterations = 2
        p.save_model()
        p._close()

        p2 = PdfParam_H5()
        # open() attempts TEs[itr] where TEs is None, raising TypeError
        with self.assertRaises(TypeError):
            p2.open(self.h5_path)


if __name__ == "__main__":
    unittest.main()
