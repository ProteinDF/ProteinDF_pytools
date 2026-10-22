# -*- coding: utf-8 -*-

import importlib.util
import io
import logging
import os
import sys
import tempfile
import unittest

import h5py
import numpy as np

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools import (
    BasisSet,
    ContractedGTO,
    Matrix,
    PdfParam,
    PdfParam_H5,
    PrimitiveGTO,
    SymmetricMatrix,
    Vector,
)


def _load_script(filename):
    path = os.path.join(os.path.dirname(__file__), "..", "scripts", filename)
    spec = importlib.util.spec_from_file_location(
        filename.replace("-", "_").replace(".py", ""), path
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _create_sample_roks_param(work_path=None, iterations=2):
    param = PdfParam()
    param.comment = "synthetic roks test system"
    param.num_of_atoms = 1
    param.num_of_AOs = 2
    param.num_of_MOs = 2
    param.method = "roks"
    param.iterations = iterations
    param.guess = "none"
    if work_path:
        param.work_path = work_path
    param._set_scf_converged(True)

    mol = bridge.AtomGroup()
    a1 = bridge.Atom(symbol="H", position=bridge.Position([0.0, 0.0, 0.0]), charge=1.0, name="H1")
    mol.set_atom("0", a1)
    param.molecule = mol

    bs = BasisSet()
    bs.name = "sto-3g-H"
    cgto = ContractedGTO("s", 1)
    cgto[0] = PrimitiveGTO(3.42525, 0.15433)
    bs.append(cgto)
    param.set_basisset("H@0", bs)

    for itr in range(1, iterations + 1):
        param.set_total_energy(itr, -0.5 * itr)

    param._data.setdefault("control", {})["file_base_name"] = {
        "eigenvalues": "eigenvalues.%s.vtr",
        "occupation_vtr": "occupation.%s.vtr",
        "Spq_matrix": "Spq.mat",
        "Hpq_matrix": "Hpq.mat",
        "Hpq2_matrix": "Hpq2.mat",
        "Fpq_matrix": "Fpq.%s.mat",
        "C_matrix": "C.%s.mat",
        "Ppq_matrix": "Ppq.%s.mat",
    }
    return param


def _create_synthetic_roks_files(param):
    """Write synthetic vector and matrix files required for ROKS archive."""
    # Occupations (roks_closed and roks_open)
    v_occ_closed = Vector(2)
    v_occ_closed[0] = 2.0
    v_occ_closed[1] = 0.0
    v_occ_closed.save(param.get_occ_path("roks_closed"))

    v_occ_open = Vector(2)
    v_occ_open[0] = 0.0
    v_occ_open[1] = 1.0
    v_occ_open.save(param.get_occ_path("roks_open"))

    for itr in range(1, param.iterations + 1):
        # Orbitals: eigenvalues, C, F (single set for 'roks')
        v_eig = Vector(2)
        v_eig[0] = -0.6 * itr
        v_eig[1] = 0.1 * itr
        v_eig.save(param.get_energy_level_path("roks", itr))

        mat_c = Matrix(2, 2)
        mat_c.set(0, 0, 1.0)
        mat_c.set(1, 1, 1.0)
        mat_c.save(param.get_c_mat_path("roks", itr))

        mat_f = SymmetricMatrix(2)
        mat_f.set(0, 0, 0.4)
        mat_f.save(param.get_f_mat_path("roks", itr))

        # Density: P (two sets: roks_closed and roks_open)
        mat_p_closed = SymmetricMatrix(2)
        mat_p_closed.set(0, 0, 1.0)
        mat_p_closed.save(param.get_density_matrix_path("roks_closed", itr))

        mat_p_open = SymmetricMatrix(2)
        mat_p_open.set(1, 1, 0.5)
        mat_p_open.save(param.get_density_matrix_path("roks_open", itr))

    s_mat = SymmetricMatrix(2)
    s_mat.set(0, 0, 1.0)
    s_mat.set(1, 1, 1.0)
    s_mat.save(param.get_s_mat_path())
    s_mat.save(param.get_h_mat_path())
    s_mat.save(param.get_h2_mat_path())


class TestRoksRunTypes(unittest.TestCase):
    def test_run_types_methods(self):
        param = PdfParam()

        param.method = "rks"
        self.assertEqual(param.run_types(), ["rks"])
        self.assertEqual(param.orbital_run_types(), ["rks"])

        param.method = "uks"
        self.assertEqual(param.run_types(), ["uks_alpha", "uks_beta"])
        self.assertEqual(param.orbital_run_types(), ["uks_alpha", "uks_beta"])

        param.method = "roks"
        self.assertEqual(param.run_types(), ["roks_closed", "roks_open"])
        self.assertEqual(param.orbital_run_types(), ["roks"])


class TestRoksArchiveH5(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_roks_archive_and_readback(self):
        """Synthetic ROKS data is archived by pdf-archive-h5.py and read back correctly."""
        work_path = os.path.join(self.tmpdir.name, "fl_Work")
        os.makedirs(work_path, exist_ok=True)
        param = _create_sample_roks_param(work_path=work_path, iterations=2)
        _create_synthetic_roks_files(param)

        mpac_path = os.path.join(self.tmpdir.name, "roks.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "pdfresults.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]
        mod.main()

        self.assertTrue(os.path.exists(h5_path))

        # Check HDF5 structure directly
        with h5py.File(h5_path, "r") as h5_raw:
            # occ: roks_closed, roks_open
            self.assertIn("occ/roks_closed", h5_raw)
            self.assertIn("occ/roks_open", h5_raw)
            self.assertNotIn("occ/roks", h5_raw)

            # energy_level: roks_1, roks_2
            self.assertIn("energy_level/roks_1", h5_raw)
            self.assertIn("energy_level/roks_2", h5_raw)
            self.assertNotIn("energy_level/roks_closed_1", h5_raw)

            # C: roks_1, roks_2
            self.assertIn("C/roks_1", h5_raw)
            self.assertIn("C/roks_2", h5_raw)

            # F: roks_1, roks_2
            self.assertIn("F/roks_1", h5_raw)
            self.assertIn("F/roks_2", h5_raw)

            # P: roks_closed_1, roks_open_1, roks_closed_2, roks_open_2
            self.assertIn("P/roks_closed_1", h5_raw)
            self.assertIn("P/roks_open_1", h5_raw)
            self.assertIn("P/roks_closed_2", h5_raw)
            self.assertIn("P/roks_open_2", h5_raw)
            self.assertNotIn("P/roks_1", h5_raw)

        # Read back with PdfParam_H5
        h5 = PdfParam_H5()
        h5.open(h5_path)
        self.assertEqual(h5.iterations, 2)

        # Occupations
        occ_closed = h5.get_occ_vector("roks_closed")
        self.assertAlmostEqual(occ_closed[0], 2.0)
        self.assertAlmostEqual(occ_closed[1], 0.0)

        occ_open = h5.get_occ_vector("roks_open")
        self.assertAlmostEqual(occ_open[0], 0.0)
        self.assertAlmostEqual(occ_open[1], 1.0)

        # HOMO level for roks should resolve to 1 (from open orbital)
        self.assertEqual(h5.get_HOMO_level("roks"), 1)

        # Eigenvalues
        e1 = h5.get_energy_level("roks", 1)
        self.assertAlmostEqual(e1[0], -0.6)
        self.assertAlmostEqual(e1[1], 0.1)

        e2 = h5.get_energy_level("roks", 2)
        self.assertAlmostEqual(e2[0], -1.2)
        self.assertAlmostEqual(e2[1], 0.2)

        # Matrices
        for itr in range(1, h5.iterations + 1):
            c = h5.get_c_matrix("roks", itr)
            self.assertIsInstance(c, Matrix)
            self.assertEqual(c.rows, 2)
            self.assertEqual(c.cols, 2)

            f = h5.get_f_matrix("roks", itr)
            self.assertIsInstance(f, SymmetricMatrix)
            self.assertEqual(f.rows, 2)
            self.assertEqual(f.cols, 2)

            p_closed = h5.get_density_matrix("roks_closed", itr)
            self.assertIsInstance(p_closed, SymmetricMatrix)
            self.assertEqual(p_closed.rows, 2)
            self.assertEqual(p_closed.cols, 2)

            p_open = h5.get_density_matrix("roks_open", itr)
            self.assertIsInstance(p_open, SymmetricMatrix)
            self.assertEqual(p_open.rows, 2)
            self.assertEqual(p_open.cols, 2)

        h5._close()

    def test_roks_missing_files_warning(self):
        """Warning message contains run type name and file path when file fails to load."""
        work_path = os.path.join(self.tmpdir.name, "fl_Work_empty")
        os.makedirs(work_path, exist_ok=True)
        param = _create_sample_roks_param(work_path=work_path, iterations=1)
        # S, H, H2 matrices are saved so setup doesn't fail early
        s_mat = SymmetricMatrix(2)
        s_mat.save(param.get_s_mat_path())
        s_mat.save(param.get_h_mat_path())
        s_mat.save(param.get_h2_mat_path())

        mpac_path = os.path.join(self.tmpdir.name, "roks_empty.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "pdfresults_empty.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]

        with self.assertLogs("pdf_archive_h5", level="WARNING") as cm:
            mod.main()

        output_log = "\n".join(cm.output)
        self.assertIn("roks_closed", output_log)
        self.assertIn("roks_open", output_log)
        self.assertIn("roks", output_log)
        self.assertIn("occupation.roks_closed.vtr", output_log)
        self.assertIn("eigenvalues.roks1.vtr", output_log)


if __name__ == "__main__":
    unittest.main()
