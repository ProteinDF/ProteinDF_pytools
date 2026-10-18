#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import importlib
import importlib.util
import io
import math
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam

try:
    import pandas
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    import sklearn
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _create_sample_param():
    """Create minimal synthetic PdfParam representing water (H2O: O, H, H)."""
    param = PdfParam()
    param.comment = "synthetic water system"
    param.num_of_atoms = 3
    param.num_of_AOs = 5
    param.num_of_MOs = 5
    param.method = "rks"
    param.iterations = 1
    param.guess = "none"
    param._set_scf_converged(True)

    mol = bridge.AtomGroup()
    a_o = bridge.Atom(
        symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="O1"
    )
    a_h1 = bridge.Atom(
        symbol="H", position=bridge.Position([0.0, 1.0, 0.0]), charge=1.0, name="H1"
    )
    a_h2 = bridge.Atom(
        symbol="H", position=bridge.Position([1.0, 0.0, 0.0]), charge=1.0, name="H2"
    )
    mol.set_atom("0", a_o)
    mol.set_atom("1", a_h1)
    mol.set_atom("2", a_h2)
    param.molecule = mol

    # BasisSet for O: 1 s shell
    bs_o = BasisSet()
    bs_o.name = "sto-3g-O"
    cgto_o = ContractedGTO("s", 1)
    cgto_o[0] = PrimitiveGTO(71.6168, 0.15433)
    bs_o.append(cgto_o)

    # BasisSet for H: 1 s shell
    bs_h = BasisSet()
    bs_h.name = "sto-3g-H"
    cgto_h = ContractedGTO("s", 1)
    cgto_h[0] = PrimitiveGTO(3.42525, 0.15433)
    bs_h.append(cgto_h)

    param.set_basisset("O", bs_o)
    param.set_basisset("H", bs_h)
    param.set_total_energy(1, -75.0)

    param._data.setdefault("control", {})["file_base_name"] = {
        "eigenvalues": "eigenvalues.%s.vtr",
        "occupation_vtr": "occ.%s.vtr",
        "Spq_matrix": "s.mat",
        "Hpq_matrix": "h.mat",
        "Hpq2_matrix": "h2.mat",
        "Fpq_matrix": "f.%s.mat",
        "C_matrix": "c.%s.mat",
        "Ppq_matrix": "p.%s.mat",
    }
    return param


def _create_sample_system(tmpdir_path):
    """Create a consistent synthetic regression system for 3 atoms (O, H, H).

    Known true charges: q_true = [0.25, -0.125, -0.125] with total charge Q = 0.0.
    Returns: dict of paths (mpac, db, esp, mat, vtr, charge_csv).
    """
    param = _create_sample_param()
    mpac_path = os.path.join(tmpdir_path, "water.mpac")
    pdf.save_pdfparam(param, mpac_path)

    db_path = os.path.join(tmpdir_path, "water.db")
    PdfArchive_Sqlite3(db_path, pdfparam=param)

    # 4x4 design matrix for 3 real atoms + 1 Lagrange multiplier constraint
    # X = [[4, 1, 1, 1],
    #      [1, 4, 1, 1],
    #      [1, 1, 4, 1],
    #      [1, 1, 1, 0]]
    mat = pdf.Matrix(4, 4)
    for i in range(3):
        for j in range(3):
            mat.set(i, j, 4.0 if i == j else 1.0)
        mat.set(i, 3, 1.0)
        mat.set(3, i, 1.0)
    mat.set(3, 3, 0.0)
    mat_path = os.path.join(tmpdir_path, "esp_design.mat")
    mat.save(mat_path)

    # Target vector corresponding to q_true = [0.25, -0.125, -0.125, 0.0]
    # y = X * q_ext = [0.75, -0.375, -0.375, 0.0]
    vtr = pdf.Vector([0.75, -0.375, -0.375, 0.0])
    vtr_path = os.path.join(tmpdir_path, "esp_target.vtr")
    vtr.save(vtr_path)

    # ESP data at grid point (0.0, 0.0, 2.0) matching exact ESP from q_true
    atoms, _ = pdf.PopUtils._get_atomlist(param.molecule)
    q_true = [0.25, -0.125, -0.125]
    for i in range(3):
        atoms[i].charge = q_true[i]

    grid_pt = [0.0, 0.0, 2.0]
    pos = bridge.Position(grid_pt[0], grid_pt[1], grid_pt[2])
    pos *= 1.0 / 0.5291772108  # Angstrom to a.u.
    esp_utils = pdf.PopEspUtils()
    exact_esp = esp_utils._calc_esp(pos, atoms)
    esp_data = {"grids": [grid_pt], "ESP": [exact_esp]}
    esp_path = os.path.join(tmpdir_path, "grid-esp.mpac")
    bridge.save_msgpack(esp_data, esp_path)

    # Charges CSV for pdf-pop-rrms.py
    charge_csv_path = os.path.join(tmpdir_path, "charges.csv")
    with open(charge_csv_path, "w") as f:
        for q in q_true:
            f.write(f"{q}\n")

    return {
        "mpac_path": mpac_path,
        "db_path": db_path,
        "mat_path": mat_path,
        "vtr_path": vtr_path,
        "esp_path": esp_path,
        "charge_csv_path": charge_csv_path,
    }


class TestPdfPopRidge(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self._orig_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.sys_data = _create_sample_system(self.tmpdir.name)
        self.mod = _load_script("pdf-pop-ridge.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_ridge_param_success(self):
        """Ridge regression calculates partial charges from param file."""
        out_csv = os.path.join(self.tmpdir.name, "out.csv")
        out_vtr = os.path.join(self.tmpdir.name, "model.vtr")
        sys.argv = [
            "pdf-pop-ridge.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--alpha",
            "0.001",
            "--csv",
            out_csv,
            "--model_vector_path",
            out_vtr,
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("total charges:  0.000", output)
        self.assertIn("RRMS=", output)
        self.assertIn("O1", output)
        self.assertIn("H1", output)
        self.assertIn("H2", output)

        # Check CSV output
        self.assertTrue(os.path.exists(out_csv))
        with open(out_csv) as f:
            lines = [line.strip() for line in f if line.strip()]
        self.assertEqual(len(lines), 3)

        # Check model vector output
        self.assertTrue(os.path.exists(out_vtr))
        saved_vtr = pdf.Vector()
        saved_vtr.load(out_vtr)
        self.assertEqual(len(saved_vtr), 3)
        self.assertAlmostEqual(saved_vtr[0], 0.25, delta=0.01)
        self.assertAlmostEqual(saved_vtr[1], -0.125, delta=0.01)
        self.assertAlmostEqual(saved_vtr[2], -0.125, delta=0.01)

    def test_ridge_alpha_regularization(self):
        """Small alpha approximates OLS solution and large alpha shrinks coefficients."""
        mat = pdf.Matrix()
        mat.load(self.sys_data["mat_path"])
        vtr = pdf.Vector()
        vtr.load(self.sys_data["vtr_path"])

        # Weak regularization -> close to exact least squares solution
        ridge_weak = self.mod.Ridge(alpha=1e-8)
        ridge_weak.fit_atomic_charges(mat, vtr)
        coef_weak = list(ridge_weak.coef[:-1])
        self.assertAlmostEqual(coef_weak[0], 0.25, delta=1e-4)
        self.assertAlmostEqual(coef_weak[1], -0.125, delta=1e-4)
        self.assertAlmostEqual(coef_weak[2], -0.125, delta=1e-4)

        # Strong regularization -> coefficients shrink towards zero
        ridge_strong = self.mod.Ridge(alpha=100.0)
        ridge_strong.fit_atomic_charges(mat, vtr)
        coef_strong = list(ridge_strong.coef[:-1])
        norm_weak = sum(c * c for c in coef_weak)
        norm_strong = sum(c * c for c in coef_strong)
        self.assertLess(norm_strong, norm_weak)
        self.assertLess(norm_strong, 0.01)

    def test_ridge_db_success(self):
        """Ridge regression calculates partial charges from database file."""
        sys.argv = [
            "pdf-pop-ridge.py",
            "-d",
            self.sys_data["db_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--alpha",
            "0.001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertIn("total charges:  0.000", out.getvalue())

    def test_ridge_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-ridge.py",
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_ridge_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-ridge.py",
            "-p",
            self.sys_data["mpac_path"],
            "-d",
            self.sys_data["db_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_ridge_nonexistent_param_fails(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = [
            "pdf-pop-ridge.py",
            "-p",
            os.path.join(self.tmpdir.name, "missing.mpac"),
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()

    def test_ridge_nonexistent_design_matrix_fails(self):
        """Nonexistent design matrix file raises ValueError."""
        sys.argv = [
            "pdf-pop-ridge.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            os.path.join(self.tmpdir.name, "missing.mat"),
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(ValueError):
            self.mod.main()

    def test_ridge_dimension_mismatch_fails(self):
        """Mismatch between atom count and target vector size raises AssertionError."""
        bad_vtr = pdf.Vector([0.5, 0.0])  # Size 2 instead of 4
        bad_vtr_path = os.path.join(self.tmpdir.name, "bad.vtr")
        bad_vtr.save(bad_vtr_path)

        sys.argv = [
            "pdf-pop-ridge.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            bad_vtr_path,
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises((AssertionError, IndexError)):
            self.mod.main()


@unittest.skipUnless(HAS_PANDAS, "pandas is required")
class TestPdfPopLasso(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self._orig_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.sys_data = _create_sample_system(self.tmpdir.name)
        self.mod = _load_script("pdf-pop-lasso.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_lasso_fit_th_regularization(self):
        """Lasso fit_TH converges to exact solution with small alpha and shrinks with large alpha.
        """
        mat = pdf.Matrix()
        mat.load(self.sys_data["mat_path"])
        vtr = pdf.Vector()
        vtr.load(self.sys_data["vtr_path"])

        # Small alpha -> close to known solution
        lasso_weak = self.mod.Lasso(alpha=1e-6, max_iter=2000, tol=1e-5)
        lasso_weak.fit_TH(mat, vtr)
        self.assertTrue(lasso_weak.is_converged)
        coef_weak = list(lasso_weak.coef[:-1])
        self.assertAlmostEqual(coef_weak[0], 0.25, delta=1e-3)
        self.assertAlmostEqual(coef_weak[1], -0.125, delta=1e-3)
        self.assertAlmostEqual(coef_weak[2], -0.125, delta=1e-3)

        # Strong alpha -> coefficients shrink towards zero
        lasso_strong = self.mod.Lasso(alpha=5.0, max_iter=2000, tol=1e-5)
        lasso_strong.fit_TH(mat, vtr)
        coef_strong = list(lasso_strong.coef[:-1])
        norm_weak = sum(c * c for c in coef_weak)
        norm_strong = sum(c * c for c in coef_strong)
        self.assertLess(norm_strong, norm_weak)

    def test_lasso_param_success(self):
        """Lasso regression calculates partial charges from param file."""
        out_csv = os.path.join(self.tmpdir.name, "out.csv")
        sys.argv = [
            "pdf-pop-lasso.py",
            "-p",
            self.sys_data["mpac_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "-m",
            "1000",
            "-t",
            "0.0001",
            "--alpha",
            "0.001",
            "--csv",
            out_csv,
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("charge=", output)
        self.assertIn("0.000", output)
        self.assertIn("RRMS=", output)
        self.assertTrue(os.path.exists(out_csv))
        with open(out_csv) as f:
            lines = [line.strip() for line in f if line.strip()]
        self.assertEqual(len(lines), 3)

    def test_lasso_param_default_m_fails(self):
        """Lasso regression runs with default max_iterations."""
        sys.argv = [
            "pdf-pop-lasso.py",
            "-p",
            self.sys_data["mpac_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--alpha",
            "0.001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertIn("charge=0.000", out.getvalue())

    @unittest.skipUnless(HAS_SKLEARN, "scikit-learn is required")
    def test_lasso_param_scikit_learn(self):
        """Lasso regression with --use_scikit_learn calculates charges."""
        sys.argv = [
            "pdf-pop-lasso.py",
            "-p",
            self.sys_data["mpac_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "-s",
            "-m",
            "1000",
            "-t",
            "0.0001",
            "--alpha",
            "0.001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("use scikit-learn module", output)
        self.assertIn("RRMS=", output)

    def test_lasso_db_success(self):
        """Lasso regression calculates partial charges from database file."""
        sys.argv = [
            "pdf-pop-lasso.py",
            "-d",
            self.sys_data["db_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "-m",
            "1000",
            "-t",
            "0.0001",
            "--alpha",
            "0.001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertIn("charge=0.000", out.getvalue())

    def test_lasso_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-lasso.py",
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_lasso_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-lasso.py",
            "-p",
            self.sys_data["mpac_path"],
            "-d",
            self.sys_data["db_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_lasso_nonexistent_param_fails(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = [
            "pdf-pop-lasso.py",
            "-p",
            os.path.join(self.tmpdir.name, "missing.mpac"),
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()

    def test_lasso_dimension_mismatch_fails(self):
        """Mismatch between design matrix and target vector raises AssertionError."""
        bad_vtr = pdf.Vector([0.5, 0.0])  # Size 2 instead of 4
        bad_vtr_path = os.path.join(self.tmpdir.name, "bad.vtr")
        bad_vtr.save(bad_vtr_path)

        sys.argv = [
            "pdf-pop-lasso.py",
            "-p",
            self.sys_data["mpac_path"],
            self.sys_data["mat_path"],
            bad_vtr_path,
            "-m",
            "1000",
            "-t",
            "0.0001",
        ]
        with self.assertRaises(AssertionError):
            self.mod.main()


@unittest.skipUnless(HAS_PANDAS, "pandas is required")
class TestPdfPopElasticNet(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self._orig_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.sys_data = _create_sample_system(self.tmpdir.name)
        self.mod = _load_script("pdf-pop-elasticnet.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_elasticnet_fit_th_regularization(self):
        """ElasticNet fit_TH approximates exact solution with small alpha and shrinks.
        """
        mat = pdf.Matrix()
        mat.load(self.sys_data["mat_path"])
        vtr = pdf.Vector()
        vtr.load(self.sys_data["vtr_path"])

        # Small alpha -> close to known solution
        enet_weak = self.mod.ElasticNet(alpha=1e-6, l1_ratio=0.5, max_iter=2000, tol=1e-5)
        enet_weak.fit_TH(mat, vtr)
        self.assertTrue(enet_weak.is_converged)
        coef_weak = list(enet_weak.coef[:-1])
        self.assertAlmostEqual(coef_weak[0], 0.25, delta=1e-3)
        self.assertAlmostEqual(coef_weak[1], -0.125, delta=1e-3)
        self.assertAlmostEqual(coef_weak[2], -0.125, delta=1e-3)

        # Strong alpha -> coefficients shrink
        enet_strong = self.mod.ElasticNet(alpha=5.0, l1_ratio=0.5, max_iter=2000, tol=1e-5)
        enet_strong.fit_TH(mat, vtr)
        coef_strong = list(enet_strong.coef[:-1])
        norm_weak = sum(c * c for c in coef_weak)
        norm_strong = sum(c * c for c in coef_strong)
        self.assertLess(norm_strong, norm_weak)

    def test_elasticnet_param_success(self):
        """ElasticNet regression calculates partial charges and RRMS from param file."""
        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-p",
            self.sys_data["mpac_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "-m",
            "1000",
            "-t",
            "0.0001",
            "--alpha",
            "0.001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("charge=-0.000", output)
        self.assertIn("RRMS=", output)

    def test_elasticnet_db_success(self):
        """ElasticNet regression calculates partial charges from database file."""
        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-d",
            self.sys_data["db_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "-m",
            "1000",
            "-t",
            "0.0001",
            "--alpha",
            "0.001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertIn("charge=-0.000", out.getvalue())

    def test_elasticnet_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-elasticnet.py",
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_elasticnet_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-p",
            self.sys_data["mpac_path"],
            "-d",
            self.sys_data["db_path"],
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_elasticnet_nonexistent_param_fails(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-p",
            os.path.join(self.tmpdir.name, "missing.mpac"),
            self.sys_data["mat_path"],
            self.sys_data["vtr_path"],
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()

    def test_elasticnet_dimension_mismatch_fails(self):
        """Mismatch between design matrix and target vector raises AssertionError."""
        bad_vtr = pdf.Vector([0.5, 0.0])  # Size 2 instead of 4
        bad_vtr_path = os.path.join(self.tmpdir.name, "bad.vtr")
        bad_vtr.save(bad_vtr_path)

        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-p",
            self.sys_data["mpac_path"],
            self.sys_data["mat_path"],
            bad_vtr_path,
            "-m",
            "1000",
            "-t",
            "0.0001",
        ]
        with self.assertRaises(AssertionError):
            self.mod.main()


@unittest.skipUnless(HAS_SKLEARN, "scikit-learn is required")
class TestPdfPopClasso(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self._orig_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.sys_data = _create_sample_system(self.tmpdir.name)
        self.mod = _load_script("pdf-pop-classo.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_classo_param_success(self):
        """Constrained Lasso calculates partial charges satisfying constraint."""
        out_csv = os.path.join(self.tmpdir.name, "out.csv")
        out_vtr = os.path.join(self.tmpdir.name, "model.vtr")
        sys.argv = [
            "pdf-pop-classo.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--max-alpha",
            "0.01",
            "--min-alpha",
            "0.01",
            "--csv",
            out_csv,
            "--model_vector_path",
            out_vtr,
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("total charges:", output)
        self.assertIn("0.000", output)
        self.assertIn("RRMS=", output)
        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(os.path.exists(out_vtr))

        saved_vtr = pdf.Vector()
        saved_vtr.load(out_vtr)
        self.assertEqual(len(saved_vtr), 3)
        self.assertAlmostEqual(sum(saved_vtr), 0.0, places=4)

    def test_classo_constraint_satisfied(self):
        """ConstrainedLasso.fit strictly satisfies the linear constraint."""
        cl = self.mod.ConstrainedLasso()
        X = bridge.Matrix(2, 2)
        X.set(0, 0, 10.0)
        X.set(1, 1, 1.0)
        y = bridge.Vector([10.0, 0.0])
        C = bridge.Matrix(1, 2)
        C.set(0, 0, 1.0)
        C.set(0, 1, 1.0)
        b = bridge.Vector([1.0])
        cl.lambda_max = 0.01
        cl.lambda_min = 0.01
        cl.fit(X, y, C, b)
        self.assertAlmostEqual(sum(cl.coef), 1.0, places=5)

    def test_classo_db_success(self):
        """Constrained Lasso calculates partial charges from database file."""
        sys.argv = [
            "pdf-pop-classo.py",
            "-d",
            self.sys_data["db_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--max-alpha",
            "0.01",
            "--min-alpha",
            "0.01",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertIn("total charges:  0.000", out.getvalue())

    def test_classo_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-classo.py",
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_classo_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-classo.py",
            "-p",
            self.sys_data["mpac_path"],
            "-d",
            self.sys_data["db_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_classo_nonexistent_param_fails(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = [
            "pdf-pop-classo.py",
            "-p",
            os.path.join(self.tmpdir.name, "missing.mpac"),
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()

    def test_classo_nonexistent_design_matrix_fails(self):
        """Nonexistent design matrix file raises ValueError."""
        sys.argv = [
            "pdf-pop-classo.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            os.path.join(self.tmpdir.name, "missing.mat"),
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(ValueError):
            self.mod.main()


class TestPdfPopRrms(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self._orig_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.sys_data = _create_sample_system(self.tmpdir.name)
        self._orig_argv = sys.argv
        self._orig_pytools = sys.modules.get("proteindf_pytools")

    def tearDown(self):
        sys.argv = self._orig_argv
        if self._orig_pytools is not None:
            sys.modules["proteindf_pytools"] = self._orig_pytools
        else:
            sys.modules.pop("proteindf_pytools", None)
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_rrms_import_error(self):
        """pdf-pop-rrms.py can be loaded without ModuleNotFoundError."""
        sys.modules.pop("proteindf_pytools", None)
        sys.modules.pop("pdf_pop_rrms", None)
        mod = _load_script("pdf-pop-rrms.py")
        self.assertIsNotNone(mod)

    def test_rrms_calc_known_charges(self):
        """calc_RRMS yields 0.0 when exact ESPs match charges."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")

        param = _create_sample_param()
        atoms, _ = pdf.PopUtils._get_atomlist(param.molecule)
        q_true = [0.25, -0.125, -0.125]
        for i in range(3):
            atoms[i].charge = q_true[i]

        grid_pt = [0.0, 0.0, 2.0]
        pos = bridge.Position(grid_pt[0], grid_pt[1], grid_pt[2])
        pos *= 1.0 / 0.5291772108
        exact_esp = mod.calc_esp(pos, atoms)

        rrms = mod.calc_RRMS([pos], [exact_esp], atoms)
        self.assertAlmostEqual(rrms, 0.0, places=5)

    def test_rrms_load_charges(self):
        """load_charges correctly parses csv charge file."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")
        charges = mod.load_charges(self.sys_data["charge_csv_path"])
        self.assertEqual(len(charges), 3)
        self.assertAlmostEqual(charges[0], 0.25)
        self.assertAlmostEqual(charges[1], -0.125)
        self.assertAlmostEqual(charges[2], -0.125)

    def test_rrms_param_success(self):
        """RRMS is calculated from param file and charge csv."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")
        sys.argv = [
            "pdf-pop-rrms.py",
            "-p",
            self.sys_data["mpac_path"],
            "-e",
            self.sys_data["esp_path"],
            "-c",
            self.sys_data["charge_csv_path"],
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            mod.main()
        self.assertIn("RRMS=", out.getvalue())

    def test_rrms_db_success(self):
        """RRMS is calculated from database file."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")
        sys.argv = [
            "pdf-pop-rrms.py",
            "-d",
            self.sys_data["db_path"],
            "-e",
            self.sys_data["esp_path"],
            "-c",
            self.sys_data["charge_csv_path"],
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            mod.main()
        self.assertIn("RRMS=", out.getvalue())

    def test_rrms_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")
        sys.argv = [
            "pdf-pop-rrms.py",
            "-e",
            self.sys_data["esp_path"],
            "-c",
            self.sys_data["charge_csv_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_rrms_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")
        sys.argv = [
            "pdf-pop-rrms.py",
            "-p",
            self.sys_data["mpac_path"],
            "-d",
            self.sys_data["db_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_rrms_dimension_mismatch_fails(self):
        """Fewer charges in CSV than atoms raises AssertionError."""
        sys.modules["proteindf_pytools"] = pdf
        mod = _load_script("pdf-pop-rrms.py")
        mod.load_ESPs = lambda p: ([], [])
        bad_csv = os.path.join(self.tmpdir.name, "bad.csv")
        with open(bad_csv, "w") as f:
            f.write("0.5\n")

        sys.argv = [
            "pdf-pop-rrms.py",
            "-p",
            self.sys_data["mpac_path"],
            "-c",
            bad_csv,
        ]
        with self.assertRaises(AssertionError):
            mod.main()


class TestPdfPopResp(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self._orig_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.sys_data = _create_sample_system(self.tmpdir.name)
        self.mod = _load_script("pdf-pop-resp.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def test_resp_param_success_hyperbolic(self):
        """Hyperbolic RESP calculates charges from param file."""
        out_csv = os.path.join(self.tmpdir.name, "out.csv")
        out_vtr = os.path.join(self.tmpdir.name, "model.vtr")
        sys.argv = [
            "pdf-pop-resp.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--alpha",
            "0.0001",
            "--csv",
            out_csv,
            "--model_vector_path",
            out_vtr,
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("total charges:", output)
        self.assertIn("0.000", output)
        self.assertIn("RRMS=", output)
        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(os.path.exists(out_vtr))

        saved_vtr = pdf.Vector()
        saved_vtr.load(out_vtr)
        self.assertEqual(len(saved_vtr), 3)
        self.assertAlmostEqual(saved_vtr[0], 0.25, delta=0.01)
        self.assertAlmostEqual(saved_vtr[1], -0.125, delta=0.01)
        self.assertAlmostEqual(saved_vtr[2], -0.125, delta=0.01)

    def test_resp_param_success_quadratic(self):
        """Quadratic RESP calculates charges from param file."""
        sys.argv = [
            "pdf-pop-resp.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--use-quadratic",
            "--alpha",
            "0.0001",
            "--model_vector_path",
            os.path.join(self.tmpdir.name, "resp_quad.vtr"),
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("use_quadratc:  True", output)
        self.assertIn("total charges:", output)
        self.assertIn("0.000", output)
        self.assertIn("RRMS=", output)

    def test_resp_fit_known_charges(self):
        """Resp fit_q and fit_h converge to known charges when alpha is small."""
        mat = pdf.Matrix()
        mat.load(self.sys_data["mat_path"])
        vtr = pdf.Vector()
        vtr.load(self.sys_data["vtr_path"])

        resp_h = self.mod.Resp(alpha=1e-8, beta=0.1, max_iter=100)
        resp_h.fit_h(mat, vtr)
        self.assertTrue(resp_h.is_converged)
        coef_h = list(resp_h.coef[:-1])
        self.assertAlmostEqual(coef_h[0], 0.25, delta=1e-4)
        self.assertAlmostEqual(coef_h[1], -0.125, delta=1e-4)
        self.assertAlmostEqual(coef_h[2], -0.125, delta=1e-4)

        resp_q = self.mod.Resp(alpha=1e-8, beta=0.1, max_iter=100)
        resp_q.fit_q(mat, vtr)
        self.assertTrue(resp_q.is_converged)
        coef_q = list(resp_q.coef[:-1])
        self.assertAlmostEqual(coef_q[0], 0.25, delta=1e-4)
        self.assertAlmostEqual(coef_q[1], -0.125, delta=1e-4)
        self.assertAlmostEqual(coef_q[2], -0.125, delta=1e-4)

    def test_resp_db_success(self):
        """RESP charges are calculated from database file."""
        sys.argv = [
            "pdf-pop-resp.py",
            "-d",
            self.sys_data["db_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
            "--alpha",
            "0.0001",
        ]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertIn("total charges:  0.000", out.getvalue())

    def test_resp_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-resp.py",
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_resp_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.argv = [
            "pdf-pop-resp.py",
            "-p",
            self.sys_data["mpac_path"],
            "-d",
            self.sys_data["db_path"],
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_resp_nonexistent_param_fails(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = [
            "pdf-pop-resp.py",
            "-p",
            os.path.join(self.tmpdir.name, "missing.mpac"),
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()

    def test_resp_nonexistent_design_matrix_fails(self):
        """Nonexistent design matrix file raises ValueError."""
        sys.argv = [
            "pdf-pop-resp.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            os.path.join(self.tmpdir.name, "missing.mat"),
            "--target_vector_path",
            self.sys_data["vtr_path"],
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(ValueError):
            self.mod.main()

    def test_resp_dimension_mismatch_fails(self):
        """Mismatch between design matrix and target vector raises AssertionError."""
        bad_vtr = pdf.Vector([0.5, 0.0])  # Size 2 instead of 4
        bad_vtr_path = os.path.join(self.tmpdir.name, "bad.vtr")
        bad_vtr.save(bad_vtr_path)

        sys.argv = [
            "pdf-pop-resp.py",
            "-p",
            self.sys_data["mpac_path"],
            "--design_matrix_path",
            self.sys_data["mat_path"],
            "--target_vector_path",
            bad_vtr_path,
            "-e",
            self.sys_data["esp_path"],
        ]
        with self.assertRaises(AssertionError):
            self.mod.main()


if __name__ == "__main__":
    unittest.main()
