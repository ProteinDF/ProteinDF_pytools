#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import io
import logging
import os
import sys
import tempfile
import unittest
import unittest.mock as mock
from contextlib import redirect_stdout

# Mock rich if it is not available in the test environment
try:
    import rich
    import rich.logging
    import rich.progress
except ImportError:
    class FakeRichHandler(logging.Handler):
        def __init__(self, *args, **kwargs):
            super().__init__()

        def emit(self, record):
            pass

    rich_mock = mock.MagicMock()
    rich_mock.__path__ = []
    rich_logging_mock = mock.MagicMock()
    rich_logging_mock.RichHandler = FakeRichHandler
    rich_mock.logging = rich_logging_mock
    sys.modules["rich"] = rich_mock
    sys.modules["rich.logging"] = rich_logging_mock
    sys.modules["rich.progress"] = mock.MagicMock()

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.matrix import Matrix, SymmetricMatrix
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5
from proteindf_tools.vector import Vector

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _create_sample_param(work_path=None, iterations=2):
    """Create minimal PdfParam representing an H atom with STO-3G basis."""
    param = PdfParam()
    param.comment = "synthetic test system"
    param.num_of_atoms = 1
    param.num_of_AOs = 2
    param.num_of_MOs = 2
    param.method = "rks"
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
        "occupation_vtr": "occ.%s.vtr",
        "Spq_matrix": "s.mat",
        "Hpq_matrix": "h.mat",
        "Hpq2_matrix": "h2.mat",
        "Fpq_matrix": "f.%s.mat",
        "C_matrix": "c.%s.mat",
        "Ppq_matrix": "p.%s.mat",
    }
    return param


def _create_synthetic_files(param):
    """Write synthetic vector and matrix files required for archive."""
    v_occ = Vector(2)
    v_occ[0] = 2.0
    v_occ[1] = 0.0
    v_occ.save(param.get_occ_path("rks"))

    for itr in range(1, param.iterations + 1):
        v_eig = Vector(2)
        v_eig[0] = -0.5 * itr
        v_eig[1] = 0.2 * itr
        v_eig.save(param.get_energy_level_path("rks", itr))

        mat_c = Matrix(2, 2)
        mat_c.set(0, 0, 1.0)
        mat_c.set(1, 1, 1.0)
        mat_c.save(param.get_c_mat_path("rks", itr))

        mat_p = Matrix(2, 2)
        mat_p.save(param.get_density_matrix_path("rks", itr))

        mat_f = Matrix(2, 2)
        mat_f.save(param.get_f_mat_path("rks", itr))

    s_mat = SymmetricMatrix(2)
    s_mat.set(0, 0, 1.0)
    s_mat.set(1, 1, 1.0)
    s_mat.save(param.get_s_mat_path())
    s_mat.save(param.get_h_mat_path())
    s_mat.save(param.get_h2_mat_path())


class TestArchiveScripts(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        root_logger = logging.getLogger()
        for h in list(root_logger.handlers):
            h.close()
            root_logger.removeHandler(h)

    def tearDown(self):
        # Remove and close any handlers created by logging in scripts
        root_logger = logging.getLogger()
        for h in list(root_logger.handlers):
            h.close()
            root_logger.removeHandler(h)
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    # -------------------------------------------------------------------------
    # 1. pdf-archive.py tests
    # -------------------------------------------------------------------------

    def test_pdf_archive_with_param_option(self):
        """pdf-archive.py with -p <path> reads parameter and generates DB."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)
        mpac_path = os.path.join(self.tmpdir.name, "custom.mpac")
        pdf.save_pdfparam(param, mpac_path)

        db_path = os.path.join(self.tmpdir.name, "result.db")
        mod = _load_script("pdf-archive.py")
        sys.argv = ["pdf-archive.py", "-p", mpac_path, "-o", db_path]
        mod.main()

        self.assertTrue(os.path.exists(db_path))
        archive = PdfArchive_Sqlite3(db_path)
        self.assertEqual(archive.iterations, 2)
        e1 = archive.get_energylevel("rks", 1)
        self.assertAlmostEqual(e1[0], -0.5)

    def test_pdf_archive_default_param(self):
        """pdf-archive.py without -p reads pdfparam.mpac in current working directory."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)
        mpac_path = os.path.join(self.tmpdir.name, "pdfparam.mpac")
        pdf.save_pdfparam(param, mpac_path)

        db_path = os.path.join(self.tmpdir.name, "result_default.db")
        mod = _load_script("pdf-archive.py")
        sys.argv = ["pdf-archive.py", "-o", db_path]
        mod.main()

        self.assertTrue(os.path.exists(db_path))
        archive = PdfArchive_Sqlite3(db_path)
        self.assertEqual(archive.iterations, 2)

    def test_pdf_archive_nonexistent_param(self):
        """pdf-archive.py handles nonexistent mpac without creating DB."""
        db_path = os.path.join(self.tmpdir.name, "result_nonexistent.db")
        mod = _load_script("pdf-archive.py")
        sys.argv = ["pdf-archive.py", "-p", "nonexistent.mpac", "-o", db_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertFalse(os.path.exists(db_path))
        self.assertIn("FileNotFoundError", buf.getvalue())

    # -------------------------------------------------------------------------
    # 2. pdf-archive-h5.py tests
    # -------------------------------------------------------------------------

    def test_pdf_archive_h5_with_param_option(self):
        """pdf-archive-h5.py with -p <path> reads parameter and generates HDF5."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)
        mpac_path = os.path.join(self.tmpdir.name, "custom.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "result.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]
        mod.main()

        self.assertTrue(os.path.exists(h5_path))
        h5 = PdfParam_H5()
        h5.open(h5_path)
        self.assertEqual(h5.iterations, 2)
        e1 = h5.get_energy_level("rks", 1)
        self.assertAlmostEqual(e1[0], -0.5)
        tes = h5.get_TEs()
        self.assertAlmostEqual(tes[0], -0.5)
        self.assertAlmostEqual(tes[1], -1.0)
        h5._close()

    def test_pdf_archive_h5_default_param(self):
        """pdf-archive-h5.py without -p reads pdfparam.mpac in current working directory."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)
        mpac_path = os.path.join(self.tmpdir.name, "pdfparam.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "result_default.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-o", h5_path]
        mod.main()

        self.assertTrue(os.path.exists(h5_path))
        h5 = PdfParam_H5()
        h5.open(h5_path)
        self.assertEqual(h5.iterations, 2)
        h5._close()

    def test_pdf_archive_h5_nonexistent_param(self):
        """pdf-archive-h5.py handles nonexistent mpac gracefully without crashing."""
        h5_path = os.path.join(self.tmpdir.name, "result_nonexistent.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", "nonexistent.mpac", "-o", h5_path]
        mod.main()

        self.assertFalse(os.path.exists(h5_path))

    # -------------------------------------------------------------------------
    # 3. pdf-info-iteration.py tests
    # -------------------------------------------------------------------------

    def test_pdf_info_iteration_param(self):
        """pdf-info-iteration.py outputs iterations from mpac parameter file."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=3)
        mpac_path = os.path.join(self.tmpdir.name, "info.mpac")
        pdf.save_pdfparam(param, mpac_path)

        mod = _load_script("pdf-info-iteration.py")
        sys.argv = ["pdf-info-iteration.py", "-p", mpac_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertEqual(buf.getvalue().strip(), "3")

    def test_pdf_info_iteration_db(self):
        """pdf-info-iteration.py outputs iterations from HDF5 results file."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=4)
        _create_synthetic_files(param)
        mpac_path = os.path.join(self.tmpdir.name, "temp.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "info.h5")
        mod_arch = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]
        mod_arch.main()

        mod = _load_script("pdf-info-iteration.py")
        sys.argv = ["pdf-info-iteration.py", "-d", h5_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertEqual(buf.getvalue().strip(), "4")

    def test_pdf_info_iteration_nonexistent(self):
        """pdf-info-iteration.py raises FileNotFoundError when parameter file is missing."""
        mod = _load_script("pdf-info-iteration.py")
        sys.argv = ["pdf-info-iteration.py", "-p", "nonexistent.mpac"]
        with self.assertRaises(FileNotFoundError):
            mod.main()

    # -------------------------------------------------------------------------
    # 4. pdf-mat-info.py tests
    # -------------------------------------------------------------------------

    def test_pdf_mat_info_std(self):
        """pdf-mat-info.py displays row and column size of standard Matrix."""
        mat_path = os.path.join(self.tmpdir.name, "std.mat")
        m = Matrix(2, 3)
        m.save(mat_path)

        mod = _load_script("pdf-mat-info.py")
        sys.argv = ["pdf-mat-info.py", mat_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        output = buf.getvalue()
        self.assertIn("standard dens matrix", output)
        self.assertIn("size: 2 x 3", output)

    def test_pdf_mat_info_sym(self):
        """pdf-mat-info.py displays row and column size of SymmetricMatrix."""
        sym_path = os.path.join(self.tmpdir.name, "sym.mat")
        s = SymmetricMatrix(3)
        s.save(sym_path)

        mod = _load_script("pdf-mat-info.py")
        sys.argv = ["pdf-mat-info.py", sym_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        output = buf.getvalue()
        self.assertIn("symmetric dens matrix", output)
        self.assertIn("size: 3 x 3", output)

    def test_pdf_mat_info_nonexistent(self):
        """pdf-mat-info.py prints error message without crashing for nonexistent file."""
        mod = _load_script("pdf-mat-info.py")
        sys.argv = ["pdf-mat-info.py", os.path.join(self.tmpdir.name, "nonexistent.mat")]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertIn("cannot load file.", buf.getvalue())

    def test_pdf_mat_info_corrupted(self):
        """pdf-mat-info.py prints error message without crashing for corrupted file."""
        corrupt_path = os.path.join(self.tmpdir.name, "corrupt.mat")
        with open(corrupt_path, "wb") as f:
            f.write(b"not a valid matrix binary")

        mod = _load_script("pdf-mat-info.py")
        sys.argv = ["pdf-mat-info.py", corrupt_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertIn("cannot load file.", buf.getvalue())

    # -------------------------------------------------------------------------
    # 5. pdf-mat-show.py tests
    # -------------------------------------------------------------------------

    def test_pdf_mat_show_std(self):
        """pdf-mat-show.py displays contents of standard Matrix."""
        mat_path = os.path.join(self.tmpdir.name, "std_show.mat")
        m = Matrix(2, 2)
        m.set(0, 0, 1.234)
        m.save(mat_path)

        mod = _load_script("pdf-mat-show.py")
        sys.argv = ["pdf-mat-show.py", mat_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        output = buf.getvalue()
        self.assertIn("standard dens matrix", output)
        self.assertIn("size: 2 x 2", output)
        self.assertIn("1.234", output)

    def test_pdf_mat_show_sym(self):
        """pdf-mat-show.py displays contents of SymmetricMatrix."""
        sym_path = os.path.join(self.tmpdir.name, "sym_show.mat")
        s = SymmetricMatrix(2)
        s.set(0, 0, 5.678)
        s.save(sym_path)

        mod = _load_script("pdf-mat-show.py")
        sys.argv = ["pdf-mat-show.py", sym_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        output = buf.getvalue()
        self.assertIn("symmetric dens matrix", output)
        self.assertIn("size: 2 x 2", output)
        self.assertIn("5.678", output)

    def test_pdf_mat_show_nonexistent(self):
        """pdf-mat-show.py prints error message without crashing for nonexistent file."""
        mod = _load_script("pdf-mat-show.py")
        sys.argv = ["pdf-mat-show.py", os.path.join(self.tmpdir.name, "nonexistent.mat")]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertIn("cannot load file.", buf.getvalue())

    def test_pdf_mat_show_corrupted(self):
        """pdf-mat-show.py prints error message without crashing for corrupted file."""
        corrupt_path = os.path.join(self.tmpdir.name, "corrupt.mat")
        with open(corrupt_path, "wb") as f:
            f.write(b"not a valid matrix binary")

        mod = _load_script("pdf-mat-show.py")
        sys.argv = ["pdf-mat-show.py", corrupt_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            mod.main()

        self.assertIn("cannot load file.", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
