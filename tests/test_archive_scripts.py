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
from contextlib import redirect_stderr, redirect_stdout

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
    rich_progress_mock = mock.MagicMock()
    rich_progress_mock.track = lambda seq, *args, **kwargs: seq
    rich_mock.progress = rich_progress_mock
    sys.modules["rich"] = rich_mock
    sys.modules["rich.logging"] = rich_logging_mock
    sys.modules["rich.progress"] = rich_progress_mock

import numpy as np
import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.matrix import Matrix, SymmetricMatrix
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

        mat_p = SymmetricMatrix(2)
        mat_p.set(0, 0, 1.0)
        mat_p.save(param.get_density_matrix_path("rks", itr))

        mat_f = SymmetricMatrix(2)
        mat_f.set(0, 0, 0.5)
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
    # pdf-archive-h5.py tests
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
        for itr in range(1, h5.iterations + 1):
            c = h5.get_c_matrix("rks", itr)
            self.assertIsInstance(c, Matrix)
            self.assertEqual(c.rows, 2)
            self.assertEqual(c.cols, 2)
            p = h5.get_density_matrix("rks", itr)
            self.assertIsInstance(p, SymmetricMatrix)
            self.assertEqual(p.rows, 2)
            self.assertEqual(p.cols, 2)
            f = h5.get_f_matrix("rks", itr)
            self.assertIsInstance(f, SymmetricMatrix)
            self.assertEqual(f.rows, 2)
            self.assertEqual(f.cols, 2)
        h5._close()

    def test_pdf_archive_h5_pop_mulliken_default_no_pop(self):
        """pdf-archive-h5.py does not calculate or archive Mulliken population by default."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)

        mpac_path = os.path.join(self.tmpdir.name, "pop_default.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "pop_default.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]

        with mock.patch.object(mod.subprocess, "run") as mock_run:
            mod.main()
            mock_run.assert_not_called()

        self.assertTrue(os.path.exists(h5_path))
        h5 = PdfParam_H5()
        h5.open(h5_path)
        self.assertIsNone(h5.get_pop_mulliken_atom("rks", 1))
        self.assertIsNone(h5.get_pop_mulliken_atom("rks", 2))
        h5._close()

    def test_pdf_archive_h5_pop_mulliken_converged(self):
        """pdf-archive-h5.py --pop calculates Mulliken population for converged solution only."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)

        mpac_path = os.path.join(self.tmpdir.name, "pop_conv.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "pop_conv.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path, "--pop"]

        dummy_output = (
            " Mulliken Population\n"
            " total atomic charge =    2.00\n"
            "          Net charge =    0.00\n\n"
            " Gross atom population\n"
            "  index Atom    gross         net\n"
            "      0   H       0.850000     0.150000\n"
            "      1   H       1.150000    -0.150000\n\n"
            " Mulliken population analysis (orbial population)\n"
        )

        def fake_subprocess_run(cmd, capture_output=True, text=True):
            self.assertEqual(cmd[0], "pdf")
            self.assertEqual(cmd[1], "pop-mulliken")
            self.assertIn("-p", cmd)
            self.assertIn("-P", cmd)
            self.assertIn("-S", cmd)
            self.assertNotIn("-i", cmd)
            p_file = cmd[cmd.index("-P") + 1]
            s_file = cmd[cmd.index("-S") + 1]
            mpac_file = cmd[cmd.index("-p") + 1]
            self.assertTrue(os.path.exists(p_file))
            self.assertTrue(os.path.exists(s_file))
            self.assertTrue(os.path.exists(mpac_file))
            return mock.Mock(returncode=0, stdout=dummy_output, stderr="")

        with mock.patch.object(mod.subprocess, "run", side_effect=fake_subprocess_run) as mock_run:
            mod.main()
            self.assertEqual(mock_run.call_count, 1)

        self.assertTrue(os.path.exists(h5_path))
        h5 = PdfParam_H5()
        h5.open(h5_path)
        # Iteration 1 must not be calculated
        self.assertIsNone(h5.get_pop_mulliken_atom("rks", 1))
        # Iteration 2 (converged) must be calculated and saved
        pop_loaded2 = h5.get_pop_mulliken_atom("rks", 2)
        self.assertIsInstance(pop_loaded2, Vector)
        self.assertEqual(len(pop_loaded2), 2)
        self.assertAlmostEqual(pop_loaded2[0], 0.15)
        self.assertAlmostEqual(pop_loaded2[1], -0.15)
        h5._close()

    def test_pdf_archive_h5_pop_mulliken_unconverged(self):
        """pdf-archive-h5.py --pop calculates pop for last iteration when unconverged."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=3)
        param._set_scf_converged(False)
        _create_synthetic_files(param)

        mpac_path = os.path.join(self.tmpdir.name, "pop_unconv.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "pop_unconv.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path, "--pop"]

        dummy_output = (
            " Gross atom population\n"
            "  index Atom    gross         net\n"
            "      0   H       0.800000     0.200000\n"
            "      1   H       1.200000    -0.200000\n"
        )

        def fake_subprocess_run(cmd, capture_output=True, text=True):
            self.assertEqual(cmd[0], "pdf")
            self.assertEqual(cmd[1], "pop-mulliken")
            self.assertIn("-p", cmd)
            self.assertIn("-P", cmd)
            self.assertIn("-S", cmd)
            self.assertNotIn("-i", cmd)
            return mock.Mock(returncode=0, stdout=dummy_output, stderr="")

        with mock.patch.object(mod.subprocess, "run", side_effect=fake_subprocess_run) as mock_run:
            mod.main()
            self.assertEqual(mock_run.call_count, 1)

        self.assertTrue(os.path.exists(h5_path))
        h5 = PdfParam_H5()
        h5.open(h5_path)
        self.assertIsNone(h5.get_pop_mulliken_atom("rks", 1))
        # Iteration 2 (last iteration for unconverged) must be calculated
        pop_loaded2 = h5.get_pop_mulliken_atom("rks", 2)
        self.assertIsInstance(pop_loaded2, Vector)
        self.assertEqual(len(pop_loaded2), 2)
        self.assertAlmostEqual(pop_loaded2[0], 0.20)
        self.assertAlmostEqual(pop_loaded2[1], -0.20)
        self.assertIsNone(h5.get_pop_mulliken_atom("rks", 3))
        h5._close()

    def test_pdf_archive_h5_pop_mulliken_failure_does_not_save(self):
        """pdf-archive-h5.py does not save pop if pop-mulliken fails (non-zero exit)."""
        param = _create_sample_param(work_path=self.tmpdir.name, iterations=2)
        _create_synthetic_files(param)

        mpac_path = os.path.join(self.tmpdir.name, "pop_fail.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "pop_fail.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path, "--pop"]

        with mock.patch.object(
            mod.subprocess,
            "run",
            return_value=mock.Mock(returncode=1, stdout="", stderr="error"),
        ) as mock_run:
            mod.main()
            self.assertEqual(mock_run.call_count, 1)

        self.assertTrue(os.path.exists(h5_path))
        h5 = PdfParam_H5()
        h5.open(h5_path)
        self.assertIsNone(h5.get_pop_mulliken_atom("rks", 2))
        h5._close()

    def test_parse_mulliken_pop(self):
        """parse_mulliken_pop extracts net charges from Gross atom population table."""
        mod = _load_script("pdf-archive-h5.py")
        output = (
            " Mulliken Population\n"
            " total atomic charge =    8.00\n"
            "          Net charge =    0.00\n\n"
            " Gross atom population\n"
            "  index Atom    gross         net\n"
            "      0   O       8.345678    -0.345678\n"
            "      1   H       0.827161     0.172839\n"
            "      2   H       0.827161     0.172839\n\n"
            " Mulliken population analysis (orbial population)\n"
        )
        vec = mod.parse_mulliken_pop(output)
        self.assertIsInstance(vec, Vector)
        self.assertEqual(len(vec), 3)
        self.assertAlmostEqual(vec[0], -0.345678)
        self.assertAlmostEqual(vec[1], 0.172839)
        self.assertAlmostEqual(vec[2], 0.172839)

        # Invalid or empty output
        self.assertIsNone(mod.parse_mulliken_pop(""))
        self.assertIsNone(mod.parse_mulliken_pop("No population table here"))

    def test_pdf_archive_h5_input_density_matrix(self):
        """pdf-archive-h5.py saves iteration 1 input density matrix P_in/rks_1."""
        work_dir = os.path.join(self.tmpdir.name, "fl_Work")
        os.makedirs(work_dir, exist_ok=True)
        param = _create_sample_param(work_path=work_dir, iterations=2)
        param._data.setdefault("control", {})["file_base_name"]["Ppq_matrix"] = "Ppq.%s.mat"
        _create_synthetic_files(param)

        # Write distinct input density matrix ~Ppq.rks1.mat
        mat_p_in = SymmetricMatrix(2)
        mat_p_in.set(0, 0, 3.14)
        mat_p_in.set(1, 0, 1.41)
        mat_p_in.set(1, 1, 2.71)
        in_path = param.get_input_density_matrix_path("rks", 1)
        self.assertTrue(os.path.basename(in_path).startswith("~Ppq"))
        mat_p_in.save(in_path)

        mpac_path = os.path.join(self.tmpdir.name, "param.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "archive_with_pin.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]
        mod.main()

        self.assertTrue(os.path.exists(h5_path))
        import h5py

        with h5py.File(h5_path, "r") as h5_raw:
            self.assertIn("P_in/rks_1", h5_raw)
            self.assertIn("P/rks_1", h5_raw)
            self.assertIn("P/rks_2", h5_raw)
            self.assertNotIn("P_in/rks_2", h5_raw)
            self.assertTrue(np.allclose(h5_raw["P_in/rks_1"][:], mat_p_in.data))

        h5 = PdfParam_H5()
        h5.open(h5_path)
        loaded_p_in = h5.get_input_density_matrix("rks", 1)
        self.assertIsInstance(loaded_p_in, SymmetricMatrix)
        self.assertTrue(np.allclose(loaded_p_in.data, mat_p_in.data))
        for itr in range(1, h5.iterations + 1):
            p = h5.get_density_matrix("rks", itr)
            self.assertIsInstance(p, SymmetricMatrix)
            self.assertEqual(p.rows, 2)
        h5._close()

    def test_pdf_archive_h5_input_density_matrix_missing(self):
        """pdf-archive-h5.py continues with a warning when ~Ppq.rks1.mat is missing."""
        work_dir = os.path.join(self.tmpdir.name, "fl_Work")
        os.makedirs(work_dir, exist_ok=True)
        param = _create_sample_param(work_path=work_dir, iterations=2)
        param._data.setdefault("control", {})["file_base_name"]["Ppq_matrix"] = "Ppq.%s.mat"
        _create_synthetic_files(param)

        # Do not create ~Ppq.rks1.mat
        in_path = param.get_input_density_matrix_path("rks", 1)
        self.assertFalse(os.path.exists(in_path))

        mpac_path = os.path.join(self.tmpdir.name, "param_no_pin.mpac")
        pdf.save_pdfparam(param, mpac_path)

        h5_path = os.path.join(self.tmpdir.name, "archive_no_pin.h5")
        mod = _load_script("pdf-archive-h5.py")
        sys.argv = ["pdf-archive-h5.py", "-p", mpac_path, "-o", h5_path]

        with self.assertLogs(mod.logger, level="WARNING") as log:
            mod.main()
            self.assertTrue(
                any("input P not found" in record.getMessage() for record in log.records)
            )

        self.assertTrue(os.path.exists(h5_path))
        import h5py

        with h5py.File(h5_path, "r") as h5_raw:
            self.assertNotIn("P_in/rks_1", h5_raw)
            self.assertIn("P/rks_1", h5_raw)
            self.assertIn("P/rks_2", h5_raw)

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
        for itr in range(1, h5.iterations + 1):
            c = h5.get_c_matrix("rks", itr)
            self.assertIsInstance(c, Matrix)
            self.assertEqual(c.rows, 2)
            self.assertEqual(c.cols, 2)
            p = h5.get_density_matrix("rks", itr)
            self.assertIsInstance(p, SymmetricMatrix)
            self.assertEqual(p.rows, 2)
            self.assertEqual(p.cols, 2)
            f = h5.get_f_matrix("rks", itr)
            self.assertIsInstance(f, SymmetricMatrix)
            self.assertEqual(f.rows, 2)
            self.assertEqual(f.cols, 2)
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
