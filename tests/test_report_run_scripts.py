#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cProfile
import importlib.util
import io
import logging
import os
import pstats
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

try:
    import rich
    import rich.logging

    HAS_RICH = True
except ImportError:
    HAS_RICH = False

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _create_sample_param(iterations=2, te=-75.0):
    """Create minimal synthetic PdfParam representing water (H2O: O, H, H)."""
    param = PdfParam()
    param.comment = "synthetic water system"
    param.num_of_atoms = 3
    param.num_of_AOs = 3
    param.num_of_MOs = 3
    param.method = "rks"
    param.iterations = iterations
    param.guess = "none"
    param._set_scf_converged(True)

    mol = bridge.AtomGroup()
    a_o = bridge.Atom(symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="O1")
    a_h1 = bridge.Atom(symbol="H", position=bridge.Position([0.0, 1.0, 0.0]), charge=1.0, name="H1")
    a_h2 = bridge.Atom(symbol="H", position=bridge.Position([1.0, 0.0, 0.0]), charge=1.0, name="H2")
    mol.set_atom("0", a_o)
    mol.set_atom("1", a_h1)
    mol.set_atom("2", a_h2)
    param.molecule = mol

    bs_o = BasisSet()
    bs_o.name = "sto-3g-O"
    cgto_o = ContractedGTO("s", 1)
    cgto_o[0] = PrimitiveGTO(71.6168, 0.15433)
    bs_o.append(cgto_o)

    bs_h = BasisSet()
    bs_h.name = "sto-3g-H"
    cgto_h = ContractedGTO("s", 1)
    cgto_h[0] = PrimitiveGTO(3.42525, 0.15433)
    bs_h.append(cgto_h)

    param.set_basisset("O", bs_o)
    param.set_basisset("O@O1", bs_o)
    param.set_basisset("H", bs_h)
    param.set_basisset("H@H1", bs_h)
    param.set_basisset("H@H2", bs_h)

    for itr in range(1, iterations + 1):
        param.set_total_energy(itr, te - 0.1 * (itr - 1))

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


def _create_sample_db(db_path, iterations=2, te=-75.0):
    """Create synthetic SQLite3 DB archive."""
    param = _create_sample_param(iterations=iterations, te=te)
    arc = PdfArchive_Sqlite3(db_path, pdfparam=param)
    arc.set_occupations("rks", pdf.Vector([2.0, 0.0, 0.0]))
    for itr in range(1, iterations + 1):
        arc.set_energylevel("rks", itr, pdf.Vector([-0.5 - 0.01 * itr, 0.2, 0.5]))
    return arc


def _create_sample_h5(h5_path, iterations=2, te=-75.0):
    """Create synthetic HDF5 archive."""
    param = _create_sample_param(iterations=iterations, te=te)
    h = PdfParam_H5(param)
    h.open(h5_path)
    h.save_model()
    tes = pdf.Vector([te - 0.1 * i for i in range(iterations)])
    h.save_TEs(tes)
    h.save_occ("rks", pdf.Vector([2.0, 0.0, 0.0]))
    for itr in range(1, iterations + 1):
        h.save_energy_level("rks", itr, pdf.Vector([-0.5 - 0.01 * itr, 0.2, 0.5]))
    h._close()
    return param


class BaseScriptTestCase(unittest.TestCase):
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
        root_logger = logging.getLogger()
        for h in list(root_logger.handlers):
            h.close()
            root_logger.removeHandler(h)
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()


# =============================================================================
# 1. pdf-test.py
# =============================================================================
class TestPdfTest(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-test.py")
        self.db1 = os.path.join(self.tmpdir.name, "result1.db")
        self.db2 = os.path.join(self.tmpdir.name, "result2.db")
        self.db_diff = os.path.join(self.tmpdir.name, "result_diff.db")
        _create_sample_db(self.db1, te=-75.0)
        _create_sample_db(self.db2, te=-75.0)
        _create_sample_db(self.db_diff, te=-76.0)

    def test_same_data(self):
        """Comparing identical SQLite3 results exits with 0 (OK)."""
        sys.argv = ["pdf-test.py", self.db1, self.db2]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 0)

    def test_different_energy(self):
        """Comparing results with different total energy exits with 1 (NG)."""
        sys.argv = ["pdf-test.py", self.db1, self.db_diff]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 1)

    def test_compare_info_option(self):
        """--compare-info checks system information between archives."""
        sys.argv = ["pdf-test.py", "--compare-info", self.db1, self.db2]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 0)

    def test_compare_energy_option(self):
        """--compare-energy detects total energy discrepancy."""
        sys.argv = ["pdf-test.py", "--compare-energy", self.db1, self.db_diff]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 1)

    def test_nonexistent_file(self):
        """Comparing with nonexistent file exits with 'file not found'."""
        sys.argv = ["pdf-test.py", "nonexistent.db", self.db1]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertIn("file not found", str(cm.exception))


# =============================================================================
# 2. pdf-test-h5.py
# =============================================================================
class TestPdfTestH5(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-test-h5.py")
        self.h5_1 = os.path.join(self.tmpdir.name, "result1.h5")
        self.h5_2 = os.path.join(self.tmpdir.name, "result2.h5")
        self.h5_diff = os.path.join(self.tmpdir.name, "result_diff.h5")
        _create_sample_h5(self.h5_1, te=-75.0)
        _create_sample_h5(self.h5_2, te=-75.0)
        _create_sample_h5(self.h5_diff, te=-76.0)

    def test_same_data(self):
        """Comparing identical HDF5 results exits with 0 (OK)."""
        sys.argv = ["pdf-test-h5.py", self.h5_1, self.h5_2]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 0)

    def test_different_energy(self):
        """Comparing HDF5 results with different total energy exits with 1 (NG)."""
        sys.argv = ["pdf-test-h5.py", self.h5_1, self.h5_diff]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 1)

    def test_compare_info_option(self):
        """--compare-info checks system information between HDF5 files."""
        sys.argv = ["pdf-test-h5.py", "--compare-info", self.h5_1, self.h5_2]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 0)

    def test_compare_energy_option(self):
        """--compare-energy detects total energy discrepancy in HDF5."""
        sys.argv = ["pdf-test-h5.py", "--compare-energy", self.h5_1, self.h5_diff]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 1)

    def test_nonexistent_file(self):
        """Comparing with nonexistent HDF5 file exits with 'file not found'."""
        sys.argv = ["pdf-test-h5.py", "nonexistent.h5", self.h5_1]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertIn("file not found", str(cm.exception))


# =============================================================================
# 3. pdf-report.py
# =============================================================================
class TestPdfReport(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-report.py")
        self.db = os.path.join(self.tmpdir.name, "pdfresults.db")
        _create_sample_db(self.db, iterations=2, te=-75.0)

    def test_report_generation(self):
        """pdf-report.py generates convergence dat and png files in output directory."""
        out_dir = os.path.join(self.tmpdir.name, "report_dir")
        sys.argv = ["pdf-report.py", self.db, out_dir]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        expected_files = [
            "TE_hist.dat",
            "TE_hist.png",
            "eigvals_hist.dat",
            "eigvals_hist.png",
            "eigvals_last.png",
        ]
        for f in expected_files:
            self.assertTrue(
                os.path.exists(os.path.join(out_dir, f)),
                f"Expected file {f} was not generated",
            )
        self.assertIn("ProteinDF report:", buf.getvalue())

    def test_profile_option(self):
        """pdf-report.py with -p generates profiling stat file."""
        out_dir = os.path.join(self.tmpdir.name, "report_prof_dir")
        stat_file = os.path.join(self.tmpdir.name, "profile.stat")
        sys.argv = ["pdf-report.py", self.db, out_dir, "-p", stat_file]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        self.assertTrue(os.path.exists(stat_file))

    def test_nonexistent_db(self):
        """pdf-report.py raises error or stops when db file does not exist."""
        out_dir = os.path.join(self.tmpdir.name, "report_err")
        sys.argv = ["pdf-report.py", "nonexistent.db", out_dir]
        with self.assertRaises(Exception):
            self.mod.main()

    def test_cannot_create_output_dir(self):
        """pdf-report.py fails if output directory path exists as a regular file."""
        file_as_dir = os.path.join(self.tmpdir.name, "file_as_dir")
        with open(file_as_dir, "w") as f:
            f.write("blocking file")
        sys.argv = ["pdf-report.py", self.db, file_as_dir]
        with self.assertRaises(FileExistsError):
            self.mod.main()


# =============================================================================
# 4. pdf-report-h5.py
# =============================================================================
@unittest.skipUnless(HAS_RICH, "rich is not installed")
class TestPdfReportH5(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-report-h5.py")
        self.h5 = os.path.join(self.tmpdir.name, "pdfresults.h5")
        _create_sample_h5(self.h5, iterations=2, te=-75.0)

    def test_report_generation(self):
        """pdf-report-h5.py generates report plots and data from HDF5."""
        out_dir = os.path.join(self.tmpdir.name, "report_h5_dir")
        sys.argv = ["pdf-report-h5.py", self.h5, out_dir]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        expected_files = [
            "TE_hist.dat",
            "TE_hist.png",
            "eigvals_hist.rks.dat",
            "eigvals_hist.rks.png",
            "eigvals_last.rks.png",
        ]
        for f in expected_files:
            self.assertTrue(
                os.path.exists(os.path.join(out_dir, f)),
                f"Expected file {f} was not generated",
            )
        self.assertIn("ProteinDF report:", buf.getvalue())

    def test_iteration_option(self):
        """pdf-report-h5.py with -i uses specified iteration for report."""
        out_dir = os.path.join(self.tmpdir.name, "report_h5_iter")
        sys.argv = ["pdf-report-h5.py", self.h5, out_dir, "-i", "1"]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        self.assertTrue(os.path.exists(os.path.join(out_dir, "TE_hist.png")))

    @unittest.expectedFailure
    def test_profile_option(self):
        """pdf-report-h5.py with -p runs cProfile profiling.

        scripts/pdf-report-h5.py:108: cProfile.runctx("report()", ...) で
        PdfReport インスタンスを関数として呼び出そうとするため
        TypeError: 'PdfReport' object is not callable が発生する。
        修正されたら外す。
        """
        out_dir = os.path.join(self.tmpdir.name, "report_h5_prof")
        stat_file = os.path.join(self.tmpdir.name, "pdf-report.stat")
        sys.argv = ["pdf-report-h5.py", self.h5, out_dir, "-p", stat_file]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        self.assertTrue(os.path.exists(stat_file))

    def test_nonexistent_h5(self):
        """pdf-report-h5.py raises error or stops when HDF5 file does not exist."""
        out_dir = os.path.join(self.tmpdir.name, "report_err")
        sys.argv = ["pdf-report-h5.py", "nonexistent.h5", out_dir]
        with self.assertRaises(Exception):
            self.mod.main()

    def test_cannot_create_output_dir(self):
        """pdf-report-h5.py fails if output directory path exists as a regular file."""
        file_as_dir = os.path.join(self.tmpdir.name, "file_as_dir")
        with open(file_as_dir, "w") as f:
            f.write("blocking file")
        sys.argv = ["pdf-report-h5.py", self.h5, file_as_dir]
        with self.assertRaises(FileExistsError):
            self.mod.main()


# =============================================================================
# 5. pdf-reg-harris.py
# =============================================================================
class TestPdfRegHarris(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-reg-harris.py")
        self.param = _create_sample_param()
        self.param_path = os.path.join(self.tmpdir.name, "pdfparam.mpac")
        pdf.save_pdfparam(self.param, self.param_path)

        # 3 orbitals: 0: O@O1, 1: H@H1, 2: H@H2
        self.dens_mat = pdf.SymmetricMatrix(3)
        self.dens_mat.set(0, 0, 1.25)
        self.dens_mat.set(1, 0, 0.10)
        self.dens_mat.set(1, 1, 0.45)
        self.dens_mat.set(2, 0, 0.15)
        self.dens_mat.set(2, 1, 0.05)
        self.dens_mat.set(2, 2, 0.55)
        self.dens_path = os.path.join(self.tmpdir.name, "density.mat")
        self.dens_mat.save(self.dens_path)

    def test_reg_harris_success(self):
        """pdf-reg-harris.py extracts target atom density matrix and registers to DB."""
        out_mpac = os.path.join(self.tmpdir.name, "custom_harris.mpac")
        sys.argv = [
            "pdf-reg-harris.py",
            self.dens_path,
            "O@O1",
            "-o",
            out_mpac,
            "-p",
            self.param_path,
        ]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        self.assertTrue(os.path.exists(out_mpac))
        db = bridge.load_msgpack(out_mpac)
        self.assertIn("basis_sets", db)
        self.assertIn("density_matrix", db)
        self.assertIn("O", db["basis_sets"])
        self.assertIn("O", db["density_matrix"])
        dens_o = db["density_matrix"]["O"]
        # Target orbital for O@O1 is index 0; value at (0, 0) is 1.25
        self.assertEqual(dens_o["row"], 1)
        self.assertEqual(dens_o["col"], 1)
        self.assertAlmostEqual(dens_o["data"][0], 1.25)

    @unittest.expectedFailure
    def test_reg_harris_default_output(self):
        """pdf-reg-harris.py without -o writes output to default harris.mpac.

        scripts/pdf-reg-harris.py:57: args.db_path[0] でデフォルト文字列
        'harris.mpac' の1文字目('h')がファイル名として扱われてしまう不具合があるため失敗する。
        修正されたら外す。
        """
        sys.argv = ["pdf-reg-harris.py", self.dens_path, "O@O1", "-p", self.param_path]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        self.assertTrue(os.path.exists("harris.mpac"))

    def test_nonexistent_matrix(self):
        """pdf-reg-harris.py stops when density matrix file does not exist."""
        out_mpac = os.path.join(self.tmpdir.name, "out.mpac")
        sys.argv = [
            "pdf-reg-harris.py",
            "nonexistent.mat",
            "O@O1",
            "-o",
            out_mpac,
            "-p",
            self.param_path,
        ]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(AssertionError):
                self.mod.main()

    def test_mismatched_atomlabel(self):
        """pdf-reg-harris.py stops when atom label does not match basisset."""
        out_mpac = os.path.join(self.tmpdir.name, "out.mpac")
        sys.argv = [
            "pdf-reg-harris.py",
            self.dens_path,
            "UNKNOWN@X",
            "-o",
            out_mpac,
            "-p",
            self.param_path,
        ]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(Exception):
                self.mod.main()


# =============================================================================
# 6. pdfcmd.py
# =============================================================================
class TestPdfCmd(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdfcmd.py")
        # Create a mock executable in temporary directory
        self.subcmd_path = os.path.join(self.tmpdir.name, "pdf-testrun")
        with open(self.subcmd_path, "w") as f:
            f.write("#!/bin/sh\nexit 0\n")
        os.chmod(self.subcmd_path, 0o755)

        self._orig_path = os.environ.get("PATH", "")
        os.environ["PATH"] = self.tmpdir.name + os.pathsep + self._orig_path

    def tearDown(self):
        os.environ["PATH"] = self._orig_path
        super().tearDown()

    def test_cmd_execution_mocked(self):
        """pdfcmd correctly resolves command and passes args to Process."""
        with mock.patch("proteindf_tools.Process") as mock_process:
            mock_proc_inst = mock.MagicMock()
            mock_process.return_value = mock_proc_inst
            mock_proc_inst.cmd.return_value = mock_proc_inst
            mock_proc_inst.commit.return_value = 0

            sys.argv = ["pdfcmd.py", "testrun", "--arg1", "val1"]
            rc = self.mod.main()

            self.assertEqual(rc, 0)
            mock_proc_inst.cmd.assert_called_once_with(f"{self.subcmd_path} --arg1 val1")

    def test_output_option(self):
        """--OUTPUT option specifies the logfile path for Process."""
        with mock.patch("proteindf_tools.Process") as mock_process:
            mock_proc_inst = mock.MagicMock()
            mock_process.return_value = mock_proc_inst
            mock_proc_inst.cmd.return_value = mock_proc_inst
            mock_proc_inst.commit.return_value = 0

            sys.argv = ["pdfcmd.py", "testrun", "--OUTPUT", "custom.log"]
            rc = self.mod.main()

            self.assertEqual(rc, 0)
            mock_process.assert_called_once_with(logfile_path="custom.log")

    @unittest.expectedFailure
    def test_nonexistent_command(self):
        """pdfcmd exits gracefully with non-zero exit code on nonexistent command.

        scripts/pdfcmd.py:70-73: コマンドが見つからない場合(yn == False)に
        処理を中断しないため、存在しないコマンドを実行しようとして
        FileNotFoundError などの例外で異常終了する。
        修正されたら外す。
        """
        sys.argv = ["pdfcmd.py", "nonexistent_subcmd"]
        buf = io.StringIO()
        with redirect_stderr(buf):
            rc = self.mod.main()
        self.assertNotEqual(rc, 0)
        self.assertIn("cannot command:", buf.getvalue())


# =============================================================================
# 7. pdf-env.py
# =============================================================================
class TestPdfEnv(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-env.py")

    def test_env_output(self):
        """pdf-env.py prints environment variables to stdout."""
        os.environ["PDF_TEST_ENV_VAR"] = "test_value_12345"
        sys.argv = ["pdf-env.py"]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        output = buf.getvalue()
        self.assertIn("PDF_TEST_ENV_VAR : test_value_12345", output)

    def test_env_options(self):
        """pdf-env.py runs with -v and -d options without error."""
        sys.argv = ["pdf-env.py", "-v", "-d"]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        self.assertTrue(len(buf.getvalue()) > 0)


# =============================================================================
# 8. pdf-show-profile.py
# =============================================================================
class TestPdfShowProfile(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-show-profile.py")
        self.stat_file = os.path.join(self.tmpdir.name, "sample.stat")
        cProfile.run("sum(range(50))", self.stat_file)

    def test_show_profile_success(self):
        """pdf-show-profile.py reads stats file and prints profiling information."""
        sys.argv = ["pdf-show-profile.py", self.stat_file]
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.mod.main()

        output = buf.getvalue()
        self.assertIn("Ordered by: cumulative time", output)
        self.assertIn("Ordered by: internal time", output)

    def test_nonexistent_stats_file(self):
        """pdf-show-profile.py raises FileNotFoundError for nonexistent stats file."""
        sys.argv = ["pdf-show-profile.py", "nonexistent.stat"]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()


# =============================================================================
# 9. pdf-plot-mo-tracer.py
# =============================================================================
class TestPdfPlotMoTracer(BaseScriptTestCase):
    def setUp(self):
        super().setUp()
        self.mod = _load_script("pdf-plot-mo-tracer.py")
        self.dir1 = os.path.join(self.tmpdir.name, "calc1")
        self.dir2 = os.path.join(self.tmpdir.name, "calc2")
        os.makedirs(self.dir1)
        os.makedirs(self.dir2)

        _create_sample_h5(os.path.join(self.dir1, "pdfresults.h5"), iterations=1, te=-75.0)
        _create_sample_h5(os.path.join(self.dir2, "pdfresults.h5"), iterations=1, te=-75.2)

        p1 = _create_sample_param(iterations=1, te=-75.0)
        p2 = _create_sample_param(iterations=1, te=-75.2)
        pdf.save_pdfparam(p1, os.path.join(self.dir1, "pdfparam.mpac"))
        pdf.save_pdfparam(p2, os.path.join(self.dir2, "pdfparam.mpac"))

    def _mock_run_pdf(self, args):
        """Mock info-mo-overlap call by writing a synthetic CSC matrix."""
        out_idx = args.index("-o") + 1
        target_csc = args[out_idx]
        m = pdf.Matrix(3, 3)
        m.set(0, 0, 0.95)
        m.set(1, 1, 0.90)
        m.set(2, 2, 0.85)
        m.save(target_csc)

    def test_plot_mo_tracer_success(self):
        """pdf-plot-mo-tracer.py mocks info-mo-overlap and creates mo-tracer PNG."""
        out_png = os.path.join(self.tmpdir.name, "tracer.png")
        csc_mat = os.path.join(self.tmpdir.name, "CSC.mat")
        sys.argv = [
            "pdf-plot-mo-tracer.py",
            self.dir1,
            self.dir2,
            "-m",
            csc_mat,
            "-o",
            out_png,
        ]
        with mock.patch("proteindf_tools.run_pdf", side_effect=self._mock_run_pdf):
            buf = io.StringIO()
            with redirect_stdout(buf):
                self.mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(os.path.exists(csc_mat))

    def test_output_option(self):
        """-o option specifies output image filename."""
        custom_png = os.path.join(self.tmpdir.name, "custom_mo.png")
        sys.argv = [
            "pdf-plot-mo-tracer.py",
            self.dir1,
            self.dir2,
            "-o",
            custom_png,
        ]
        with mock.patch("proteindf_tools.run_pdf", side_effect=self._mock_run_pdf):
            buf = io.StringIO()
            with redirect_stdout(buf):
                self.mod.main()

        self.assertTrue(os.path.exists(custom_png))

    def test_info_mo_overlap_missing(self):
        """pdf-plot-mo-tracer.py stops when info-mo-overlap executable is not found."""
        sys.argv = ["pdf-plot-mo-tracer.py", self.dir1, self.dir2]
        with mock.patch(
            "proteindf_tools.run_pdf",
            side_effect=FileNotFoundError("info-mo-overlap not found"),
        ):
            with self.assertRaises(FileNotFoundError):
                self.mod.main()


if __name__ == "__main__":
    unittest.main()
