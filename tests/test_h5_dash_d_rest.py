#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import io
import logging
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5
from proteindf_tools.vector import Vector

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _is_valid_png(filepath):
    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return False
    with open(filepath, "rb") as f:
        return f.read(8) == PNG_SIGNATURE


def _create_synthetic_system(tmpdir, iterations=2, te=-75.0):
    """Create identical synthetic H2O + X systems in .db and .h5 formats."""
    mol = bridge.AtomGroup()
    a_o = bridge.Atom(
        symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="O1"
    )
    a_h1 = bridge.Atom(
        symbol="H", position=bridge.Position([0.0, 1.43, -1.1]), charge=1.0, name="H1"
    )
    a_h2 = bridge.Atom(
        symbol="H", position=bridge.Position([0.0, -1.43, -1.1]), charge=1.0, name="H2"
    )
    a_x = bridge.Atom(
        symbol="X", position=bridge.Position([0.0, 0.0, 3.0]), charge=0.5, name="X1"
    )
    mol.set_atom("0", a_o)
    mol.set_atom("1", a_h1)
    mol.set_atom("2", a_h2)
    mol.set_atom("3", a_x)

    bs_o = BasisSet("sto-3g-O", 1)
    cg_o = ContractedGTO("s", 1)
    cg_o[0] = PrimitiveGTO(130.7, 0.15)
    bs_o[0] = cg_o

    bs_h = BasisSet("sto-3g-H", 1)
    cg_h = ContractedGTO("s", 1)
    cg_h[0] = PrimitiveGTO(3.4, 0.15)
    bs_h[0] = cg_h

    occ_vec = Vector([2.0, 2.0, 0.0, 0.0])
    elevels = [
        Vector([-19.1 - 0.01 * itr, -1.2 - 0.01 * itr, 0.3 + 0.01 * itr, 0.9 + 0.01 * itr])
        for itr in range(1, iterations + 1)
    ]
    tes = [te - 0.1 * i for i in range(iterations)]

    # SQLite archive
    param = PdfParam()
    param.comment = "synthetic H2O+X system"
    param.num_of_atoms = 4
    param.num_of_AOs = 4
    param.num_of_MOs = 4
    param.method = "rks"
    param.iterations = iterations
    param.guess = "none"
    param._set_scf_converged(True)
    param.molecule = mol
    param.set_basisset("O@O1", bs_o)
    param.set_basisset("H@H1", bs_h)
    param.set_basisset("H@H2", bs_h)
    for itr in range(1, iterations + 1):
        param.set_total_energy(itr, tes[itr - 1])

    db_path = os.path.join(tmpdir, "water.db")
    db = PdfArchive_Sqlite3(db_path, pdfparam=param)
    db.set_occupations("rks", occ_vec)
    for itr in range(1, iterations + 1):
        db.set_energylevel("rks", itr, elevels[itr - 1])

    # HDF5 archive
    h5_path = os.path.join(tmpdir, "water.h5")
    p5 = PdfParam_H5()
    p5.open(h5_path)
    p5.comment = "synthetic H2O+X system"
    p5.num_of_atoms = 4
    p5.num_of_AOs = 4
    p5.num_of_MOs = 4
    p5.method = "rks"
    p5.iterations = iterations
    p5.guess = "none"
    p5._set_scf_converged(True)
    p5.molecule = mol
    p5.set_basisset("O@O1", bs_o)
    p5.set_basisset("H@H1", bs_h)
    p5.set_basisset("H@H2", bs_h)
    p5.save_model()
    p5.save_TEs(Vector(tes))
    p5.save_occ("rks", occ_vec)
    for itr in range(1, iterations + 1):
        p5.save_energy_level("rks", itr, elevels[itr - 1])
    p5._close()

    return {"db": db_path, "h5": h5_path}


class TestH5DashDRest(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        root_logger = logging.getLogger()
        for h in list(root_logger.handlers):
            h.close()
            root_logger.removeHandler(h)
        self.data = _create_synthetic_system(self.tmpdir.name)

    def tearDown(self):
        root_logger = logging.getLogger()
        for h in list(root_logger.handlers):
            h.close()
            root_logger.removeHandler(h)
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        plt.close("all")
        self.tmpdir.cleanup()

    # =========================================================================
    # 1. pdf-plot-elevel.py
    # =========================================================================
    def test_plot_elevel_h5_matches_db(self):
        """pdf-plot-elevel.py generates identical dat and valid PNG for .h5 and .db."""
        mod = _load_script("pdf-plot-elevel.py")
        png_db = os.path.join(self.tmpdir.name, "elevel_db.png")
        png_h5 = os.path.join(self.tmpdir.name, "elevel_h5.png")

        # Run with .db
        sys.argv = ["pdf-plot-elevel.py", "-d", self.data["db"], "-o", png_db, "-v"]
        out_db = io.StringIO()
        with redirect_stdout(out_db):
            mod.main()
        with open("tmp.dat", "r") as f:
            dat_db = f.read()

        # Run with .h5
        sys.argv = ["pdf-plot-elevel.py", "-d", self.data["h5"], "-o", png_h5, "-v"]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()
        with open("tmp.dat", "r") as f:
            dat_h5 = f.read()

        self.assertEqual(dat_db, dat_h5)
        self.assertTrue(_is_valid_png(png_db))
        self.assertTrue(_is_valid_png(png_h5))
        self.assertIn("HOMO: 1", out_db.getvalue())
        self.assertIn("HOMO: 1", out_h5.getvalue())

    def test_plot_elevel_default_h5(self):
        """pdf-plot-elevel.py -d without argument defaults to pdfresults.h5."""
        shutil.copyfile(self.data["h5"], os.path.join(self.tmpdir.name, "pdfresults.h5"))
        mod = _load_script("pdf-plot-elevel.py")
        out_png = os.path.join(self.tmpdir.name, "elevel_default.png")

        sys.argv = ["pdf-plot-elevel.py", "-d", "-o", out_png]
        mod.main()
        self.assertTrue(_is_valid_png(out_png))

    def test_plot_elevel_h5_multiple_runtypes_raises(self):
        """pdf-plot-elevel.py raises ValueError on HDF5 with multiple runtypes (UKS)."""
        h5_uks = os.path.join(self.tmpdir.name, "uks.h5")
        p = PdfParam_H5()
        p.open(h5_uks)
        p.method = "uks"
        p.iterations = 1
        p.num_of_atoms = 1
        p.num_of_AOs = 2
        p.num_of_MOs = 2
        p.save_model()
        p.save_TEs(Vector([-1.0]))
        for rt in ("uks_alpha", "uks_beta"):
            p.save_occ(rt, Vector([1.0, 0.0]))
            p.save_energy_level(rt, 1, Vector([-0.5, 0.2]))
        p._close()

        mod = _load_script("pdf-plot-elevel.py")
        sys.argv = ["pdf-plot-elevel.py", "-d", h5_uks]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        msg = str(cm.exception)
        self.assertIn("uks_alpha", msg)
        self.assertIn("uks_beta", msg)

    def test_plot_elevel_h5_zero_runtypes_raises(self):
        """pdf-plot-elevel.py raises ValueError on HDF5 without energy_level."""
        h5_empty = os.path.join(self.tmpdir.name, "empty.h5")
        p = PdfParam_H5()
        p.open(h5_empty)
        p.method = "rks"
        p.iterations = 1
        p.num_of_atoms = 1
        p.num_of_AOs = 2
        p.num_of_MOs = 2
        p.save_model()
        p.save_TEs(Vector([-1.0]))
        p._close()

        mod = _load_script("pdf-plot-elevel.py")
        sys.argv = ["pdf-plot-elevel.py", "-d", h5_empty]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("ENERGY_LEVEL", str(cm.exception))

    # =========================================================================
    # 2. pdf-report.py
    # =========================================================================
    def test_report_h5_matches_db(self):
        """pdf-report.py produces identical output files and stdout for .h5 and .db."""
        mod = _load_script("pdf-report.py")
        dir_db = os.path.join(self.tmpdir.name, "report_db")
        dir_h5 = os.path.join(self.tmpdir.name, "report_h5")

        sys.argv = ["pdf-report.py", self.data["db"], dir_db]
        out_db = io.StringIO()
        with redirect_stdout(out_db):
            mod.main()

        sys.argv = ["pdf-report.py", self.data["h5"], dir_h5]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertEqual(out_db.getvalue(), out_h5.getvalue())
        self.assertIn("ProteinDF report:", out_h5.getvalue())

        for dat_file in ("TE_hist.dat", "eigvals_hist.dat"):
            with open(os.path.join(dir_db, dat_file)) as f1, open(
                os.path.join(dir_h5, dat_file)
            ) as f2:
                self.assertEqual(f1.read(), f2.read())

        for png_file in ("TE_hist.png", "eigvals_hist.png", "eigvals_last.png"):
            self.assertTrue(_is_valid_png(os.path.join(dir_db, png_file)))
            self.assertTrue(_is_valid_png(os.path.join(dir_h5, png_file)))

    def test_report_assert_isinstance_accepts_both(self):
        """plot_convergence_check and get_rst accept both PdfArchive and PdfParam_H5."""
        mod = _load_script("pdf-report.py")
        db_entry = pdf.open_archive(self.data["db"])
        h5_entry = pdf.open_archive(self.data["h5"])
        try:
            self.assertEqual(mod.get_rst(db_entry), mod.get_rst(h5_entry))
            with self.assertRaises(AssertionError):
                mod.get_rst("invalid")
            with self.assertRaises(AssertionError):
                mod.plot_convergence_check("invalid", self.tmpdir.name)
        finally:
            h5_entry._close()

    # =========================================================================
    # 3. pdf-test.py
    # =========================================================================
    def test_pdf_test_h5_and_db_comparisons(self):
        """pdf-test.py compares .h5/.h5, .db/.db, and .h5/.db identically."""
        mod = _load_script("pdf-test.py")

        diff_dir = os.path.join(self.tmpdir.name, "diff")
        os.makedirs(diff_dir)
        data_same = _create_synthetic_system(diff_dir, te=-75.0)
        diff_dir2 = os.path.join(self.tmpdir.name, "diff2")
        os.makedirs(diff_dir2)
        data_diff = _create_synthetic_system(diff_dir2, te=-76.0)

        # Same energy: .db vs .db, .h5 vs .h5, .h5 vs .db -> all exit 0
        for p1, p2 in [
            (self.data["db"], data_same["db"]),
            (self.data["h5"], data_same["h5"]),
            (self.data["h5"], data_same["db"]),
        ]:
            sys.argv = ["pdf-test.py", p1, p2]
            buf = io.StringIO()
            with redirect_stdout(buf), redirect_stderr(buf):
                with self.assertRaises(SystemExit) as cm:
                    mod.main()
            self.assertEqual(cm.exception.code, 0, msg=f"Failed for {p1} vs {p2}")

        # Different energy: .db vs .db, .h5 vs .h5, .h5 vs .db -> all exit 1
        for p1, p2 in [
            (self.data["db"], data_diff["db"]),
            (self.data["h5"], data_diff["h5"]),
            (self.data["h5"], data_diff["db"]),
        ]:
            sys.argv = ["pdf-test.py", p1, p2]
            buf = io.StringIO()
            with redirect_stdout(buf), redirect_stderr(buf):
                with self.assertRaises(SystemExit) as cm:
                    mod.main()
            self.assertEqual(cm.exception.code, 1, msg=f"Failed for {p1} vs {p2}")

        # -d (--debug) option preserves its debug flag meaning
        sys.argv = ["pdf-test.py", "-d", self.data["h5"], data_same["h5"]]
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            with self.assertRaises(SystemExit) as cm:
                mod.main()
        self.assertEqual(cm.exception.code, 0)

    # =========================================================================
    # 4. pdf-info-iteration.py, pdf-plot-dos.py, pdf-plot-gvalues.py defaults
    # =========================================================================
    def test_info_iteration_default_h5(self):
        """pdf-info-iteration.py -d without path uses pdfresults.h5 by default."""
        shutil.copyfile(self.data["h5"], os.path.join(self.tmpdir.name, "pdfresults.h5"))
        mod = _load_script("pdf-info-iteration.py")
        sys.argv = ["pdf-info-iteration.py", "-d"]
        out = io.StringIO()
        with redirect_stdout(out):
            mod.main()
        self.assertEqual(out.getvalue().strip(), "2")

    def test_plot_dos_default_h5(self):
        """pdf-plot-dos.py without -d uses pdfresults.h5 by default."""
        shutil.copyfile(self.data["h5"], os.path.join(self.tmpdir.name, "pdfresults.h5"))
        mod = _load_script("pdf-plot-dos.py")
        out_png = os.path.join(self.tmpdir.name, "dos_default.png")
        sys.argv = ["pdf-plot-dos.py", "-o", out_png]
        out = io.StringIO()
        with redirect_stdout(out):
            mod.main()
        self.assertTrue(_is_valid_png(out_png))

    def test_plot_gvalues_default_h5(self):
        """pdf-plot-gvalues.py without -d uses pdfresults.h5 by default."""
        shutil.copyfile(self.data["h5"], os.path.join(self.tmpdir.name, "pdfresults.h5"))
        mod = _load_script("pdf-plot-gvalues.py")
        out_png = os.path.join(self.tmpdir.name, "gvalues_default.png")

        def fake_run_pdf(args):
            vtr = Vector([0.1, 0.2, 0.3, 0.4])
            vtr.save("g-values.vtr")

        sys.argv = ["pdf-plot-gvalues.py", "-o", out_png]
        out = io.StringIO()
        with redirect_stdout(out), mock.patch("proteindf_tools.run_pdf", side_effect=fake_run_pdf):
            mod.main()
        self.assertTrue(_is_valid_png(out_png))


if __name__ == "__main__":
    unittest.main()
