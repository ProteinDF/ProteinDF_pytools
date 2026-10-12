#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import contextlib
import io
import tempfile
import unittest
from unittest.mock import MagicMock

try:
    import h5py
    import numpy as np

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

import matplotlib.pyplot as plt

from proteindf_tools.pdfparam_hdf5 import PdfParam_H5
from proteindf_tools.pdfparam_object import PdfParamObject
from proteindf_tools.report import PdfReport
from proteindf_tools.vector import Vector

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def is_valid_png(filepath):
    if not os.path.isfile(filepath):
        return False
    if os.path.getsize(filepath) == 0:
        return False
    with open(filepath, "rb") as f:
        sig = f.read(8)
    return sig == PNG_SIGNATURE


@unittest.skipUnless(HAVE_H5PY, "h5py or numpy not available")
class TestPdfReport(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = os.path.join(self.tmpdir.name, "report")
        self.h5_path = os.path.join(self.tmpdir.name, "test.h5")

    def tearDown(self):
        plt.close("all")
        self.tmpdir.cleanup()

    def _create_h5_rks(self, iterations=2, num_of_AOs=4, num_of_MOs=4):
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.num_of_atoms = 1
        p.num_of_AOs = num_of_AOs
        p.num_of_MOs = num_of_MOs
        p.method = "rks"
        p.iterations = iterations
        p.save_model()
        if iterations > 0:
            p.save_TEs(Vector([-100.0, -101.5][:iterations]))
            p.save_occ("rks", Vector([2.0, 2.0, 0.0, 0.0]))
            p.save_energy_level("rks", 1, Vector([-1.0, -0.5, 0.2, 0.8]))
            if iterations >= 2:
                p.save_energy_level("rks", 2, Vector([-1.1, -0.6, 0.1, 0.7]))
        p._close()

        p_loaded = PdfParam_H5()
        p_loaded.open(self.h5_path)
        return p_loaded

    def _create_h5_uks(self, iterations=2, num_of_AOs=4, num_of_MOs=4):
        p = PdfParam_H5()
        p.open(self.h5_path)
        p.num_of_atoms = 1
        p.num_of_AOs = num_of_AOs
        p.num_of_MOs = num_of_MOs
        p.method = "uks"
        p.iterations = iterations
        p.save_model()
        if iterations > 0:
            p.save_TEs(Vector([-100.0, -101.5][:iterations]))
            for rt in ["uks_alpha", "uks_beta"]:
                p.save_occ(rt, Vector([1.0, 1.0, 0.0, 0.0]))
                p.save_energy_level(rt, 1, Vector([-1.0, -0.5, 0.2, 0.8]))
                if iterations >= 2:
                    p.save_energy_level(rt, 2, Vector([-1.1, -0.6, 0.1, 0.7]))
        p._close()

        p_loaded = PdfParam_H5()
        p_loaded.open(self.h5_path)
        return p_loaded

    # -------------------------------------------------------------------------
    # 1. PdfReport instantiation & _make_work_dir
    # -------------------------------------------------------------------------
    def test_init_invalid_type_raises_assertion_error(self):
        """PdfReport(not PdfParamObject) raises AssertionError."""
        with self.assertRaises(AssertionError):
            PdfReport("invalid_param")

        with self.assertRaises(AssertionError):
            PdfReport(None)

    def test_make_work_dir(self):
        """_make_work_dir creates workdir if not existing, or succeeds if existing."""
        p = MagicMock(spec=PdfParamObject)
        rep = PdfReport(p, self.workdir)

        self.assertFalse(os.path.exists(self.workdir))
        rep._make_work_dir()
        self.assertTrue(os.path.isdir(self.workdir))

        # calling again when directory already exists should not raise error
        rep._make_work_dir()
        self.assertTrue(os.path.isdir(self.workdir))

    # -------------------------------------------------------------------------
    # 2. _get_total_energy_vector & _get_rst
    # -------------------------------------------------------------------------
    def test_get_total_energy_vector(self):
        """_get_total_energy_vector returns list of TEs for all iterations."""
        p = self._create_h5_rks(iterations=2)
        try:
            rep = PdfReport(p, self.workdir)
            tes = rep._get_total_energy_vector()
            self.assertEqual(len(tes), 2)
            self.assertAlmostEqual(tes[0], -100.0)
            self.assertAlmostEqual(tes[1], -101.5)
        finally:
            p._close()

    def test_get_rst(self):
        """_get_rst generates reStructuredText table containing AOs, MOs, iterations, TE."""
        p = self._create_h5_rks(iterations=2, num_of_AOs=4, num_of_MOs=4)
        try:
            rep = PdfReport(p, self.workdir)
            rst = rep._get_rst()
            self.assertIn("ProteinDF report:", rst)
            self.assertRegex(rst, r"the number of AOs\s+4")
            self.assertRegex(rst, r"the number of MOs\s+4")
            self.assertRegex(rst, r"the number of iterations\s+2")
            self.assertRegex(rst, r"Total energy\s+-101\.50000000")
        finally:
            p._close()

    def test_get_rst_iterations_zero(self):
        """_get_rst raises IndexError when iterations == 0 due to TEs[-1] access."""
        p = self._create_h5_rks(iterations=0)
        try:
            rep = PdfReport(p, self.workdir)
            tes = rep._get_total_energy_vector()
            self.assertEqual(tes, [])
            with self.assertRaises(IndexError):
                rep._get_rst()
        finally:
            p._close()

    # -------------------------------------------------------------------------
    # 3. _plot_convergence_check & _plot_convergence_energy_level
    # -------------------------------------------------------------------------
    def test_plot_convergence_check(self):
        """_plot_convergence_check creates TE_hist.dat and TE_hist.png."""
        p = self._create_h5_rks(iterations=2)
        try:
            rep = PdfReport(p, self.workdir)
            rep._plot_convergence_check(2)

            dat_path = os.path.join(self.workdir, "TE_hist.dat")
            png_path = os.path.join(self.workdir, "TE_hist.png")

            self.assertTrue(os.path.isfile(dat_path))
            self.assertTrue(is_valid_png(png_path))

            with open(dat_path, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
            self.assertEqual(len(lines), 2)
            self.assertIn("1,", lines[0])
            self.assertIn("-100.0000000000", lines[0])
            self.assertIn("2,", lines[1])
            self.assertIn("-101.5000000000", lines[1])
        finally:
            p._close()

    def test_plot_convergence_energy_level(self):
        """_plot_convergence_energy_level creates eigvals_hist.<rt>.dat/.png and eigvals_last.<rt>.png."""
        p = self._create_h5_rks(iterations=2)
        try:
            rep = PdfReport(p, self.workdir)
            rep._plot_convergence_energy_level("rks", 2)

            dat_path = os.path.join(self.workdir, "eigvals_hist.rks.dat")
            png_hist_path = os.path.join(self.workdir, "eigvals_hist.rks.png")
            png_last_path = os.path.join(self.workdir, "eigvals_last.rks.png")

            self.assertTrue(os.path.isfile(dat_path))
            self.assertTrue(is_valid_png(png_hist_path))
            self.assertTrue(is_valid_png(png_last_path))

            with open(dat_path, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
            # 2 iterations * 4 levels = 8 lines
            self.assertEqual(len(lines), 8)
            # Energy in dat is in eV (multiplied by 27.2116)
            # itr 1, level 0: -1.0 * 27.2116 = -27.2116000000
            parts = [p.strip() for p in lines[0].split(",")]
            self.assertEqual(parts[0], "1")
            self.assertEqual(parts[1], "0")
            self.assertAlmostEqual(float(parts[2]), -27.2116, places=4)
        finally:
            p._close()

    # -------------------------------------------------------------------------
    # 4. report(iteration=0) & method="uks"
    # -------------------------------------------------------------------------
    def test_report_rks(self):
        """report(iteration=0) runs full reporting pipeline for RKS."""
        p = self._create_h5_rks(iterations=2)
        try:
            rep = PdfReport(p, self.workdir)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rep.report(iteration=0)
            output = buf.getvalue()
            self.assertIn("ProteinDF report:", output)

            self.assertTrue(os.path.isfile(os.path.join(self.workdir, "TE_hist.dat")))
            self.assertTrue(is_valid_png(os.path.join(self.workdir, "TE_hist.png")))
            self.assertTrue(os.path.isfile(os.path.join(self.workdir, "eigvals_hist.rks.dat")))
            self.assertTrue(is_valid_png(os.path.join(self.workdir, "eigvals_hist.rks.png")))
            self.assertTrue(is_valid_png(os.path.join(self.workdir, "eigvals_last.rks.png")))
        finally:
            p._close()

    def test_report_uks(self):
        """report(iteration=0) generates outputs for both run_types when method is UKS."""
        p = self._create_h5_uks(iterations=2)
        try:
            rep = PdfReport(p, self.workdir)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rep.report(iteration=0)
            output = buf.getvalue()
            self.assertIn("ProteinDF report:", output)

            self.assertTrue(os.path.isfile(os.path.join(self.workdir, "TE_hist.dat")))
            self.assertTrue(is_valid_png(os.path.join(self.workdir, "TE_hist.png")))

            for rt in ["uks_alpha", "uks_beta"]:
                self.assertTrue(os.path.isfile(os.path.join(self.workdir, f"eigvals_hist.{rt}.dat")))
                self.assertTrue(is_valid_png(os.path.join(self.workdir, f"eigvals_hist.{rt}.png")))
                self.assertTrue(is_valid_png(os.path.join(self.workdir, f"eigvals_last.{rt}.png")))
        finally:
            p._close()

    # -------------------------------------------------------------------------
    # 5. _make_plot_pop_mulliken_data & _plot_pop_mulliken
    # -------------------------------------------------------------------------
    def test_make_plot_pop_mulliken_data(self):
        """_make_plot_pop_mulliken_data writes population vector to output file."""
        p = MagicMock(spec=PdfParamObject)
        rep = PdfReport(p, self.workdir)
        rep._make_work_dir()
        out_path = os.path.join(self.workdir, "mulliken.dat")

        pop_vtr = [0.15, -0.25, 0.10]
        rep._make_plot_pop_mulliken_data(pop_vtr, out_path)

        self.assertTrue(os.path.isfile(out_path))
        with open(out_path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
        self.assertEqual(lines, ["1, 0.15", "2, -0.25", "3, 0.1"])

    @unittest.expectedFailure
    def test_plot_pop_mulliken_with_pdfparam_h5(self):
        """_plot_pop_mulliken generates mulliken.dat and mulliken.png.

        Expected failure: PdfParam_H5.get_pop_mulliken_atom calls nonexistent
        self._load_vector, raising AttributeError: 'PdfParam_H5' object has no attribute '_load_vector'.
        """
        p = self._create_h5_rks(iterations=1)
        try:
            rep = PdfReport(p, self.workdir)
            rep._plot_pop_mulliken(1)
            self.assertTrue(os.path.isfile(os.path.join(self.workdir, "mulliken.dat")))
            self.assertTrue(is_valid_png(os.path.join(self.workdir, "mulliken.png")))
        finally:
            p._close()

    @unittest.expectedFailure
    def test_plot_pop_mulliken_uks_run_type(self):
        """_plot_pop_mulliken uses appropriate run_type for UKS method.

        Expected failure: _plot_pop_mulliken hardcodes run_type = 'rks' regardless
        of self._pdfparam.method.
        """
        p = MagicMock(spec=PdfParamObject)
        p.method = "uks"
        p.get_pop_mulliken_atom.return_value = Vector([0.1, -0.1])
        rep = PdfReport(p, self.workdir)
        rep._make_work_dir()
        rep._plot_pop_mulliken(1)

        # Should query for UKS run_type (e.g. uks_alpha / uks_beta), not hardcoded "rks"
        call_args = p.get_pop_mulliken_atom.call_args[0]
        self.assertIn(call_args[0], ["uks_alpha", "uks_beta"])


if __name__ == "__main__":
    unittest.main()
