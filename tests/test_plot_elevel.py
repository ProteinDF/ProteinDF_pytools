#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt

import proteindf_bridge as bridge
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
import proteindf_tools.pdfgraph as pg
import proteindf_tools.vector as vec

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "pdf-plot-elevel.py"
spec = importlib.util.spec_from_file_location("pdf_plot_elevel", SCRIPT_PATH)
pdf_plot_elevel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pdf_plot_elevel)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


class TestPdfPlotElevel(unittest.TestCase):
    def setUp(self):
        self._temp_dir = tempfile.TemporaryDirectory()
        self.temp_dir = self._temp_dir.name

    def tearDown(self):
        plt.close("all")
        self._temp_dir.cleanup()

    def test_plot_elevel_e_mode(self):
        v = vec.Vector(5)
        for i in range(5):
            v.set(i, -10.0 + i)
        vec_path = os.path.join(self.temp_dir, "test.vec")
        v.save(vec_path)

        out_png = os.path.join(self.temp_dir, "out.png")

        orig_cwd = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            with mock.patch.object(
                sys,
                "argv",
                [
                    "pdf-plot-elevel.py",
                    "-e",
                    vec_path,
                    "-o",
                    out_png,
                    "-l",
                    "2",
                ],
            ):
                pdf_plot_elevel.main()
        finally:
            os.chdir(orig_cwd)

        self.assertTrue(os.path.isfile(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)
        with open(out_png, "rb") as f:
            self.assertEqual(f.read(8), PNG_SIGNATURE)

    def test_df_matrix_graph_no_stdout(self):
        m = bridge.Matrix(2, 2)
        m.set(0, 0, 1.0)
        m.set(1, 1, 2.0)
        g = pg.DfMatrixGraph(m)
        g.should_write_values = True
        out_png = os.path.join(self.temp_dir, "matrix.png")

        buf = io.StringIO()
        with mock.patch("sys.stdout", buf):
            g.save(out_png)

        self.assertEqual(buf.getvalue(), "")

    def _create_synth_h5(self, h5_path, iterations=3, runtypes=("rks",), add_occ=True):
        from proteindf_tools.pdfparam_hdf5 import PdfParam_H5

        p = PdfParam_H5()
        p.open(h5_path)
        p.iterations = iterations
        p.num_of_atoms = 1
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(vec.Vector([-1.0] * iterations))
        for rt in runtypes:
            elevel = vec.Vector(4)
            for i in range(4):
                elevel.set(i, -10.0 + i)
            p.save_energy_level(rt, iterations, elevel)
            if add_occ:
                occ = vec.Vector(4)
                occ.set(0, 2.0)
                occ.set(1, 2.0)
                occ.set(2, 0.0)
                occ.set(3, 0.0)
                p.save_occ(rt, occ)
        p._close()

    def test_plot_elevel_d_mode_single_runtype(self):
        h5_path = os.path.join(self.temp_dir, "synth.h5")
        self._create_synth_h5(h5_path, iterations=3, runtypes=["rks"])
        out_png = os.path.join(self.temp_dir, "out_d.png")

        orig_cwd = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            with mock.patch.object(
                sys,
                "argv",
                ["pdf-plot-elevel.py", "-d", h5_path, "-o", out_png],
            ):
                pdf_plot_elevel.main()
        finally:
            os.chdir(orig_cwd)

        self.assertTrue(os.path.isfile(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)
        with open(out_png, "rb") as f:
            self.assertEqual(f.read(8), PNG_SIGNATURE)

    def test_plot_elevel_d_mode_multiple_runtypes(self):
        h5_path = os.path.join(self.temp_dir, "multi.h5")
        self._create_synth_h5(h5_path, iterations=3, runtypes=["uks_alpha", "uks_beta"])
        out_png = os.path.join(self.temp_dir, "out_multi.png")

        orig_cwd = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            with mock.patch.object(
                sys,
                "argv",
                ["pdf-plot-elevel.py", "-d", h5_path, "-o", out_png],
            ):
                with self.assertRaises(ValueError) as cm:
                    pdf_plot_elevel.main()
        finally:
            os.chdir(orig_cwd)

        msg = str(cm.exception)
        self.assertIn("uks_alpha", msg)
        self.assertIn("uks_beta", msg)
        self.assertTrue("UKS" in msg or "ROKS" in msg)

    def test_plot_elevel_d_mode_zero_runtypes(self):
        h5_path = os.path.join(self.temp_dir, "empty.h5")
        self._create_synth_h5(h5_path, iterations=3, runtypes=[])
        out_png = os.path.join(self.temp_dir, "out_zero.png")

        orig_cwd = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            with mock.patch.object(
                sys,
                "argv",
                ["pdf-plot-elevel.py", "-d", h5_path, "-o", out_png],
            ):
                with self.assertRaises(ValueError) as cm:
                    pdf_plot_elevel.main()
        finally:
            os.chdir(orig_cwd)

        msg = str(cm.exception)
        self.assertIn("ENERGY_LEVEL", msg)
        self.assertNotIn("UKS", msg)
        self.assertNotIn("ROKS", msg)

    def test_plot_elevel_d_mode_db_raises_error(self):
        db_path = os.path.join(self.temp_dir, "test.db")
        out_png = os.path.join(self.temp_dir, "out_db.png")

        with mock.patch.object(
            sys,
            "argv",
            ["pdf-plot-elevel.py", "-d", db_path, "-o", out_png],
        ):
            with self.assertRaises(ValueError) as cm:
                pdf_plot_elevel.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

    def test_plot_elevel_single_vertical_range(self):
        data_path = os.path.join(self.temp_dir, "test.dat")
        with open(data_path, "w") as f:
            f.write("1, 0, -10.0\n")
        out_png = os.path.join(self.temp_dir, "range.png")

        with mock.patch.object(pdf_plot_elevel.pdf, "DfEnergyLevelHistoryGraphV") as mock_graph_cls:
            mock_inst = mock_graph_cls.return_value

            # Iteration 1: matches existing 0.4 - 1.6
            pdf_plot_elevel.plot_elevel_single(data_path, 1, 0, out_png)
            self.assertAlmostEqual(mock_inst.ymin, 0.4)
            self.assertAlmostEqual(mock_inst.ymax, 1.6)

            # Iteration N (N != 1): vertical range must include N
            pdf_plot_elevel.plot_elevel_single(data_path, 5, 0, out_png)
            self.assertAlmostEqual(mock_inst.ymin, 4.4)
            self.assertAlmostEqual(mock_inst.ymax, 5.6)
            self.assertTrue(mock_inst.ymin <= 5 <= mock_inst.ymax)

    def test_plot_elevel_d_mode_does_not_call_load(self):
        h5_path = os.path.join(self.temp_dir, "load_check.h5")
        self._create_synth_h5(h5_path, iterations=3, runtypes=["rks"])
        out_png = os.path.join(self.temp_dir, "out_load.png")

        orig_cwd = os.getcwd()
        try:
            os.chdir(self.temp_dir)
            with mock.patch.object(
                vec.Vector,
                "load",
                side_effect=AssertionError("eigvals.load() must not be called in -d mode"),
            ):
                with mock.patch.object(
                    sys,
                    "argv",
                    ["pdf-plot-elevel.py", "-d", h5_path, "-o", out_png],
                ):
                    pdf_plot_elevel.main()
        finally:
            os.chdir(orig_cwd)

        self.assertTrue(os.path.isfile(out_png))


if __name__ == "__main__":
    unittest.main()
