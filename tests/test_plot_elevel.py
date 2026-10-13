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


if __name__ == "__main__":
    unittest.main()
