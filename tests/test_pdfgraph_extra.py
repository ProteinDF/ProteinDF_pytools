#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import tempfile
import unittest

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy
import proteindf_bridge as bridge
import proteindf_tools.pdfgraph as pg


class TestPdfGraphExtra(unittest.TestCase):
    def setUp(self):
        self.tmpfiles = []

    def tearDown(self):
        plt.close("all")
        for f in self.tmpfiles:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    def _create_tempfile(self, content=None, suffix=".png"):
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False, mode="w" if isinstance(content, str) else "wb") as fp:
            if content is not None:
                fp.write(content)
            path = fp.name
        self.tmpfiles.append(path)
        return path

    def _assert_valid_png(self, path):
        self.assertTrue(os.path.exists(path), f"File {path} does not exist")
        self.assertGreater(os.path.getsize(path), 0, f"File {path} is empty")
        with open(path, "rb") as fp:
            header = fp.read(8)
        self.assertEqual(header, b"\x89PNG\r\n\x1a\n", f"File {path} is not a valid PNG")

    # ----------------------------------------------------------------------
    # 1. MidpointNormalize / MidpointLogNorm tests
    # ----------------------------------------------------------------------
    def test_midpoint_normalize(self):
        norm = pg.MidpointNormalize(vmin=-10.0, vmax=20.0, midpoint=0.0)
        self.assertAlmostEqual(float(norm(0.0)), 0.5)
        self.assertAlmostEqual(float(norm(-10.0)), 0.0)
        self.assertAlmostEqual(float(norm(20.0)), 1.0)
        self.assertAlmostEqual(float(norm(-5.0)), 0.25)
        self.assertAlmostEqual(float(norm(10.0)), 0.75)
        # Out of bounds: numpy.interp clamps to y[0] and y[-1]
        self.assertAlmostEqual(float(norm(-20.0)), 0.0)
        self.assertAlmostEqual(float(norm(30.0)), 1.0)

        # Array input
        arr = numpy.array([-10.0, -5.0, 0.0, 10.0, 20.0])
        res = norm(arr)
        numpy.testing.assert_allclose(res, [0.0, 0.25, 0.5, 0.75, 1.0])

    def test_midpoint_log_norm(self):
        norm = pg.MidpointLogNorm(vmin=1.0, vmax=100.0, midpoint=10.0)
        self.assertAlmostEqual(float(norm(10.0)), 0.5)
        self.assertAlmostEqual(float(norm(1.0)), 0.0)
        self.assertAlmostEqual(float(norm(100.0)), 1.0)
        # Out of bounds
        self.assertAlmostEqual(float(norm(0.1)), 0.0)
        self.assertAlmostEqual(float(norm(200.0)), 1.0)

    # ----------------------------------------------------------------------
    # 2. DfMatrixGraph / DfVectorGraph tests
    # ----------------------------------------------------------------------
    def test_df_matrix_graph_matrix(self):
        m = bridge.Matrix(2, 3)
        values = [[-1.0, 0.0, 1.0], [2.0, -2.0, 0.5]]
        for r in range(2):
            for c in range(3):
                m.set(r, c, values[r][c])

        g = pg.DfMatrixGraph(m)
        self.assertTrue(g.is_diverging)
        self.assertFalse(g.should_write_values)

        # Configuration: cmap, vmin/vmax, is_diverging, should_write_values
        g.cmap = cm.viridis
        self.assertEqual(g.cmap, cm.viridis)
        g.vmin = -3.0
        g.vmax = 3.0
        self.assertEqual(g.vmin, -3.0)
        self.assertEqual(g.vmax, 3.0)
        g.is_diverging = False
        self.assertFalse(g.is_diverging)
        g.should_write_values = True
        self.assertTrue(g.should_write_values)

        out_path = self._create_tempfile()
        g.save(out_path)
        self._assert_valid_png(out_path)

    def test_df_matrix_graph_symmetric_matrix(self):
        sm = bridge.SymmetricMatrix(2)
        sm.set(0, 0, 1.0)
        sm.set(0, 1, -2.0)
        sm.set(1, 1, 3.0)

        g = pg.DfMatrixGraph(sm)
        out_path = self._create_tempfile()
        g.save(out_path)
        self._assert_valid_png(out_path)

    def test_df_vector_graph(self):
        v = bridge.Vector(4)
        for i, val in enumerate([1.5, -2.5, 3.0, 0.0]):
            v.set(i, val)

        g = pg.DfVectorGraph(v)
        out_path = self._create_tempfile()
        g.save(out_path)
        self._assert_valid_png(out_path)

    # ----------------------------------------------------------------------
    # 3. DfDosGraph / DfDistanceVsElementGraph / DfEnergyLevelTraceGraph / DfGvalsGraph
    # ----------------------------------------------------------------------
    def test_df_dos_graph(self):
        csv_data = "0.0,1.0,2.0\n1.0,2.0,3.0\n2.0,3.0,4.0\n"
        csv_path = self._create_tempfile(content=csv_data, suffix=".csv")

        dos = pg.DfDosGraph()
        dos.load_data(csv_path)
        self.assertEqual(dos.num_of_series, 2)
        dos.add_draw_series(0, label="series_0", fmt="b-")
        dos.add_draw_series(1)

        out_path = self._create_tempfile()
        dos.save(out_path)
        self._assert_valid_png(out_path)

    def test_df_distance_vs_element_graph(self):
        csv_data = "1.0,2.0\n2.0,3.0\n3.0,1.0\n"
        csv_path = self._create_tempfile(content=csv_data, suffix=".csv")

        # Linear
        dve = pg.DfDistanceVsElementGraph(log=False)
        dve.load(csv_path)
        out_path = self._create_tempfile()
        dve.save(out_path)
        self._assert_valid_png(out_path)

        # Logarithmic
        dve_log = pg.DfDistanceVsElementGraph(log=True)
        dve_log.load(csv_path)
        out_path_log = self._create_tempfile()
        dve_log.save(out_path_log)
        self._assert_valid_png(out_path_log)

    def test_df_energy_level_trace_graph(self):
        elevel0 = [-15.0, -10.0, -5.0, 1.0]
        elevel1 = [-14.5, -9.8, -4.5, 1.5]
        connections = [
            [0, 0, 0.95],
            [1, 1, 0.9],
            [2, 2, 1.0],
            [3, 3, 0.8],
        ]

        elt = pg.DfEnergyLevelTraceGraph()
        elt.set_data(elevel0, elevel1, connections)
        elt.reference_line_color = "blue"
        elt.target_line_color = "gray"
        elt.connection_color = "red"
        self.assertEqual(elt.reference_line_color, "blue")
        self.assertEqual(elt.target_line_color, "gray")
        self.assertEqual(elt.connection_color, "red")

        out_path = self._create_tempfile()
        elt.save(out_path)
        self._assert_valid_png(out_path)

    def test_df_gvals_graph(self):
        x_occ = [-10.0, -8.0, -5.0]
        y_occ = [2.001, 2.002, 2.005]
        x_vir = [1.0, 3.0]
        y_vir = [2.010, 2.015]

        gv = pg.DfGvalsGraph()
        gv.set_data(x_occ, y_occ, x_vir, y_vir)

        out_path = self._create_tempfile()
        gv.save(out_path)
        self._assert_valid_png(out_path)

    # ----------------------------------------------------------------------
    # 4. BarGraph and legacy charts / graphs
    # ----------------------------------------------------------------------
    def test_line_chart(self):
        lc = pg.LineChart()
        lc.set_data([0.0, 1.0, 2.0], [10.0, 20.0, 15.0])
        out_path = self._create_tempfile()
        lc.savefig(out_path)
        self._assert_valid_png(out_path)



if __name__ == "__main__":
    unittest.main()
