#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import tempfile
import unittest
import matplotlib.pyplot as plt

from proteindf_tools.pdfgraph import (
    DfEnergyLevelHistoryGraphH,
    DfEnergyLevelHistoryGraphV,
    DfGraph,
    DfLineChart,
    DfPopulationGraph,
    DfTotalEnergyHistGraph,
)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


class BaseGraphTestCase(unittest.TestCase):
    def setUp(self):
        self._temp_dir = tempfile.TemporaryDirectory()
        self.temp_dir = self._temp_dir.name

    def tearDown(self):
        plt.close("all")
        self._temp_dir.cleanup()

    def assert_valid_png(self, path):
        self.assertTrue(os.path.isfile(path))
        self.assertGreater(os.path.getsize(path), 0)
        with open(path, "rb") as f:
            self.assertEqual(f.read(8), PNG_SIGNATURE)


class TestDfGraph(BaseGraphTestCase):
    def test_default_properties(self):
        g = DfGraph()
        self.assertEqual(g.title, "")
        self.assertEqual(g.xlabel, "")
        self.assertEqual(g.ylabel, "")
        self.assertFalse(g.is_draw_grid)
        self.assertIsNone(g.xmin)
        self.assertIsNone(g.xmax)
        self.assertIsNone(g.ymin)
        self.assertIsNone(g.ymax)
        self.assertEqual(len(g.size), 2)
        self.assertIsNone(g.aspect)
        self.assertEqual(g.xscale, "linear")
        self.assertEqual(g.yscale, "linear")
        self.assertIsNone(g.xticks)
        self.assertIsNone(g.yticks)
        self.assertIsNone(g.xticklabels)
        self.assertIsNone(g.yticklabels)
        self.assertFalse(g.legend)

    def test_set_properties(self):
        g = DfGraph()
        g.title = "Sample Title"
        self.assertEqual(g.title, "Sample Title")
        g.xlabel = "X-axis"
        self.assertEqual(g.xlabel, "X-axis")
        g.ylabel = "Y-axis"
        self.assertEqual(g.ylabel, "Y-axis")
        g.is_draw_grid = True
        self.assertTrue(g.is_draw_grid)
        g.xmin = 0.5
        self.assertEqual(g.xmin, 0.5)
        g.xmax = 10.5
        self.assertEqual(g.xmax, 10.5)
        g.ymin = -2.0
        self.assertEqual(g.ymin, -2.0)
        g.ymax = 8.0
        self.assertEqual(g.ymax, 8.0)
        g.size = (8.0, 6.0)
        self.assertEqual(tuple(g.size), (8.0, 6.0))
        g.aspect = "equal"
        self.assertEqual(g.aspect, "equal")
        g.xscale = "log"
        self.assertEqual(g.xscale, "log")
        g.yscale = "log"
        self.assertEqual(g.yscale, "log")
        g.xticks = [1, 5, 10]
        self.assertEqual(g.xticks, [1, 5, 10])
        g.yticks = [0, 2, 4]
        self.assertEqual(g.yticks, [0, 2, 4])
        g.xticklabels = ["one", "five", "ten"]
        self.assertEqual(g.xticklabels, ["one", "five", "ten"])
        g.yticklabels = ["zero", "two", "four"]
        self.assertEqual(g.yticklabels, ["zero", "two", "four"])
        g.legend = True
        self.assertTrue(g.legend)

    def test_invalid_scale_values(self):
        g = DfGraph()
        g.xscale = "log"
        self.assertEqual(g.xscale, "log")
        g.xscale = "invalid_scale"
        self.assertEqual(g.xscale, "linear")
        g.xscale = 123
        self.assertEqual(g.xscale, "linear")

        g.yscale = "log"
        self.assertEqual(g.yscale, "log")
        g.yscale = "invalid_scale"
        self.assertEqual(g.yscale, "linear")
        g.yscale = None
        self.assertEqual(g.yscale, "linear")

    def test_invalid_types_raise_assertion_error(self):
        g = DfGraph()
        with self.assertRaises(AssertionError):
            g.size = [8.0, 6.0]

        with self.assertRaises(AssertionError):
            g.xticks = (1, 2, 3)

        with self.assertRaises(AssertionError):
            g.yticks = "not_a_list"

    def test_save(self):
        g = DfGraph()
        g.title = "Save Test"
        g.xlabel = "X"
        g.ylabel = "Y"
        g.is_draw_grid = True
        g.xmin = 0.0
        g.xmax = 5.0
        g.ymin = 0.0
        g.ymax = 10.0
        g.xticks = [0, 2, 4]
        g.yticks = [0, 5, 10]
        g.xticklabels = ["A", "B", "C"]
        g.yticklabels = ["Low", "Mid", "High"]
        g.aspect = "auto"
        g.legend = True
        out_path = os.path.join(self.temp_dir, "dfgraph.png")
        g.save(out_path)
        self.assert_valid_png(out_path)


class TestDfTotalEnergyHistGraph(BaseGraphTestCase):
    def test_default_properties(self):
        g = DfTotalEnergyHistGraph()
        self.assertEqual(g.title, "total energy history")
        self.assertEqual(g.xlabel, "iteration")
        self.assertEqual(g.ylabel, "energy / a.u.")
        self.assertTrue(g.is_draw_grid)

    def test_load_data_and_save(self):
        data_path = os.path.join(self.temp_dir, "energy.dat")
        with open(data_path, "w") as f:
            f.write("1, -120.500\n")
            f.write("2, -120.750\n")
            f.write("3, -120.755\n")

        g = DfTotalEnergyHistGraph()
        g.load_data(data_path)
        self.assertEqual(g._x, [1.0, 2.0, 3.0])
        self.assertEqual(g._y, [-120.500, -120.750, -120.755])

        out_path = os.path.join(self.temp_dir, "energy.png")
        g.save(out_path)
        self.assert_valid_png(out_path)


class TestDfPopulationGraph(BaseGraphTestCase):
    def test_default_properties(self):
        g = DfPopulationGraph()
        self.assertEqual(g.title, "Population")
        self.assertEqual(g.xlabel, "Order of atom list")
        self.assertEqual(g.ylabel, "charge")

    def test_load_data_and_save(self):
        data_path = os.path.join(self.temp_dir, "population.dat")
        with open(data_path, "w") as f:
            f.write("0, 0.45\n")
            f.write("1, -0.45\n")
            f.write("2, 0.12\n")

        g = DfPopulationGraph()
        g.load_data(data_path)
        self.assertEqual(g._x, [0.0, 1.0, 2.0])
        self.assertEqual(g._y, [0.45, -0.45, 0.12])

        out_path = os.path.join(self.temp_dir, "population.png")
        g.save(out_path)
        self.assert_valid_png(out_path)


class TestDfLineChart(BaseGraphTestCase):
    def test_add_data_and_save(self):
        chart = DfLineChart()
        chart.add_data([1.0, 2.0, 3.0], [10.0, 20.0, 30.0])
        chart.add_data([1.0, 2.0, 3.0], [30.0, 25.0, 20.0])
        self.assertEqual(len(chart._data), 2)
        self.assertEqual(chart._data[0], ([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]))

        out_path = os.path.join(self.temp_dir, "line_chart.png")
        chart.save(out_path)
        self.assert_valid_png(out_path)


class TestDfEnergyLevelHistoryGraph(BaseGraphTestCase):
    def test_default_properties_h(self):
        gh = DfEnergyLevelHistoryGraphH()
        self.assertEqual(gh.title, "Energy Level")
        self.assertEqual(gh.xlabel, "iteration")
        self.assertEqual(gh.ylabel, "energy / eV")
        self.assertTrue(gh.is_draw_grid)
        self.assertIsNone(gh.xmin)
        self.assertIsNone(gh.xmax)
        self.assertEqual(gh.ymin, -20.0)
        self.assertEqual(gh.ymax, 5.0)
        self.assertEqual(gh._HOMO_level, -1)
        self.assertEqual(gh._LUMO_level, -1)
        self.assertEqual(gh._max_itr, 0)
        self.assertEqual(gh._select_iterations, [])

    def test_default_properties_v(self):
        gv = DfEnergyLevelHistoryGraphV()
        self.assertEqual(gv.title, "Energy Level")
        self.assertEqual(gv.xlabel, "energy / eV")
        self.assertEqual(gv.ylabel, "iteration")
        self.assertTrue(gv.is_draw_grid)
        self.assertEqual(gv.xmin, -20.0)
        self.assertEqual(gv.xmax, 5.0)
        self.assertIsNone(gv.ymin)
        self.assertIsNone(gv.ymax)
        self.assertEqual(gv._HOMO_level, -1)
        self.assertEqual(gv._LUMO_level, -1)

    def test_load_data_and_save_h(self):
        data_path = os.path.join(self.temp_dir, "eigvals_hist.dat")
        with open(data_path, "w") as f:
            f.write("1, 0, -15.0\n")
            f.write("1, 1, -10.0\n")
            f.write("1, 2, -2.0\n")
            f.write("\n")
            f.write("2, 0, -15.1\n")
            f.write("2, 1, -9.9\n")
            f.write("2, 2, -1.9\n")

        gh = DfEnergyLevelHistoryGraphH()
        gh.load_data(data_path)
        self.assertEqual(gh._max_itr, 2)
        self.assertEqual(len(gh._data), 6)
        self.assertEqual(gh._data[0], [1, 0.0, -15.0])

        gh.set_HOMO_level(1)
        gh.set_LUMO_level(2)
        self.assertEqual(gh._HOMO_level, 1)
        self.assertEqual(gh._LUMO_level, 2)

        # 全反復の描画
        out_path = os.path.join(self.temp_dir, "eigvals_h.png")
        gh.save(out_path)
        self.assert_valid_png(out_path)

        # 反復指定の描画
        gh.select_iterations([2])
        self.assertEqual(gh._select_iterations, {2})
        out_path2 = os.path.join(self.temp_dir, "eigvals_h_itr2.png")
        gh.save(out_path2)
        self.assert_valid_png(out_path2)

    def test_load_data_and_save_v(self):
        data_path = os.path.join(self.temp_dir, "eigvals_hist.dat")
        with open(data_path, "w") as f:
            f.write("1, 0, -15.0\n")
            f.write("1, 1, -10.0\n")
            f.write("1, 2, -2.0\n")

        gv = DfEnergyLevelHistoryGraphV()
        gv.load_data(data_path)
        gv.set_HOMO_level(0)
        gv.set_LUMO_level(1)
        gv.select_iterations([1])

        out_path = os.path.join(self.temp_dir, "eigvals_v.png")
        gv.save(out_path)
        self.assert_valid_png(out_path)


class TestDataFileErrorHandling(BaseGraphTestCase):
    def test_nonexistent_data_file(self):
        nonexistent = os.path.join(self.temp_dir, "nonexistent.dat")

        g_energy = DfTotalEnergyHistGraph()
        with self.assertRaises(IOError):
            g_energy.load_data(nonexistent)

        g_pop = DfPopulationGraph()
        with self.assertRaises(IOError):
            g_pop.load_data(nonexistent)

        g_energy_h = DfEnergyLevelHistoryGraphH()
        with self.assertRaises(IOError):
            g_energy_h.load_data(nonexistent)

    def test_empty_data_file(self):
        empty_path = os.path.join(self.temp_dir, "empty.dat")
        open(empty_path, "w").close()

        g_energy = DfTotalEnergyHistGraph()
        g_energy.load_data(empty_path)
        self.assertEqual(g_energy._x, [])
        self.assertEqual(g_energy._y, [])
        out_energy = os.path.join(self.temp_dir, "empty_energy.png")
        g_energy.save(out_energy)
        self.assert_valid_png(out_energy)

        g_pop = DfPopulationGraph()
        g_pop.load_data(empty_path)
        self.assertEqual(g_pop._x, [])
        self.assertEqual(g_pop._y, [])
        out_pop = os.path.join(self.temp_dir, "empty_pop.png")
        g_pop.save(out_pop)
        self.assert_valid_png(out_pop)

        g_energy_h = DfEnergyLevelHistoryGraphH()
        g_energy_h.load_data(empty_path)
        self.assertEqual(g_energy_h._data, [])
        out_h = os.path.join(self.temp_dir, "empty_h.png")
        g_energy_h.save(out_h)
        self.assert_valid_png(out_h)

    def test_insufficient_columns(self):
        bad_path = os.path.join(self.temp_dir, "bad.dat")
        with open(bad_path, "w") as f:
            f.write("1\n")

        g_energy = DfTotalEnergyHistGraph()
        with self.assertRaises(IndexError):
            g_energy.load_data(bad_path)

        g_pop = DfPopulationGraph()
        with self.assertRaises(IndexError):
            g_pop.load_data(bad_path)

        g_energy_h = DfEnergyLevelHistoryGraphH()
        with self.assertRaises(ValueError):
            g_energy_h.load_data(bad_path)


if __name__ == "__main__":
    unittest.main()
