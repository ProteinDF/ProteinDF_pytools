# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import importlib.util
import os
import sys
import tempfile
import unittest
import yaml

from proteindf_tools.vector import Vector
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def is_valid_png(filepath):
    if not os.path.isfile(filepath):
        return False
    with open(filepath, "rb") as f:
        header = f.read(8)
    return header == b"\x89PNG\r\n\x1a\n"


class TestPlotEnergy(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        self.mod = _load_script("pdf-plot-energy.py")

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def _create_h5(self, h5_path, tes):
        p = PdfParam_H5()
        p.open(h5_path)
        p.num_of_atoms = 1
        p.num_of_AOs = 1
        p.num_of_MOs = 1
        p.method = "rks"
        p.iterations = len(tes)
        p.save_model()
        if len(tes) > 0:
            p.save_TEs(Vector(tes))
        p._close()

    def test_plot_energy_default(self):
        """Default run without options creates energy.png from pdfresults.h5."""
        h5_path = "pdfresults.h5"
        self._create_h5(h5_path, [-100.0, -100.5, -100.8])

        sys.argv = ["pdf-plot-energy.py"]
        self.mod.main()

        out_png = "energy.png"
        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(is_valid_png(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_energy_custom_output(self):
        """-o option specifies the output PNG file path."""
        h5_path = os.path.join(self.tmpdir.name, "custom.h5")
        self._create_h5(h5_path, [-50.0, -51.0])

        out_png = os.path.join(self.tmpdir.name, "my_energy.png")
        sys.argv = ["pdf-plot-energy.py", h5_path, "-o", out_png]
        self.mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(is_valid_png(out_png))
        self.assertFalse(os.path.exists(os.path.join(self.tmpdir.name, "energy.png")))

    def test_plot_energy_yaml_only(self):
        """--yaml option writes YAML and suppresses default PNG creation."""
        h5_path = os.path.join(self.tmpdir.name, "results.h5")
        energies = [-100.1234, -101.5678, -101.9999]
        self._create_h5(h5_path, energies)

        out_yaml = os.path.join(self.tmpdir.name, "energy.yaml")
        sys.argv = ["pdf-plot-energy.py", h5_path, "--yaml", out_yaml]
        self.mod.main()

        self.assertTrue(os.path.exists(out_yaml))
        self.assertFalse(os.path.exists("energy.png"))

        with open(out_yaml, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)

        self.assertEqual(data["schema"], "pdf-plot/1")
        self.assertEqual(data["kind"], "total_energy")
        self.assertEqual(data["title"], "total energy history")
        self.assertEqual(data["source"], h5_path)
        self.assertEqual(
            data["axes"],
            {
                "x": {"label": "iteration", "unit": None},
                "y": {"label": "energy / a.u.", "unit": "a.u."},
            },
        )
        self.assertEqual(len(data["series"]), 1)
        series0 = data["series"][0]
        self.assertEqual(series0["name"], "total_energy")
        points = series0["points"]
        self.assertEqual(len(points), len(energies))
        for i, (pt, expected_e) in enumerate(zip(points, energies)):
            self.assertEqual(pt["iteration"], i + 1)
            self.assertAlmostEqual(pt["energy"], expected_e, places=7)

    def test_plot_energy_both_png_and_yaml(self):
        """When both -o and --yaml are specified, both files are generated."""
        h5_path = os.path.join(self.tmpdir.name, "results.h5")
        self._create_h5(h5_path, [-10.0, -10.5])

        out_png = os.path.join(self.tmpdir.name, "both.png")
        out_yaml = os.path.join(self.tmpdir.name, "both.yaml")
        sys.argv = [
            "pdf-plot-energy.py",
            h5_path,
            "-o",
            out_png,
            "--yaml",
            out_yaml,
        ]
        self.mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(is_valid_png(out_png))
        self.assertTrue(os.path.exists(out_yaml))

    def test_plot_energy_skip_unreadable_iteration(self):
        """Unreadable/None iterations are skipped with a warning log."""
        h5_path = os.path.join(self.tmpdir.name, "partial.h5")
        self._create_h5(h5_path, [-100.0, -101.0, -102.0])

        p = PdfParam_H5()
        p.open(h5_path)
        # Simulate iteration 2 returning None
        p._data["TEs"][2] = None

        points = self.mod.read_total_energies(p)
        self.assertEqual(len(points), 2)
        self.assertEqual(points[0], (1, -100.0))
        self.assertEqual(points[1], (3, -102.0))
