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

import glob
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


class TestPlotElevelHistory(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        self.mod = _load_script("pdf-plot-elevel-history.py")

    def tearDown(self):
        sys.argv = self._orig_argv
        os.chdir(self._orig_cwd)
        self.tmpdir.cleanup()

    def _create_h5(self, h5_path, method="rks", iterations=3, levels_per_itr=None, occ_dict=None):
        p = PdfParam_H5()
        p.open(h5_path)
        p.num_of_atoms = 1
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.method = method
        p.iterations = iterations
        p.save_model()
        p.save_TEs(Vector([-1.0] * iterations))

        if levels_per_itr is None:
            # default: 4 levels per iteration
            levels_per_itr = {}
            rts = ["rks"] if method == "rks" else ["uks_alpha", "uks_beta"]
            for rt in rts:
                levels_per_itr[rt] = {}
                for itr in range(1, iterations + 1):
                    levels_per_itr[rt][itr] = [-10.0 + i + itr * 0.1 for i in range(4)]

        for rt, itr_map in levels_per_itr.items():
            for itr, vals in itr_map.items():
                if vals is not None:
                    p.save_energy_level(rt, itr, Vector(vals))

        if occ_dict is not None:
            for rt, occ in occ_dict.items():
                p.save_occ(rt, Vector(occ))
        else:
            rts = ["rks"] if method == "rks" else ["uks_alpha", "uks_beta"]
            for rt in rts:
                occ = [2.0, 2.0, 0.0, 0.0] if method == "rks" else [1.0, 1.0, 0.0, 0.0]
                p.save_occ(rt, Vector(occ))

        p._close()

    def test_plot_elevel_history_default(self):
        """Default run without options creates eigvals_hist.<run_type>.png from pdfresults.h5."""
        h5_path = "pdfresults.h5"
        self._create_h5(h5_path, method="rks", iterations=3)

        sys.argv = ["pdf-plot-elevel-history.py"]
        self.mod.main()

        out_png = "eigvals_hist.rks.png"
        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(is_valid_png(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

        # Ensure no intermediate .dat files are left
        dat_files = glob.glob("*.dat")
        self.assertEqual(dat_files, [])

    def test_plot_elevel_history_custom_output_single(self):
        """-o option specifies the output PNG file path for single run_type."""
        h5_path = os.path.join(self.tmpdir.name, "custom.h5")
        self._create_h5(h5_path, method="rks", iterations=2)

        out_png = os.path.join(self.tmpdir.name, "custom.png")
        sys.argv = ["pdf-plot-elevel-history.py", h5_path, "-o", out_png]
        self.mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(is_valid_png(out_png))
        self.assertFalse(os.path.exists("eigvals_hist.rks.png"))
        self.assertEqual(glob.glob("*.dat"), [])

    def test_plot_elevel_history_multiple_runtypes_custom_output(self):
        """When multiple run_types exist, -o uses <stem>.<run_type><ext>."""
        h5_path = os.path.join(self.tmpdir.name, "multi.h5")
        self._create_h5(h5_path, method="uks", iterations=2)

        out_png = os.path.join(self.tmpdir.name, "multi.png")
        sys.argv = ["pdf-plot-elevel-history.py", h5_path, "-o", out_png]
        self.mod.main()

        alpha_png = os.path.join(self.tmpdir.name, "multi.uks_alpha.png")
        beta_png = os.path.join(self.tmpdir.name, "multi.uks_beta.png")
        self.assertTrue(os.path.exists(alpha_png))
        self.assertTrue(is_valid_png(alpha_png))
        self.assertTrue(os.path.exists(beta_png))
        self.assertTrue(is_valid_png(beta_png))
        self.assertFalse(os.path.exists(out_png))
        self.assertEqual(glob.glob("*.dat"), [])

    def test_plot_elevel_history_multiple_runtypes_default_output(self):
        """When multiple run_types exist and -o is omitted, default names are used."""
        h5_path = "pdfresults.h5"
        self._create_h5(h5_path, method="uks", iterations=2)

        sys.argv = ["pdf-plot-elevel-history.py"]
        self.mod.main()

        alpha_png = "eigvals_hist.uks_alpha.png"
        beta_png = "eigvals_hist.uks_beta.png"
        self.assertTrue(os.path.exists(alpha_png))
        self.assertTrue(is_valid_png(alpha_png))
        self.assertTrue(os.path.exists(beta_png))
        self.assertTrue(is_valid_png(beta_png))
        self.assertEqual(glob.glob("*.dat"), [])

    def test_plot_elevel_history_yaml_only_single(self):
        """--yaml writes YAML and suppresses PNG creation for single run_type."""
        h5_path = os.path.join(self.tmpdir.name, "results.h5")
        levels_map = {
            "rks": {
                1: [-10.5, -5.2, 0.1, 1.2],
                2: [-10.4, -5.1, 0.2, 1.3],
            }
        }
        self._create_h5(
            h5_path,
            method="rks",
            iterations=2,
            levels_per_itr=levels_map,
            occ_dict={"rks": [2.0, 2.0, 0.0, 0.0]},
        )

        out_yaml = os.path.join(self.tmpdir.name, "out.yaml")
        sys.argv = ["pdf-plot-elevel-history.py", h5_path, "--yaml", out_yaml]
        self.mod.main()

        self.assertTrue(os.path.exists(out_yaml))
        self.assertFalse(os.path.exists("eigvals_hist.rks.png"))
        self.assertEqual(glob.glob("*.png"), [])
        self.assertEqual(glob.glob("*.dat"), [])

        with open(out_yaml, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)

        self.assertEqual(data["schema"], "pdf-plot/1")
        self.assertEqual(data["kind"], "elevel_history")
        self.assertEqual(data["title"], "energy level history")
        self.assertEqual(data["source"], h5_path)
        self.assertEqual(data["run_type"], "rks")
        self.assertEqual(data["HOMO_level"], 1)
        self.assertEqual(
            data["axes"],
            {
                "x": {"label": "iteration", "unit": None},
                "y": {"label": "energy / eV", "unit": "eV"},
            },
        )
        self.assertEqual(len(data["series"]), 1)
        series0 = data["series"][0]
        self.assertEqual(series0["name"], "eigvals_hist.rks")
        points = series0["points"]
        self.assertEqual(len(points), 8)

        # Check points iteration, level, and eV values
        expected_points = []
        for itr in (1, 2):
            for lvl, e in enumerate(levels_map["rks"][itr]):
                expected_points.append({"iteration": itr, "level": lvl, "energy_ev": e * 27.2116})

        for actual, expected in zip(points, expected_points):
            self.assertEqual(actual["iteration"], expected["iteration"])
            self.assertEqual(actual["level"], expected["level"])
            self.assertAlmostEqual(actual["energy_ev"], expected["energy_ev"], places=6)

    def test_plot_elevel_history_yaml_only_multiple(self):
        """--yaml with multiple run_types generates <stem>.<run_type><ext>."""
        h5_path = os.path.join(self.tmpdir.name, "results_multi.h5")
        self._create_h5(h5_path, method="uks", iterations=2)

        out_yaml = os.path.join(self.tmpdir.name, "out.yaml")
        sys.argv = ["pdf-plot-elevel-history.py", h5_path, "--yaml", out_yaml]
        self.mod.main()

        alpha_yaml = os.path.join(self.tmpdir.name, "out.uks_alpha.yaml")
        beta_yaml = os.path.join(self.tmpdir.name, "out.uks_beta.yaml")
        self.assertTrue(os.path.exists(alpha_yaml))
        self.assertTrue(os.path.exists(beta_yaml))
        self.assertFalse(os.path.exists(out_yaml))
        self.assertEqual(glob.glob("*.png"), [])

        with open(alpha_yaml, "r", encoding="utf-8") as f:
            data_alpha = yaml.safe_load(f)
        self.assertEqual(data_alpha["schema"], "pdf-plot/1")
        self.assertEqual(data_alpha["kind"], "elevel_history")
        self.assertEqual(data_alpha["run_type"], "uks_alpha")
        self.assertEqual(data_alpha["HOMO_level"], 1)

    def test_plot_elevel_history_both_png_and_yaml(self):
        """When both -o and --yaml are specified, both files are generated."""
        h5_path = os.path.join(self.tmpdir.name, "both.h5")
        self._create_h5(h5_path, method="rks", iterations=2)

        out_png = os.path.join(self.tmpdir.name, "both.png")
        out_yaml = os.path.join(self.tmpdir.name, "both.yaml")
        sys.argv = [
            "pdf-plot-elevel-history.py",
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
        self.assertEqual(glob.glob("*.dat"), [])

    def test_plot_elevel_history_specify_run_type(self):
        """-r option selects a single run_type."""
        h5_path = os.path.join(self.tmpdir.name, "select.h5")
        self._create_h5(h5_path, method="uks", iterations=2)

        out_png = os.path.join(self.tmpdir.name, "beta.png")
        out_yaml = os.path.join(self.tmpdir.name, "beta.yaml")
        sys.argv = [
            "pdf-plot-elevel-history.py",
            h5_path,
            "-r",
            "uks_beta",
            "-o",
            out_png,
            "--yaml",
            out_yaml,
        ]
        self.mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(is_valid_png(out_png))
        self.assertTrue(os.path.exists(out_yaml))
        # Since only uks_beta was selected, filenames are not suffixed with .uks_beta
        self.assertFalse(os.path.exists(os.path.join(self.tmpdir.name, "beta.uks_beta.png")))
        self.assertFalse(os.path.exists(os.path.join(self.tmpdir.name, "beta.uks_alpha.png")))

    def test_plot_elevel_history_skip_unreadable_iteration(self):
        """Unreadable/None iterations are skipped with a warning log."""
        h5_path = os.path.join(self.tmpdir.name, "partial.h5")
        levels_map = {
            "rks": {
                1: [-10.0, -5.0, 0.0, 1.0],
                2: None,  # unreadable iteration
                3: [-10.2, -5.1, -0.1, 0.9],
            }
        }
        self._create_h5(h5_path, method="rks", iterations=3, levels_per_itr=levels_map)

        p = PdfParam_H5()
        p.open(h5_path)
        points = self.mod.read_energy_levels(p, "rks")
        p._close()

        # Should have 4 points from itr 1 and 4 points from itr 3 (total 8 points)
        self.assertEqual(len(points), 8)
        itrs = {pt[0] for pt in points}
        self.assertEqual(itrs, {1, 3})


if __name__ == "__main__":
    unittest.main()
