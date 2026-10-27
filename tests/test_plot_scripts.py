#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import json
import os
import sys
import tempfile
import unittest
import unittest.mock as mock
import warnings

os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt

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


class TestPlotScripts(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        plt.close("all")
        self.tmpdir.cleanup()

    def test_plot_mat(self):
        """pdf-plot-mat.py generates a PNG file from a Matrix."""
        mat_path = os.path.join(self.tmpdir.name, "test.mat")
        out_png = os.path.join(self.tmpdir.name, "mat.png")

        mat = Matrix(3, 3)
        mat.set(0, 0, 1.0)
        mat.set(1, 1, 2.0)
        mat.set(2, 2, 3.0)
        mat.save(mat_path)

        mod = _load_script("pdf-plot-mat.py")
        sys.argv = ["pdf-plot-mat.py", mat_path, "-o", out_png]
        mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_vector(self):
        """pdf-plot-vector.py generates a PNG file from a Vector."""
        vct_path = os.path.join(self.tmpdir.name, "test.vct")
        out_png = os.path.join(self.tmpdir.name, "vct.png")

        vct = Vector([1.0, 2.0, 3.0, 4.0])
        vct.save(vct_path)

        mod = _load_script("pdf-plot-vector.py")
        sys.argv = ["pdf-plot-vector.py", vct_path, "-o", out_png]
        mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_dos(self):
        """pdf-plot-dos.py generates a PNG file from synthetic HDF5 results."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults.h5")
        out_png = os.path.join(self.tmpdir.name, "dos.png")
        csv_path = os.path.join(self.tmpdir.name, "dos.csv")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 1
        p.num_of_atoms = 1
        p.num_of_AOs = 5
        p.num_of_MOs = 5
        p.save_model()
        occ = Vector([2.0, 2.0, 1.0, 0.0, 0.0])
        p.save_occ("rks", occ)
        ene = Vector([-10.5, -1.2, -0.5, 0.2, 0.8])
        p.save_energy_level("rks", 1, ene)
        p.save_TEs(Vector([-100.0]))
        p._close()

        mod = _load_script("pdf-plot-dos.py")
        sys.argv = ["pdf-plot-dos.py", "-d", h5_path, "-o", out_png, "--csv", csv_path]
        mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_gvalues(self):
        """pdf-plot-gvalues.py generates a PNG file with mocked external component command."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults.h5")
        out_png = os.path.join(self.tmpdir.name, "g-values.png")
        csv_path = os.path.join(self.tmpdir.name, "g-values.csv")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 1
        p.num_of_atoms = 1
        p.num_of_AOs = 5
        p.num_of_MOs = 5
        p.save_model()
        occ = Vector([2.0, 2.0, 1.0, 0.0, 0.0])
        p.save_occ("rks", occ)
        ene = Vector([-10.5, -1.2, -0.5, 0.2, 0.8])
        p.save_energy_level("rks", 1, ene)
        p.save_TEs(Vector([-100.0]))
        p._close()

        called_cmds = []

        def fake_run_pdf(args):
            called_cmds.append(list(args))
            vtr = Vector([0.1, 0.2, 0.3, 0.4, 0.5])
            vtr.save("g-values.vtr")

        mod = _load_script("pdf-plot-gvalues.py")
        sys.argv = ["pdf-plot-gvalues.py", "-d", h5_path, "-o", out_png, "--csv", csv_path]
        with mock.patch("proteindf_tools.run_pdf", side_effect=fake_run_pdf):
            mod.main()

        self.assertEqual(len(called_cmds), 1)
        cmd = called_cmds[0]
        self.assertIn("component", cmd)
        self.assertIn("-i", cmd)
        self.assertEqual(cmd[cmd.index("-i") + 1], 1)
        self.assertIn("-S", cmd)
        self.assertEqual(cmd[cmd.index("-S") + 1], "CSC.mat")
        self.assertIn("-g", cmd)
        self.assertEqual(cmd[cmd.index("-g") + 1], "g-values.vtr")

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_pop(self):
        """pdf-plot-pop.py generates a PNG file from population data in an HDF5 file."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults.h5")
        out_png = os.path.join(self.tmpdir.name, "pop.png")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 1
        p.num_of_atoms = 2
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(Vector([-100.0]))

        # Save population data into HDF5 at pop/mulliken/<run_type>_<iteration>
        pop = Vector([0.45, -0.45])
        h5path = p._get_h5path_pop_mulliken_atom("rks", 1)
        p._write_vector(p._h5root, h5path, pop)
        p._close()

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "-o", out_png]
        mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_pop_missing_data(self):
        """pdf-plot-pop.py raises ValueError when population data is missing in HDF5."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults_nopop.h5")
        out_png = os.path.join(self.tmpdir.name, "pop_nopop.png")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 1
        p.num_of_atoms = 2
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(Vector([-100.0]))
        p._close()

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "-o", out_png]
        with self.assertRaises(ValueError) as ctx:
            mod.main()
        self.assertIn("Population data (Mulliken)", str(ctx.exception))
        self.assertIn("not found", str(ctx.exception))

    def test_plot_pop_cached_no_compute(self):
        """(a) When population is already in HDF5, do not run pop-mulliken."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults_cached.h5")
        out_png = os.path.join(self.tmpdir.name, "pop.png")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 2
        p.num_of_atoms = 2
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(Vector([-100.0, -100.5]))
        p.save_pop_mulliken_atom("rks", 1, Vector([0.4, -0.4]))
        p.save_pop_mulliken_atom("rks", 2, Vector([0.5, -0.5]))
        p._close()

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "-o", out_png]

        with mock.patch.object(mod.subprocess, "run") as mock_run:
            mod.main()
            mock_run.assert_not_called()

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    def test_plot_pop_compute_and_save(self):
        """(b) When population is missing in HDF5, compute with pop-mulliken and save to HDF5."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults_compute.h5")
        out_png = os.path.join(self.tmpdir.name, "pop.png")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 2
        p.num_of_atoms = 2
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(Vector([-100.0, -100.5]))

        S = SymmetricMatrix(4)
        for i in range(4):
            S.set(i, i, 1.0)
        p.save_s_matrix(S)

        P_in = SymmetricMatrix(4)
        P_in.set(0, 0, 1.0)
        p.save_input_density_matrix("rks", 1, P_in)

        P_last = SymmetricMatrix(4)
        P_last.set(1, 1, 1.0)
        p.save_density_matrix("rks", 2, P_last)
        p._close()

        dummy_output_itr1 = (
            " *** Mulliken population analysis (Gross atom population) ***\n"
            " index   atom         charge      net charge\n"
            "      1   H       0.600000     0.400000\n"
            "      2   H       1.400000    -0.400000\n"
        )
        dummy_output_itr2 = (
            " *** Mulliken population analysis (Gross atom population) ***\n"
            " index   atom         charge      net charge\n"
            "      1   H       0.500000     0.500000\n"
            "      2   H       1.500000    -0.500000\n"
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

            if fake_subprocess_run.call_count == 0:
                fake_subprocess_run.call_count += 1
                return mock.Mock(returncode=0, stdout=dummy_output_itr2, stderr="")
            else:
                fake_subprocess_run.call_count += 1
                return mock.Mock(returncode=0, stdout=dummy_output_itr1, stderr="")

        fake_subprocess_run.call_count = 0

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "-o", out_png]

        with mock.patch.object(mod.subprocess, "run", side_effect=fake_subprocess_run) as mock_run:
            mod.main()
            self.assertEqual(mock_run.call_count, 2)

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

        # Verify saved population in HDF5
        p_check = PdfParam_H5()
        p_check.open(h5_path)
        pop1 = p_check.get_pop_mulliken_atom("rks", 1)
        pop2 = p_check.get_pop_mulliken_atom("rks", 2)
        p_check._close()

        self.assertIsNotNone(pop1)
        self.assertAlmostEqual(pop1[0], 0.4)
        self.assertAlmostEqual(pop1[1], -0.4)

        self.assertIsNotNone(pop2)
        self.assertAlmostEqual(pop2[0], 0.5)
        self.assertAlmostEqual(pop2[1], -0.5)

    def test_plot_pop_readonly_fallback(self):
        """(c) When HDF5 is read-only, emit warning but plot successfully using in-memory data."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults_readonly.h5")
        out_png = os.path.join(self.tmpdir.name, "pop.png")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 2
        p.num_of_atoms = 2
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(Vector([-100.0, -100.5]))

        S = SymmetricMatrix(4)
        for i in range(4):
            S.set(i, i, 1.0)
        p.save_s_matrix(S)

        P_in = SymmetricMatrix(4)
        P_in.set(0, 0, 1.0)
        p.save_input_density_matrix("rks", 1, P_in)

        P_last = SymmetricMatrix(4)
        P_last.set(1, 1, 1.0)
        p.save_density_matrix("rks", 2, P_last)
        p._close()

        # Set read-only permissions
        os.chmod(h5_path, 0o444)

        dummy_output_itr1 = (
            " *** Mulliken population analysis (Gross atom population) ***\n"
            " index   atom         charge      net charge\n"
            "      1   H       0.600000     0.400000\n"
            "      2   H       1.400000    -0.400000\n"
        )
        dummy_output_itr2 = (
            " *** Mulliken population analysis (Gross atom population) ***\n"
            " index   atom         charge      net charge\n"
            "      1   H       0.500000     0.500000\n"
            "      2   H       1.500000    -0.500000\n"
        )

        def fake_subprocess_run(cmd, capture_output=True, text=True):
            if fake_subprocess_run.call_count == 0:
                fake_subprocess_run.call_count += 1
                return mock.Mock(returncode=0, stdout=dummy_output_itr2, stderr="")
            else:
                fake_subprocess_run.call_count += 1
                return mock.Mock(returncode=0, stdout=dummy_output_itr1, stderr="")

        fake_subprocess_run.call_count = 0

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "-o", out_png]

        try:
            with mock.patch.object(mod.subprocess, "run", side_effect=fake_subprocess_run):
                with warnings.catch_warnings(record=True) as warn_list:
                    warnings.simplefilter("always")
                    mod.main()

            # Ensure warning was emitted
            self.assertTrue(any("read-only" in str(w.message) for w in warn_list))

            # PNG should be generated successfully
            self.assertTrue(os.path.exists(out_png))
            self.assertGreater(os.path.getsize(out_png), 0)

            # Population was NOT saved to HDF5
            p_check = PdfParam_H5()
            p_check._state["h5filepath"] = h5_path
            import h5py
            p_check._h5file = h5py.File(h5_path, "r")
            p_check._h5root = p_check._h5file["/"]
            self.assertIsNone(p_check.get_pop_mulliken_atom("rks", 1))
            self.assertIsNone(p_check.get_pop_mulliken_atom("rks", 2))
            p_check._close()

        finally:
            os.chmod(h5_path, 0o644)

    def test_make_pop_diff_plot_data(self):
        """make_pop_diff_plot_data builds pdf-plot/1 schema dict with 1-based indices and diff values."""
        mod = _load_script("pdf-plot-pop.py")
        h5_path = "results.h5"
        elements = ["C", "H", "O"]
        pop_last = [1.2, -0.4, 0.8]
        pop_ref = [1.0, -0.5, 0.5]

        data = mod.make_pop_diff_plot_data(
            h5_path=h5_path,
            last_iteration=3,
            elements=elements,
            pop_vtr_last=pop_last,
            pop_vtr_ref=pop_ref,
            run_type="rks",
            reference_iteration=1,
        )

        self.assertEqual(data["schema"], "pdf-plot/1")
        self.assertEqual(data["kind"], "scatter")
        self.assertEqual(data["title"], "Mulliken charge difference")

        source = data["source"]
        self.assertEqual(source["tool"], "pdf-plot-pop.py")
        self.assertEqual(source["result"], "results.h5")
        self.assertEqual(source["run_type"], "rks")
        self.assertEqual(source["iteration"], 3)
        self.assertEqual(source["reference_iteration"], 1)

        axes = data["axes"]
        self.assertEqual(axes["x"], {"label": "atom index", "unit": None})
        self.assertEqual(
            axes["y"],
            {"label": "Mulliken charge difference (iteration 3 \u2212 iteration 1)", "unit": "e"},
        )

        atoms = data["atoms"]
        self.assertEqual(len(atoms), 3)
        for i, atom in enumerate(atoms):
            self.assertEqual(atom["index"], i + 1)
            self.assertEqual(atom["element"], elements[i])
            self.assertEqual(atom["label"], "")
            self.assertIsNone(atom["chain"])
            self.assertIsNone(atom["residue_number"])
            self.assertIsNone(atom["residue_name"])
            self.assertIsNone(atom["atom_name"])

        series = data["series"]
        self.assertEqual(len(series), 1)
        self.assertEqual(series[0]["name"], "Mulliken (diff)")
        points = series[0]["points"]
        self.assertEqual(len(points), 3)
        for i, pt in enumerate(points):
            self.assertEqual(pt["x"], i + 1)
            self.assertEqual(pt["atom"], i + 1)
            self.assertAlmostEqual(pt["y"], pop_last[i] - pop_ref[i], places=7)

        # Mismatch in lengths raises ValueError
        with self.assertRaises(ValueError):
            mod.make_pop_diff_plot_data(h5_path, 3, elements, pop_last, [1.0])
        with self.assertRaises(ValueError):
            mod.make_pop_diff_plot_data(h5_path, 3, ["C"], pop_last, pop_ref)

    def test_get_atom_elements_prioritizes_xyzc(self):
        """get_atom_elements prioritizes molecule/xyzc order over molecule.get_atom_list()."""
        mod = _load_script("pdf-plot-pop.py")

        class DummyMolecule:
            def get_atom_list(self):
                # Different order: C, H, O
                return [
                    bridge.Atom(symbol="C", position=bridge.Position(0.0, 0.0, 0.0), charge=6.0),
                    bridge.Atom(symbol="H", position=bridge.Position(0.0, 0.0, 1.0), charge=1.0),
                    bridge.Atom(symbol="O", position=bridge.Position(0.0, 0.0, 2.0), charge=8.0),
                ]

        class DummyParam:
            def __init__(self):
                self._h5root = {
                    "molecule/xyzc": [
                        [b"O", b"O1", 0.0, 0.0, 2.0, 8.0],
                        [b"C", b"C1", 0.0, 0.0, 0.0, 6.0],
                        [b"H", b"H1", 0.0, 0.0, 1.0, 1.0],
                    ]
                }
                self.molecule = DummyMolecule()

        param = DummyParam()
        elements = mod.get_atom_elements(param)
        self.assertEqual(elements, ["O", "C", "H"])

        # Also verify in JSON dict generated via make_pop_diff_plot_data
        data = mod.make_pop_diff_plot_data(
            h5_path="dummy.h5",
            last_iteration=2,
            elements=elements,
            pop_vtr_last=[0.5, 0.2, -0.7],
            pop_vtr_ref=[0.1, 0.1, -0.2],
        )
        self.assertEqual([a["element"] for a in data["atoms"]], ["O", "C", "H"])

    def test_get_atom_elements_fallback_warning(self):
        """get_atom_elements falls back to get_atom_list() with warning when molecule/xyzc is missing."""
        mod = _load_script("pdf-plot-pop.py")

        class DummyMolecule:
            def get_atom_list(self):
                return [bridge.Atom(symbol="N", position=bridge.Position(0.0, 0.0, 0.0), charge=7.0)]

        class DummyParam:
            def __init__(self):
                self._h5root = {}
                self.molecule = DummyMolecule()

        param = DummyParam()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            elements = mod.get_atom_elements(param)
            self.assertEqual(elements, ["N"])
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertIn("molecule/xyzc not found", str(w[-1].message))

    def test_plot_pop_json_cli(self):
        """pdf-plot-pop.py --json generates json directly and PNG without pop.mulliken.txt."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults.h5")
        out_png = os.path.join(self.tmpdir.name, "pop.png")
        out_json = os.path.join(self.tmpdir.name, "pop.json")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 2
        p.num_of_atoms = 3
        p.num_of_AOs = 6
        p.num_of_MOs = 6

        mol = bridge.AtomGroup()
        mol.set_atom(
            "0", bridge.Atom(symbol="C", position=bridge.Position(0.0, 0.0, 0.0), charge=6.0, name="C1")
        )
        mol.set_atom(
            "1", bridge.Atom(symbol="H", position=bridge.Position(0.0, 0.0, 1.0), charge=1.0, name="H1")
        )
        mol.set_atom(
            "2", bridge.Atom(symbol="X", position=bridge.Position(0.0, 0.0, 2.0), charge=0.0, name="X1")
        )
        p.molecule = mol
        p.save_model()
        p.save_TEs(Vector([-100.0, -100.5]))

        pop1 = Vector([0.1, -0.1, 0.0])
        pop2 = Vector([0.3, -0.25, 0.05])
        p._write_vector(p._h5root, p._get_h5path_pop_mulliken_atom("rks", 1), pop1)
        p._write_vector(p._h5root, p._get_h5path_pop_mulliken_atom("rks", 2), pop2)
        p._close()

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "--json", out_json, "-o", out_png]
        mod.main()

        # JSON and PNG should exist
        self.assertTrue(os.path.exists(out_json))
        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_json), 0)
        self.assertGreater(os.path.getsize(out_png), 0)

        # pop.mulliken.txt must NOT exist when --json is specified
        self.assertFalse(os.path.exists(os.path.join(self.tmpdir.name, "pop.mulliken.txt")))

        # Check JSON contents
        with open(out_json, "r", encoding="utf-8") as f:
            data = json.load(f)

        self.assertEqual(data["schema"], "pdf-plot/1")
        self.assertEqual(data["kind"], "scatter")
        self.assertEqual(data["source"]["iteration"], 2)
        self.assertEqual(data["source"]["reference_iteration"], 1)
        self.assertEqual(data["axes"]["y"]["label"], "Mulliken charge difference (iteration 2 \u2212 iteration 1)")

        expected_elements = ["C", "H", "X"]
        for i, atom in enumerate(data["atoms"]):
            self.assertEqual(atom["index"], i + 1)
            self.assertEqual(atom["element"], expected_elements[i])

        points = data["series"][0]["points"]
        for i, pt in enumerate(points):
            expected_diff = pop2[i] - pop1[i]
            self.assertAlmostEqual(pt["y"], expected_diff, places=7)

    def test_plot_pop_default_creates_txt(self):
        """pdf-plot-pop.py without --json still creates pop.mulliken.txt and pop.png."""
        h5_path = os.path.join(self.tmpdir.name, "pdfresults.h5")
        out_png = os.path.join(self.tmpdir.name, "pop.png")

        p = PdfParam_H5()
        p.open(h5_path)
        p.method = "rks"
        p.iterations = 1
        p.num_of_atoms = 2
        p.num_of_AOs = 4
        p.num_of_MOs = 4
        p.save_model()
        p.save_TEs(Vector([-100.0]))

        pop = Vector([0.45, -0.45])
        h5path = p._get_h5path_pop_mulliken_atom("rks", 1)
        p._write_vector(p._h5root, h5path, pop)
        p._close()

        mod = _load_script("pdf-plot-pop.py")
        sys.argv = ["pdf-plot-pop.py", "--h5", h5_path, "-o", out_png]
        mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertTrue(os.path.exists(os.path.join(self.tmpdir.name, "pop.mulliken.txt")))

    def test_plot_decaymat(self):
        """pdf-plot-decaymat.py generates a PNG file from HDF5 archive and Matrix."""
        db_path = os.path.join(self.tmpdir.name, "test.db")
        mat_path = os.path.join(self.tmpdir.name, "test.mat")
        out_png = os.path.join(self.tmpdir.name, "mat.png")

        param = PdfParam()
        param.comment = "H2 test system"
        param.num_of_atoms = 2
        param.num_of_AOs = 2
        param.num_of_MOs = 2
        param.method = "rks"
        param.iterations = 2
        param._set_scf_converged(True)

        mol = bridge.AtomGroup()
        a1 = bridge.Atom(
            symbol="H", position=bridge.Position([0.0, 0.0, 0.0]), charge=1.0, name="H1"
        )
        a2 = bridge.Atom(
            symbol="H", position=bridge.Position([0.0, 0.0, 0.74]), charge=1.0, name="H2"
        )
        mol.set_atom("0", a1)
        mol.set_atom("1", a2)
        param.molecule = mol

        bs = BasisSet()
        bs.name = "sto-3g-H"
        cgto = ContractedGTO("s", 1)
        cgto[0] = PrimitiveGTO(3.42525, 0.15433)
        bs.append(cgto)
        param.set_basisset("H", bs)
        param.set_basisset("H@H1", bs)
        param.set_basisset("H@H2", bs)
        param.set_total_energy(1, -1.10)
        param.set_total_energy(2, -1.17)

        h5_path = os.path.join(self.tmpdir.name, "test.h5")
        h = PdfParam_H5(param)
        h.open(h5_path)
        h.save_model()
        h.save_TEs(Vector([-1.10, -1.17]))
        h._close()

        mat = Matrix(2, 2)
        mat.set(0, 0, 1.0)
        mat.set(0, 1, 0.5)
        mat.set(1, 0, 0.5)
        mat.set(1, 1, 1.0)
        mat.save(mat_path)

        mod = _load_script("pdf-plot-decaymat.py")

        # .db raises ValueError
        sys.argv = ["pdf-plot-decaymat.py", "-d", db_path, mat_path, "-o", out_png]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        # .h5 succeeds
        sys.argv = ["pdf-plot-decaymat.py", "-d", h5_path, mat_path, "-o", out_png]
        mod.main()

        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)

    @unittest.skipUnless(
        os.environ.get("PDF_HOME")
        and os.path.exists(os.path.join(os.environ.get("PDF_HOME", ""), "data/basis2")),
        "PDF_HOME data/basis2 is not available",
    )
    def test_plot_basisset(self):
        """pdf-plot-basisset.py generates a PNG file when PDF_HOME data/basis2 is available."""
        mod = _load_script("pdf-plot-basisset.py")
        sys.argv = ["pdf-plot-basisset.py", "O-STO-3G.H"]
        mod.main()

        out_png = os.path.join(self.tmpdir.name, "basis_O-STO-3G.H.png")
        self.assertTrue(os.path.exists(out_png))
        self.assertGreater(os.path.getsize(out_png), 0)


if __name__ == "__main__":
    unittest.main()
