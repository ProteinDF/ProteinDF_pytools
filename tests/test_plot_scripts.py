#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import os
import sys
import tempfile
import unittest
import unittest.mock as mock

os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.matrix import Matrix
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
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

    # scripts/pdf-plot-decaymat.py:129:
    # d = distance_mat.get(atom_i, atom_j)
    # raises TypeError: '<' not supported between instances of 'NoneType' and 'NoneType'
    # (OrbInfo.__init__ does not handle PdfArchive / PdfArchive_Sqlite3, leaving orbital info empty
    # so orb_info.get_atom_id(i) returns None)
    @unittest.expectedFailure
    def test_plot_decaymat(self):
        """pdf-plot-decaymat.py generates a PNG file from PdfArchive and Matrix."""
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
        param.set_basisset("H@0", bs)
        param.set_basisset("H@1", bs)
        param.set_total_energy(1, -1.10)
        param.set_total_energy(2, -1.17)

        archive = PdfArchive_Sqlite3(db_path, pdfparam=param)

        mat = Matrix(2, 2)
        mat.set(0, 0, 1.0)
        mat.set(0, 1, 0.5)
        mat.set(1, 0, 0.5)
        mat.set(1, 1, 1.0)
        mat.save(mat_path)

        mod = _load_script("pdf-plot-decaymat.py")
        sys.argv = ["pdf-plot-decaymat.py", "-d", db_path, mat_path, "-o", out_png]
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
