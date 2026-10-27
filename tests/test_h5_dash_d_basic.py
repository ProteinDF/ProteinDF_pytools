#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.matrix import Matrix
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5
from proteindf_tools.vector import Vector

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))

TARGET_SCRIPTS = [
    "pdf-info-geom.py",
    "pdf-info-xyz.py",
    "pdf-info-orb.py",
    "pdf-atom-index.py",
    "pdf-plot-decaymat.py",
    "pdf-pop-ridge.py",
    "pdf-pop-lasso.py",
    "pdf-pop-elasticnet.py",
    "pdf-pop-resp.py",
    "pdf-pop-rrms.py",
    "pdf-pop-classo.py",
]


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _create_synthetic_system(tmpdir):
    """Create identical synthetic water systems (O, H, H) in .db and .h5 formats."""
    mol = bridge.AtomGroup()
    a_o = bridge.Atom(symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="0")
    a_h1 = bridge.Atom(
        symbol="H", position=bridge.Position([0.0, 1.43, -1.1]), charge=1.0, name="1"
    )
    a_h2 = bridge.Atom(
        symbol="H", position=bridge.Position([0.0, -1.43, -1.1]), charge=1.0, name="2"
    )
    mol.set_atom("0", a_o)
    mol.set_atom("1", a_h1)
    mol.set_atom("2", a_h2)

    bs_o = BasisSet("sto-3g-O", 1)
    cg_o = ContractedGTO("s", 1)
    cg_o[0] = PrimitiveGTO(130.7, 0.15)
    bs_o[0] = cg_o

    bs_h = BasisSet("sto-3g-H", 1)
    cg_h = ContractedGTO("s", 1)
    cg_h[0] = PrimitiveGTO(3.4, 0.15)
    bs_h[0] = cg_h

    # Synthetic parameter & SQLite archive
    param = PdfParam()
    param.comment = "synthetic H2O system"
    param.num_of_atoms = 3
    param.num_of_AOs = 3
    param.num_of_MOs = 3
    param.method = "rks"
    param.iterations = 1
    param.guess = "none"
    param._set_scf_converged(True)
    param.molecule = mol
    param.set_basisset("O@0", bs_o)
    param.set_basisset("H@1", bs_h)
    param.set_basisset("H@2", bs_h)
    param.set_total_energy(1, -75.0)

    db_path = os.path.join(tmpdir, "water.db")
    db = PdfArchive_Sqlite3(db_path, pdfparam=param)
    db.set_occupations("rks", Vector([2.0, 2.0, 0.0]))
    db.set_energylevel("rks", 1, Vector([-19.0, -1.0, 0.5]))

    # Synthetic HDF5 archive
    h5_path = os.path.join(tmpdir, "water.h5")
    p5 = PdfParam_H5()
    p5.open(h5_path)
    p5.comment = "synthetic H2O system"
    p5.num_of_atoms = 3
    p5.num_of_AOs = 3
    p5.num_of_MOs = 3
    p5.method = "rks"
    p5.iterations = 1
    p5.guess = "none"
    p5._set_scf_converged(True)
    p5.molecule = mol
    p5.set_basisset("O@0", bs_o)
    p5.set_basisset("H@1", bs_h)
    p5.set_basisset("H@2", bs_h)
    p5.save_model()
    p5.save_TEs(Vector([-75.0]))
    p5.save_occ("rks", Vector([2.0, 2.0, 0.0]))
    p5.save_energy_level("rks", 1, Vector([-19.0, -1.0, 0.5]))
    p5._close()

    # Query molecule for pdf-atom-index
    query_mol = bridge.AtomGroup()
    grp = bridge.AtomGroup()
    target_atom = bridge.Atom(
        symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="0"
    )
    grp.set_atom("0", target_atom)
    query_mol.set_group("res1", grp)
    brd_path = os.path.join(tmpdir, "query.mpac")
    bridge.save_msgpack(query_mol.get_raw_data(), brd_path)

    # Matrix for pdf-plot-decaymat
    mat = Matrix(3, 3)
    for i in range(3):
        for j in range(3):
            mat.set(i, j, 1.0 if i == j else 0.1)
    mat_path = os.path.join(tmpdir, "test.mat")
    mat.save(mat_path)

    # ESP and design matrix for pop regression scripts
    mat_design = Matrix(4, 4)
    for i in range(3):
        for j in range(3):
            mat_design.set(i, j, 4.0 if i == j else 1.0)
        mat_design.set(i, 3, 1.0)
        mat_design.set(3, i, 1.0)
    mat_design.set(3, 3, 0.0)
    mat_design_path = os.path.join(tmpdir, "esp_design.mat")
    mat_design.save(mat_design_path)

    vtr = Vector([0.75, -0.375, -0.375, 0.0])
    vtr_path = os.path.join(tmpdir, "esp_target.vtr")
    vtr.save(vtr_path)

    atoms, _ = pdf.PopUtils._get_atomlist(param.molecule)
    real_atoms = [a for a in atoms if a.symbol != "X"]
    q_true = [0.25, -0.125, -0.125]
    for i in range(3):
        real_atoms[i].charge = q_true[i]
    grid_pt = [0.0, 0.0, 2.0]
    pos = bridge.Position(grid_pt[0], grid_pt[1], grid_pt[2])
    pos *= 1.0 / 0.5291772108
    esp_utils = pdf.PopEspUtils()
    exact_esp = esp_utils._calc_esp(pos, real_atoms)
    esp_data = {"grids": [grid_pt], "ESP": [exact_esp]}
    esp_path = os.path.join(tmpdir, "grid-esp.mpac")
    bridge.save_msgpack(esp_data, esp_path)

    charge_csv_path = os.path.join(tmpdir, "charges.csv")
    with open(charge_csv_path, "w") as f:
        for q in q_true:
            f.write(f"{q}\n")

    return {
        "db": db_path,
        "h5": h5_path,
        "brd": brd_path,
        "mat": mat_path,
        "mat_design": mat_design_path,
        "vtr": vtr_path,
        "esp": esp_path,
        "charges": charge_csv_path,
    }


class TestH5DashDBasic(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        self.data = _create_synthetic_system(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_dash_d_help_and_default_config(self):
        """Verify -d help mentions HDF5/SQLite and const/default is pdfresults.h5."""
        for script_name in TARGET_SCRIPTS:
            mod = _load_script(script_name)
            sys.argv = [script_name, "-h"]
            out = io.StringIO()
            with redirect_stdout(out):
                try:
                    mod.main()
                except SystemExit:
                    pass
            help_text = out.getvalue()
            self.assertIn(
                ".h5 は HDF5、それ以外は SQLite(従来)",
                help_text,
                msg=f"Help text missing format description in {script_name}",
            )

    def test_info_geom_h5_and_db_error(self):
        """pdf-info-geom.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-info-geom.py")

        sys.argv = ["pdf-info-geom.py", "-d", self.data["db"]]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = ["pdf-info-geom.py", "-d", self.data["h5"]]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertIn("O", out_h5.getvalue())
        self.assertIn("H", out_h5.getvalue())

    def test_info_xyz_h5_and_db_error(self):
        """pdf-info-xyz.py works with .h5, raises ValueError on .db, and exits 1 in CLI."""
        mod = _load_script("pdf-info-xyz.py")

        # In-process: .db raises ValueError
        sys.argv = ["pdf-info-xyz.py", "-d", self.data["db"]]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        # CLI subprocess: exits with code 1 and prints deprecation notice to stderr
        script_path = os.path.join(SCRIPTS_DIR, "pdf-info-xyz.py")
        env = dict(os.environ)
        top_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        env["PYTHONPATH"] = top_dir + (os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else "")
        proc = subprocess.run(
            [sys.executable, script_path, "-d", self.data["db"]],
            capture_output=True,
            text=True,
            env=env,
        )
        self.assertEqual(proc.returncode, 1)
        self.assertIn("pdf-archive-h5.py", proc.stderr)
        self.assertIn(
            "SQLite 形式は廃止されました。pdf-archive-h5.py で pdfresults.h5 を作り直してください",
            proc.stderr,
        )

        # .h5 succeeds
        sys.argv = ["pdf-info-xyz.py", "-d", self.data["h5"]]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertTrue(out_h5.getvalue().startswith("3\n"))

    def test_info_orb_h5_and_db_error(self):
        """pdf-info-orb.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-info-orb.py")

        sys.argv = ["pdf-info-orb.py", "-d", self.data["db"], "0", "1"]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = ["pdf-info-orb.py", "-d", self.data["h5"], "0", "1"]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertIn("0th:", out_h5.getvalue())
        self.assertIn("shell(s;s)", out_h5.getvalue())

    def test_atom_index_h5_and_db_error(self):
        """pdf-atom-index.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-atom-index.py")

        sys.argv = ["pdf-atom-index.py", self.data["brd"], "-d", self.data["db"]]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = ["pdf-atom-index.py", self.data["brd"], "-d", self.data["h5"]]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertIn("0", out_h5.getvalue().strip())

    def test_plot_decaymat_h5_and_db_error(self):
        """pdf-plot-decaymat.py generates PNG for .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-plot-decaymat.py")
        out_png_db = os.path.join(self.tmpdir.name, "decay_db.png")
        out_png_h5 = os.path.join(self.tmpdir.name, "decay_h5.png")

        sys.argv = ["pdf-plot-decaymat.py", self.data["mat"], "-d", self.data["db"], "-o", out_png_db]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = ["pdf-plot-decaymat.py", self.data["mat"], "-d", self.data["h5"], "-o", out_png_h5]
        mod.main()
        self.assertTrue(os.path.exists(out_png_h5))

    def test_pop_ridge_h5_and_db_error(self):
        """pdf-pop-ridge.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-pop-ridge.py")
        csv_db = os.path.join(self.tmpdir.name, "ridge_db.csv")
        csv_h5 = os.path.join(self.tmpdir.name, "ridge_h5.csv")

        sys.argv = [
            "pdf-pop-ridge.py",
            "-d",
            self.data["db"],
            "--design_matrix_path",
            self.data["mat_design"],
            "--target_vector_path",
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--csv",
            csv_db,
        ]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = [
            "pdf-pop-ridge.py",
            "-d",
            self.data["h5"],
            "--design_matrix_path",
            self.data["mat_design"],
            "--target_vector_path",
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--csv",
            csv_h5,
        ]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertTrue(os.path.exists(csv_h5))
        with open(csv_h5) as f_h5:
            content = f_h5.read()
            self.assertIn("O", content)
            self.assertIn("H", content)

    def test_pop_lasso_h5_and_db_error(self):
        """pdf-pop-lasso.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-pop-lasso.py")
        csv_db = os.path.join(self.tmpdir.name, "lasso_db.csv")
        csv_h5 = os.path.join(self.tmpdir.name, "lasso_h5.csv")

        sys.argv = [
            "pdf-pop-lasso.py",
            "-d",
            self.data["db"],
            self.data["mat_design"],
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--csv",
            csv_db,
        ]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = [
            "pdf-pop-lasso.py",
            "-d",
            self.data["h5"],
            self.data["mat_design"],
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--csv",
            csv_h5,
        ]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertTrue(os.path.exists(csv_h5))
        with open(csv_h5) as f_h5:
            content = f_h5.read()
            self.assertIn("O", content)
            self.assertIn("H", content)

    def test_pop_elasticnet_h5_and_db_error(self):
        """pdf-pop-elasticnet.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-pop-elasticnet.py")
        csv_db = os.path.join(self.tmpdir.name, "enet_db.csv")
        csv_h5 = os.path.join(self.tmpdir.name, "enet_h5.csv")

        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-d",
            self.data["db"],
            self.data["mat_design"],
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "-m",
            "1000",
            "-t",
            "0.001",
            "--csv",
            csv_db,
        ]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = [
            "pdf-pop-elasticnet.py",
            "-d",
            self.data["h5"],
            self.data["mat_design"],
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "-m",
            "1000",
            "-t",
            "0.001",
            "--csv",
            csv_h5,
        ]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertTrue(os.path.exists(csv_h5))
        with open(csv_h5) as f_h5:
            content = f_h5.read()
            self.assertIn("O", content)
            self.assertIn("H", content)

    def test_pop_resp_h5_and_db_error(self):
        """pdf-pop-resp.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-pop-resp.py")
        csv_db = os.path.join(self.tmpdir.name, "resp_db.csv")
        csv_h5 = os.path.join(self.tmpdir.name, "resp_h5.csv")

        sys.argv = [
            "pdf-pop-resp.py",
            "-d",
            self.data["db"],
            "--design_matrix_path",
            self.data["mat_design"],
            "--target_vector_path",
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--csv",
            csv_db,
        ]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = [
            "pdf-pop-resp.py",
            "-d",
            self.data["h5"],
            "--design_matrix_path",
            self.data["mat_design"],
            "--target_vector_path",
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--csv",
            csv_h5,
        ]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertTrue(os.path.exists(csv_h5))
        with open(csv_h5) as f_h5:
            content = f_h5.read()
            self.assertIn("O", content)
            self.assertIn("H", content)

    def test_pop_rrms_h5_and_db_error(self):
        """pdf-pop-rrms.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-pop-rrms.py")

        sys.argv = [
            "pdf-pop-rrms.py",
            "-d",
            self.data["db"],
            "-e",
            self.data["esp"],
            "-c",
            self.data["charges"],
        ]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = [
            "pdf-pop-rrms.py",
            "-d",
            self.data["h5"],
            "-e",
            self.data["esp"],
            "-c",
            self.data["charges"],
        ]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertIn("RRMS", out_h5.getvalue())

    def test_pop_classo_h5_and_db_error(self):
        """pdf-pop-classo.py works with .h5 and raises ValueError on .db."""
        mod = _load_script("pdf-pop-classo.py")
        csv_db = os.path.join(self.tmpdir.name, "classo_db.csv")
        csv_h5 = os.path.join(self.tmpdir.name, "classo_h5.csv")

        sys.argv = [
            "pdf-pop-classo.py",
            "-d",
            self.data["db"],
            "--design_matrix_path",
            self.data["mat_design"],
            "--target_vector_path",
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--max-alpha",
            "0.01",
            "--min-alpha",
            "0.01",
            "--csv",
            csv_db,
        ]
        with self.assertRaises(ValueError) as cm:
            mod.main()
        self.assertIn("pdf-archive-h5.py", str(cm.exception))

        sys.argv = [
            "pdf-pop-classo.py",
            "-d",
            self.data["h5"],
            "--design_matrix_path",
            self.data["mat_design"],
            "--target_vector_path",
            self.data["vtr"],
            "-e",
            self.data["esp"],
            "--max-alpha",
            "0.01",
            "--min-alpha",
            "0.01",
            "--csv",
            csv_h5,
        ]
        out_h5 = io.StringIO()
        with redirect_stdout(out_h5):
            mod.main()

        self.assertTrue(os.path.exists(csv_h5))
        with open(csv_h5) as f_h5:
            content = f_h5.read()
            self.assertIn("O", content)
            self.assertIn("H", content)

    def test_dash_d_default_resolution_from_current_dir(self):
        """Invoking -d without path argument picks up pdfresults.h5 in current dir."""
        default_h5 = os.path.join(self.tmpdir.name, "pdfresults.h5")
        shutil.copyfile(self.data["h5"], default_h5)

        mod = _load_script("pdf-info-geom.py")
        sys.argv = ["pdf-info-geom.py", "-d"]
        out = io.StringIO()
        with redirect_stdout(out):
            mod.main()
        self.assertIn("O", out.getvalue())
        self.assertIn("H", out.getvalue())
