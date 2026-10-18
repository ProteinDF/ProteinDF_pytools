#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _create_sample_param():
    """Create minimal synthetic PdfParam representing water (H2O: O, H, H)."""
    param = PdfParam()
    param.comment = "synthetic water system"
    param.num_of_atoms = 3
    param.num_of_AOs = 5
    param.num_of_MOs = 5
    param.method = "rks"
    param.iterations = 1
    param.guess = "none"
    param._set_scf_converged(True)

    mol = bridge.AtomGroup()
    a_o = bridge.Atom(
        symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="O1"
    )
    a_h1 = bridge.Atom(
        symbol="H", position=bridge.Position([0.0, 1.0, 0.0]), charge=1.0, name="H1"
    )
    a_h2 = bridge.Atom(
        symbol="H", position=bridge.Position([1.0, 0.0, 0.0]), charge=1.0, name="H2"
    )
    mol.set_atom("0", a_o)
    mol.set_atom("1", a_h1)
    mol.set_atom("2", a_h2)
    param.molecule = mol

    # BasisSet for O: 1 s shell
    bs_o = BasisSet()
    bs_o.name = "sto-3g-O"
    cgto_o = ContractedGTO("s", 1)
    cgto_o[0] = PrimitiveGTO(71.6168, 0.15433)
    bs_o.append(cgto_o)

    # BasisSet for H: 1 s shell
    bs_h = BasisSet()
    bs_h.name = "sto-3g-H"
    cgto_h = ContractedGTO("s", 1)
    cgto_h[0] = PrimitiveGTO(3.42525, 0.15433)
    bs_h.append(cgto_h)

    param.set_basisset("O", bs_o)
    param.set_basisset("H", bs_h)
    param.set_total_energy(1, -75.0)

    param._data.setdefault("control", {})["file_base_name"] = {
        "eigenvalues": "eigenvalues.%s.vtr",
        "occupation_vtr": "occ.%s.vtr",
        "Spq_matrix": "s.mat",
        "Hpq_matrix": "h.mat",
        "Hpq2_matrix": "h2.mat",
        "Fpq_matrix": "f.%s.mat",
        "C_matrix": "c.%s.mat",
        "Ppq_matrix": "p.%s.mat",
    }
    return param


class TestPdfInfoOrb(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.param = _create_sample_param()
        self.mpac_path = os.path.join(self.tmpdir.name, "water.mpac")
        pdf.save_pdfparam(self.param, self.mpac_path)
        self.db_path = os.path.join(self.tmpdir.name, "water.db")
        self.archive = PdfArchive_Sqlite3(self.db_path, pdfparam=self.param)
        self.mod = _load_script("pdf-info-orb.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_info_orb_param_success(self):
        """Orbital information is output from param file."""
        sys.argv = ["pdf-info-orb.py", "-p", self.mpac_path, "0", "1"]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("0th:", output)
        self.assertIn("1th:", output)
        self.assertIn("shell(s;s)", output)

    def test_info_orb_db_success(self):
        """Orbital information is output from db file."""
        sys.argv = ["pdf-info-orb.py", "-d", self.db_path, "0"]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("0th:", output)
        self.assertIn("shell(s;s)", output)

    def test_info_orb_out_of_range(self):
        """Out of range orbital index produces no output."""
        sys.argv = ["pdf-info-orb.py", "-p", self.mpac_path, "99"]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        self.assertEqual(out.getvalue().strip(), "")

    def test_info_orb_both_options_fails(self):
        """Providing both -d and -p raises SystemExit."""
        sys.argv = ["pdf-info-orb.py", "-d", self.db_path, "-p", self.mpac_path, "0"]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_info_orb_neither_option_fails(self):
        """Providing neither -d nor -p raises SystemExit."""
        sys.argv = ["pdf-info-orb.py", "0"]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_info_orb_missing_orb_index_fails(self):
        """Missing orbital index argument raises SystemExit."""
        sys.argv = ["pdf-info-orb.py", "-p", self.mpac_path]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_info_orb_nonexistent_param(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = ["pdf-info-orb.py", "-p", os.path.join(self.tmpdir.name, "missing.mpac"), "0"]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()


class TestPdfInfoGeom(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.param = _create_sample_param()
        self.mpac_path = os.path.join(self.tmpdir.name, "water.mpac")
        pdf.save_pdfparam(self.param, self.mpac_path)
        self.db_path = os.path.join(self.tmpdir.name, "water.db")
        self.archive = PdfArchive_Sqlite3(self.db_path, pdfparam=self.param)
        self.mod = _load_script("pdf-info-geom.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_info_geom_param_success(self):
        """Coordinates of synthetic molecule are printed using -p."""
        sys.argv = ["pdf-info-geom.py", "-p", self.mpac_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("O", output)
        self.assertIn("H", output)

    def test_info_geom_db_success(self):
        """Coordinates of synthetic molecule are printed using -d."""
        sys.argv = ["pdf-info-geom.py", "-d", self.db_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue()
        self.assertIn("O", output)
        self.assertIn("H", output)

    def test_info_geom_no_options_fails(self):
        """Omitting required options raises SystemExit."""
        sys.argv = ["pdf-info-geom.py"]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_info_geom_nonexistent_param(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = ["pdf-info-geom.py", "-p", os.path.join(self.tmpdir.name, "missing.mpac")]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()


class TestPdfInfoXyz(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.param = _create_sample_param()
        self.mpac_path = os.path.join(self.tmpdir.name, "water.mpac")
        pdf.save_pdfparam(self.param, self.mpac_path)
        self.db_path = os.path.join(self.tmpdir.name, "water.db")
        self.archive = PdfArchive_Sqlite3(self.db_path, pdfparam=self.param)
        self.mod = _load_script("pdf-info-xyz.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_info_xyz_param_success(self):
        """XYZ format output is printed using -p."""
        sys.argv = ["pdf-info-xyz.py", "-p", self.mpac_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        lines = out.getvalue().strip().splitlines()
        self.assertGreaterEqual(len(lines), 3)
        self.assertEqual(lines[0].strip(), "3")
        self.assertTrue(any(line.strip().startswith("O") for line in lines[2:]))
        self.assertTrue(any(line.strip().startswith("H") for line in lines[2:]))

    def test_info_xyz_db_success(self):
        """XYZ format output is printed using -d."""
        sys.argv = ["pdf-info-xyz.py", "-d", self.db_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        lines = out.getvalue().strip().splitlines()
        self.assertGreaterEqual(len(lines), 3)
        self.assertEqual(lines[0].strip(), "3")
        self.assertTrue(any(line.strip().startswith("O") for line in lines[2:]))
        self.assertTrue(any(line.strip().startswith("H") for line in lines[2:]))

    def test_info_xyz_no_options_fails(self):
        """Omitting required options raises SystemExit."""
        sys.argv = ["pdf-info-xyz.py"]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_info_xyz_nonexistent_param(self):
        """Nonexistent parameter file raises FileNotFoundError."""
        sys.argv = ["pdf-info-xyz.py", "-p", os.path.join(self.tmpdir.name, "missing.mpac")]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()


class TestPdfAtomIndex(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.param = _create_sample_param()
        self.mpac_path = os.path.join(self.tmpdir.name, "water.mpac")
        pdf.save_pdfparam(self.param, self.mpac_path)
        self.db_path = os.path.join(self.tmpdir.name, "water.db")
        self.archive = PdfArchive_Sqlite3(self.db_path, pdfparam=self.param)

        # Create synthetic query AtomGroup with a sub-group
        self.query_mol = bridge.AtomGroup()
        grp = bridge.AtomGroup()
        target_atom = bridge.Atom(
            symbol="O", position=bridge.Position([0.0, 0.0, 0.0]), charge=8.0, name="O1"
        )
        grp.set_atom("0", target_atom)
        self.query_mol.set_group("residue1", grp)
        self.brd_path = os.path.join(self.tmpdir.name, "query.mpac")
        bridge.save_msgpack(self.query_mol.get_raw_data(), self.brd_path)

        self.mod = _load_script("pdf-atom-index.py")
        self._orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_atom_index_param_stdout(self):
        """Atom index is found and output to stdout using -p."""
        sys.argv = ["pdf-atom-index.py", self.brd_path, "-p", self.mpac_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue().strip()
        self.assertIn("0", output)

    def test_atom_index_param_output_file(self):
        """Atom index is found and saved to msgpack file specified by -o."""
        out_file = os.path.join(self.tmpdir.name, "indices.mpac")
        sys.argv = ["pdf-atom-index.py", self.brd_path, "-p", self.mpac_path, "-o", out_file]
        self.mod.main()
        self.assertTrue(os.path.exists(out_file))
        loaded = bridge.load_msgpack(out_file)
        self.assertEqual(loaded, [[0]])

    def test_atom_index_db_stdout(self):
        """Atom index is found using -d from database.

        scripts/pdf-atom-index.py:70: -d option does not define pdfparam
        (UnboundLocalError).
        """
        sys.argv = ["pdf-atom-index.py", self.brd_path, "-d", self.db_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue().strip()
        self.assertIn("0", output)

    def test_atom_index_atom_not_found(self):
        """Querying an atom not present in the system returns index -1."""
        missing_mol = bridge.AtomGroup()
        grp = bridge.AtomGroup()
        unknown_atom = bridge.Atom(
            symbol="Fe", position=bridge.Position([10.0, 10.0, 10.0]), charge=26.0, name="Fe1"
        )
        grp.set_atom("0", unknown_atom)
        missing_mol.set_group("residue1", grp)
        brd_path = os.path.join(self.tmpdir.name, "missing_atom.mpac")
        bridge.save_msgpack(missing_mol.get_raw_data(), brd_path)

        sys.argv = ["pdf-atom-index.py", brd_path, "-p", self.mpac_path]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue().strip()
        self.assertIn("-1", output)

    def test_atom_index_no_options_fails(self):
        """Omitting required -d or -p option raises SystemExit."""
        sys.argv = ["pdf-atom-index.py", self.brd_path]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_atom_index_nonexistent_brd(self):
        """Nonexistent brd path raises FileNotFoundError."""
        sys.argv = [
            "pdf-atom-index.py",
            os.path.join(self.tmpdir.name, "missing.brd"),
            "-p",
            self.mpac_path,
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()


class TestPdfInfoBasisset(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.mod = _load_script("pdf-info-basisset.py")
        self._orig_argv = sys.argv
        self._orig_pdf_home = os.environ.get("PDF_HOME")

        # Create synthetic basis2 file in tmpdir/data/basis2
        data_dir = os.path.join(self.tmpdir.name, "data")
        os.makedirs(data_dir, exist_ok=True)
        self.basis2_path = os.path.join(data_dir, "basis2")
        with open(self.basis2_path, "w") as f:
            f.write(
                "O-SYNTHETIC\n"
                " 1 0 0 0 0 0 0 0 0 0 0\n"
                " 1\n"
                "     1.000000e+00  1.000000e+00\n"
            )

    def tearDown(self):
        sys.argv = self._orig_argv
        if self._orig_pdf_home is not None:
            os.environ["PDF_HOME"] = self._orig_pdf_home
        else:
            os.environ.pop("PDF_HOME", None)
        self.tmpdir.cleanup()

    def test_info_basisset_success(self):
        """Known basis set information is printed when PDF_HOME points to basis data."""
        os.environ["PDF_HOME"] = self.tmpdir.name
        sys.argv = ["pdf-info-basisset.py", "O-SYNTHETIC"]
        out = io.StringIO()
        with redirect_stdout(out):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 0)
        output = out.getvalue()
        self.assertIn("O-SYNTHETIC", output)
        self.assertIn("1.000000e+00", output)

    def test_info_basisset_no_pdf_home_fails(self):
        """When PDF_HOME is unset or invalid, FileNotFoundError is raised."""
        os.environ.pop("PDF_HOME", None)
        sys.argv = ["pdf-info-basisset.py", "O-SYNTHETIC"]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()

    def test_info_basisset_nonexistent_name(self):
        """Nonexistent basisset name prints empty BasisSet."""
        os.environ["PDF_HOME"] = self.tmpdir.name
        sys.argv = ["pdf-info-basisset.py", "NONEXISTENT"]
        out = io.StringIO()
        with redirect_stdout(out):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 0)
        output = out.getvalue()
        self.assertNotIn("NONEXISTENT", output)

    def test_info_basisset_no_args_fails(self):
        """Omitting basisset name argument raises SystemExit."""
        sys.argv = ["pdf-info-basisset.py"]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)


class TestPdfEstimateAo(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.mod = _load_script("pdf-estimate-ao.py")
        self._orig_argv = sys.argv
        self._orig_pdf_home = os.environ.get("PDF_HOME")

        # Create synthetic basis2 file in tmpdir/data/basis2
        data_dir = os.path.join(self.tmpdir.name, "data")
        os.makedirs(data_dir, exist_ok=True)
        self.basis2_path = os.path.join(data_dir, "basis2")
        with open(self.basis2_path, "w") as f:
            f.write(
                "O-MINI.H\n"
                " 1 0 0 0 0 0 0 0 0 0 0\n"
                " 1\n"
                "     1.000000e+00  1.000000e+00\n"
                "O-MINI.O\n"
                " 2 1 0 0 0 0 0 0 0 0 0\n"
                " 1\n"
                "     2.000000e+00  1.000000e+00\n"
                " 1\n"
                "     3.000000e+00  1.000000e+00\n"
                " 1\n"
                "     4.000000e+00  1.000000e+00\n"
            )

        # Create synthetic molecule (H: 1 s -> 1 AO, O: 2 s + 1 p -> 2 + 3 = 5 AOs; total = 6)
        mol = bridge.AtomGroup()
        h = bridge.Atom(symbol="H", position=bridge.Position([0.0, 0.0, 0.0]))
        o = bridge.Atom(symbol="O", position=bridge.Position([0.0, 1.0, 0.0]))
        mol.set_atom("0", h)
        mol.set_atom("1", o)
        self.mpac_path = os.path.join(self.tmpdir.name, "sample_mol.mpac")
        bridge.save_msgpack(mol.get_raw_data(), self.mpac_path)

    def tearDown(self):
        sys.argv = self._orig_argv
        if self._orig_pdf_home is not None:
            os.environ["PDF_HOME"] = self._orig_pdf_home
        else:
            os.environ.pop("PDF_HOME", None)
        self.tmpdir.cleanup()

    def test_estimate_ao_success(self):
        """Number of AOs is correctly estimated from synthetic molecule and basis2."""
        os.environ["PDF_HOME"] = self.tmpdir.name
        sys.argv = ["pdf-estimate-ao.py", self.mpac_path, "MINI"]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue().strip()
        self.assertEqual(output, "6")

    def test_estimate_ao_verbose_defect(self):
        """Estimating AOs with -v verbose option succeeds without ValueError.

        scripts/pdf-estimate-ao.py:36: print("reading: {}}".format(...)) has
        an extra closing brace causing ValueError.
        """
        os.environ["PDF_HOME"] = self.tmpdir.name
        sys.argv = ["pdf-estimate-ao.py", self.mpac_path, "MINI", "-v"]
        out = io.StringIO()
        with redirect_stdout(out):
            self.mod.main()
        output = out.getvalue().strip()
        self.assertIn("6", output)

    def test_estimate_ao_no_args_fails(self):
        """Omitting required arguments raises SystemExit."""
        sys.argv = ["pdf-estimate-ao.py"]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                self.mod.main()
        self.assertEqual(cm.exception.code, 2)

    def test_estimate_ao_nonexistent_file(self):
        """Nonexistent input file raises FileNotFoundError."""
        os.environ["PDF_HOME"] = self.tmpdir.name
        sys.argv = [
            "pdf-estimate-ao.py",
            os.path.join(self.tmpdir.name, "missing.mpac"),
            "MINI",
        ]
        with self.assertRaises(FileNotFoundError):
            self.mod.main()


if __name__ == "__main__":
    unittest.main()
