#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

try:
    import h5py

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfcommon import open_archive
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_hdf5 import PdfParam_H5
from proteindf_tools.poputils import PopUtils
from proteindf_tools.vector import Vector


def _build_sample_molecule():
    mol = bridge.AtomGroup()

    atom_o = bridge.Atom(
        symbol="O",
        position=bridge.Position([0.0, 0.5, 1.0]),
        charge=8.0,
        name="O1",
    )
    mol.set_atom("1", atom_o)

    atom_h = bridge.Atom(
        symbol="H",
        position=bridge.Position([0.0, -0.5, 2.0]),
        charge=1.0,
        name="H2",
    )
    mol.set_atom("2", atom_h)

    atom_x = bridge.Atom(
        symbol="X",
        position=bridge.Position([1.0, 1.0, 1.0]),
        charge=0.5,
        name="X3",
    )
    mol.set_atom("3", atom_x)
    return mol


def _create_synthetic_h5(h5_path):
    p = PdfParam_H5()
    p.open(h5_path)
    p.comment = "synthetic h5"
    p.num_of_atoms = 3
    p.num_of_AOs = 5
    p.num_of_MOs = 5
    p.method = "rks"
    p.iterations = 1
    p.max_iterations = 10
    p._set_scf_converged(True)
    p.molecule = _build_sample_molecule()
    p.save_model()

    p.save_TEs(Vector([-75.12345]))
    p.save_occ("rks", Vector([2.0, 2.0, 1.0, 0.0, 0.0]))
    p.save_energy_level("rks", 1, Vector([-19.1, -1.2, -0.5, 0.3, 0.9]))
    p._close()


def _create_synthetic_db(db_path):
    param = PdfParam()
    param.comment = "synthetic db"
    param.num_of_atoms = 3
    param.num_of_AOs = 5
    param.num_of_MOs = 5
    param.method = "rks"
    param.iterations = 1
    param.max_iterations = 10
    param.guess = "none"
    param._set_scf_converged(True)
    param.molecule = _build_sample_molecule()
    param.set_total_energy(1, -75.12345)

    db = PdfArchive_Sqlite3(db_path, param)
    db.set_occupations("rks", Vector([2.0, 2.0, 1.0, 0.0, 0.0]))
    db.set_energylevel("rks", 1, Vector([-19.1, -1.2, -0.5, 0.3, 0.9]))


class TestOpenArchive(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    @unittest.skipUnless(HAVE_H5PY, "h5py not available")
    def test_open_archive_dispatch_by_extension(self):
        """Verify extension-based dispatch (.h5, .hdf5, .H5, .HDF5 accept, .db and others raise ValueError)."""
        for filename in ("test.h5", "test.hdf5", "test.H5", "test.HDF5"):
            path = os.path.join(self.tmpdir.name, filename)
            entry = open_archive(path)
            try:
                self.assertIsInstance(entry, PdfParam_H5, msg=f"Failed for {filename}")
            finally:
                entry._close()

        for filename in ("test.db", "test.sqlite", "pdfentry"):
            path = os.path.join(self.tmpdir.name, filename)
            with self.assertRaises(ValueError) as cm:
                open_archive(path)
            self.assertIn("SQLite 形式は廃止されました。pdf-archive-h5.py で pdfresults.h5 を作り直してください", str(cm.exception))
            self.assertIn("pdf-archive-h5.py", str(cm.exception))

        # Also exposed on top-level proteindf_tools module
        self.assertIs(pdf.open_archive, open_archive)

    @unittest.skipUnless(HAVE_H5PY, "h5py not available")
    def test_open_archive_h5_compat_methods(self):
        """Verify common and alias methods on synthetic HDF5 archive."""
        h5_path = os.path.join(self.tmpdir.name, "sample.h5")
        _create_synthetic_h5(h5_path)

        entry = open_archive(h5_path)
        try:
            self.assertIsInstance(entry, PdfParam_H5)

            mol_prop = entry.molecule
            mol_meth = entry.get_molecule()
            self.assertEqual(mol_prop.get_number_of_atoms(), 3)
            self.assertEqual(mol_meth.get_number_of_atoms(), 3)

            atoms_prop = mol_prop.get_atom_list()
            atoms_meth = mol_meth.get_atom_list()
            for a1, a2 in zip(atoms_prop, atoms_meth):
                self.assertEqual(a1.symbol, a2.symbol)
                self.assertEqual(a1.name, a2.name)
                self.assertAlmostEqual(a1.xyz.x, a2.xyz.x)
                self.assertAlmostEqual(a1.xyz.y, a2.xyz.y)
                self.assertAlmostEqual(a1.xyz.z, a2.xyz.z)
                self.assertAlmostEqual(a1.charge, a2.charge)

            el1 = entry.get_energy_level("rks", 1)
            el2 = entry.get_energylevel("rks", 1)
            expected_el = [-19.1, -1.2, -0.5, 0.3, 0.9]
            self.assertEqual(len(el1), len(expected_el))
            self.assertEqual(len(el2), len(expected_el))
            for v1, v2, exp in zip(el1, el2, expected_el):
                self.assertAlmostEqual(v1, exp)
                self.assertAlmostEqual(v2, exp)

            self.assertEqual(entry.get_HOMO_level("rks"), 2)
            self.assertAlmostEqual(entry.get_total_energy(1), -75.12345)
        finally:
            entry._close()

    def test_open_archive_sqlite_raises_value_error(self):
        """Verify open_archive raises ValueError with informative message for .db."""
        db_path = os.path.join(self.tmpdir.name, "sample.db")
        _create_synthetic_db(db_path)

        with self.assertRaises(ValueError) as cm:
            open_archive(db_path)
        msg = str(cm.exception)
        self.assertIn("SQLite 形式は廃止されました。pdf-archive-h5.py で pdfresults.h5 を作り直してください", msg)
        self.assertIn("pdf-archive-h5.py", msg)
        self.assertIn(db_path, msg)

    @unittest.skipUnless(HAVE_H5PY, "h5py not available")
    def test_poputils_get_atomlist_by_db(self):
        """Verify PopUtils.get_atomlist_by_db works with HDF5 and raises ValueError on .db."""
        h5_path = os.path.join(self.tmpdir.name, "pop.h5")
        db_path = os.path.join(self.tmpdir.name, "pop.db")
        _create_synthetic_h5(h5_path)
        _create_synthetic_db(db_path)

        atoms_h5, charge_h5 = PopUtils.get_atomlist_by_db(h5_path)

        # Dummy atom 'X' should be filtered out, leaving O and H
        self.assertEqual(len(atoms_h5), 2)
        self.assertEqual([a.symbol for a in atoms_h5], ["O", "H"])
        self.assertAlmostEqual(charge_h5, -0.5)

        with self.assertRaises(ValueError) as cm:
            PopUtils.get_atomlist_by_db(db_path)
        self.assertIn("pdf-archive-h5.py", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
