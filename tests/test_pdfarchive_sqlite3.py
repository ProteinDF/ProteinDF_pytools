#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import os
import tempfile
import unittest
import unittest.mock as mock

import proteindf_bridge as bridge
from proteindf_tools import pdfarchive_sqlite3
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfarchive_sqlite3 import PdfArchive_Sqlite3
from proteindf_tools.pdfparam import PdfParam


class TestPdfArchiveSqlite3(unittest.TestCase):
    def setUp(self):
        self.tmpfiles = []
        self.db_path = self._create_temp_db_path()

    def tearDown(self):
        for f in self.tmpfiles:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    def _create_temp_db_path(self):
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as fp:
            path = fp.name
        self.tmpfiles.append(path)
        return path

    def _create_sample_param(self):
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

        param.set_gradient(0, 0.0, 0.0, -0.01)
        param.set_gradient(1, 0.0, 0.0, 0.01)
        return param

    # 1. 構築テスト (未パッチ状態)
    def test_init_without_patch(self):
        # パッチなしで PdfArchive_Sqlite3 を正常に初期化できる。
        # sitecustomize.py による外部パッチの影響を排除してテストする。
        if hasattr(pdfarchive_sqlite3, "PdfArchive"):
            delattr(pdfarchive_sqlite3, "PdfArchive")
        archive = PdfArchive_Sqlite3(self.db_path)
        self.assertIsNotNone(archive)

    # 2. PdfParam からの構築と各種プロパティ・ゲッター
    def test_pdfparam_properties(self):
        param = self._create_sample_param()
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path, pdfparam=param)
            self.assertEqual(archive.num_of_atoms, 2)
            self.assertEqual(archive.num_of_AOs, 2)
            self.assertEqual(archive.num_of_MOs, 2)
            self.assertEqual(archive.method, "rks")
            self.assertEqual(archive.iterations, 2)
            self.assertTrue(archive.scf_converged)

    def test_get_molecule(self):
        param = self._create_sample_param()
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path, pdfparam=param)
            mol = archive.get_molecule()
            self.assertEqual(len(list(mol.atoms())), 2)

            atoms = dict(mol.atoms())
            self.assertIn("0", atoms)
            self.assertIn("1", atoms)

            a0 = atoms["0"]
            self.assertEqual(a0.symbol, "H")
            self.assertEqual(a0.name, "H1")
            self.assertAlmostEqual(a0.xyz.x, 0.0)
            self.assertAlmostEqual(a0.xyz.y, 0.0)
            self.assertAlmostEqual(a0.xyz.z, 0.0)
            self.assertAlmostEqual(a0.charge, 1.0)

            a1 = atoms["1"]
            self.assertEqual(a1.symbol, "H")
            self.assertEqual(a1.name, "H2")
            self.assertAlmostEqual(a1.xyz.x, 0.0)
            self.assertAlmostEqual(a1.xyz.y, 0.0)
            self.assertAlmostEqual(a1.xyz.z, 0.74)
            self.assertAlmostEqual(a1.charge, 1.0)

    def test_get_total_energy(self):
        param = self._create_sample_param()
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path, pdfparam=param)
            self.assertAlmostEqual(archive.get_total_energy(1), -1.10)
            self.assertAlmostEqual(archive.get_total_energy(2), -1.17)
            self.assertIsNone(archive.get_total_energy(99))

    def test_get_gradient_and_rms(self):
        param = self._create_sample_param()
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path, pdfparam=param)
            g0 = archive.get_gradient(0)
            g1 = archive.get_gradient(1)
            self.assertEqual(g0, [0.0, 0.0, -0.01])
            self.assertEqual(g1, [0.0, 0.0, 0.01])

            expected_rms = math.sqrt((0.01**2 + 0.01**2) / (3 * 2))
            self.assertAlmostEqual(archive.get_gradient_rms(), expected_rms, places=6)

    def test_get_basisset(self):
        # 保存された基底関数情報を取得できる。
        param = self._create_sample_param()
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path, pdfparam=param)
            name = archive.get_basisset_name("H@0")
            self.assertEqual(name, "sto-3g-H")
            bs = archive.get_basisset(name)
            self.assertEqual(bs.name, "sto-3g-H")
            self.assertEqual(len(bs), 1)

    # 3. 各種データの往復・HOMO level
    def test_occupations_roundtrip(self):
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path)
            occ = bridge.Vector([2.0, 2.0, 0.0, 0.0])
            archive.set_occupations("rks", occ)

            retrieved = archive.get_occupations("rks")
            self.assertIsNotNone(retrieved)
            self.assertEqual(list(retrieved), [2.0, 2.0, 0.0, 0.0])

    def test_energylevel_roundtrip(self):
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path)
            elevel = bridge.Vector([-0.5, -0.2, 0.1, 0.4])
            archive.set_energylevel("rks", 1, elevel)

            retrieved = archive.get_energylevel("rks", 1)
            self.assertIsNotNone(retrieved)
            self.assertEqual(list(retrieved), [-0.5, -0.2, 0.1, 0.4])

    def test_population_roundtrip(self):
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path)
            m = bridge.Matrix(2, 2)
            m.set(0, 0, 1.0)
            m.set(0, 1, 0.0)
            m.set(1, 0, 0.0)
            m.set(1, 1, 1.0)
            archive.set_population("mulliken", 1, m)

            retrieved = archive.get_population("mulliken", 1)
            self.assertIsNotNone(retrieved)
            self.assertEqual(retrieved.rows, 2)
            self.assertEqual(retrieved.cols, 2)
            self.assertAlmostEqual(retrieved.get(0, 0), 1.0)

    def test_get_homo_level(self):
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path)
            occ = bridge.Vector([2.0, 2.0, 0.0, 0.0])
            archive.set_occupations("rks", occ)
            self.assertEqual(archive.get_HOMO_level("rks"), 1)

    def test_get_homo_level_single_occ(self):
        # 占有軌道が index 0 のみの場合でも HOMO level (0) を正しく取得できる。
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive = PdfArchive_Sqlite3(self.db_path)
            occ = bridge.Vector([2.0, 0.0])
            archive.set_occupations("rks", occ)
            self.assertEqual(archive.get_HOMO_level("rks"), 0)

    # 4. 同じDBの開き直しと compare_* / __eq__
    def test_reopen_database(self):
        param = self._create_sample_param()
        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            archive1 = PdfArchive_Sqlite3(self.db_path, pdfparam=param)
            del archive1

            archive2 = PdfArchive_Sqlite3(self.db_path)
            self.assertEqual(archive2.num_of_atoms, 2)
            self.assertEqual(archive2.num_of_AOs, 2)
            self.assertEqual(archive2.method, "rks")
            self.assertAlmostEqual(archive2.get_total_energy(1), -1.10)
            self.assertAlmostEqual(archive2.get_total_energy(2), -1.17)
            self.assertEqual(archive2.get_gradient(0), [0.0, 0.0, -0.01])

    def test_compare_and_eq(self):
        param1 = self._create_sample_param()
        db_path2 = self._create_temp_db_path()
        param2 = self._create_sample_param()

        with mock.patch.object(pdfarchive_sqlite3, "PdfArchive", PdfArchive_Sqlite3, create=True):
            arc1 = PdfArchive_Sqlite3(self.db_path, pdfparam=param1)
            arc2 = PdfArchive_Sqlite3(db_path2, pdfparam=param2)

            self.assertTrue(arc1.compare_info(arc2))
            self.assertTrue(arc1.compare_energy(arc2))
            self.assertTrue(arc1.compare_grad(arc2))
            self.assertTrue(arc1 == arc2)

            # エネルギーが異なる場合
            db_path3 = self._create_temp_db_path()
            param3 = self._create_sample_param()
            param3.set_total_energy(2, -2.00)
            arc3 = PdfArchive_Sqlite3(db_path3, pdfparam=param3)

            self.assertFalse(arc1.compare_energy(arc3))
            self.assertFalse(arc1 == arc3)


if __name__ == "__main__":
    unittest.main()
