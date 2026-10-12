#!/usr/bin/env python
# -*- coding: utf-8 -*-

import contextlib
import io
import math
import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, call, patch

import proteindf_bridge as bridge
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfsim import PdfSim


class TestPdfSim(unittest.TestCase):
    def setUp(self):
        self.orig_cwd = os.getcwd()
        self.orig_argv = list(sys.argv)

    def tearDown(self):
        os.chdir(self.orig_cwd)
        sys.argv = self.orig_argv

    # ------------------------------------------------------------------
    # 1. setup
    # ------------------------------------------------------------------
    def test_setup_creates_files_and_directories(self):
        sim = PdfSim()
        param = PdfParam()
        param.get_inputfile_contents = lambda: "dummy input content\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = os.path.join(tmpdir, "work")
            sim.setup(param, workdir)

            input_path = os.path.join(workdir, "fl_Userinput")
            work_sub_path = os.path.join(workdir, "fl_Work")

            self.assertTrue(os.path.exists(input_path))
            with open(input_path, "r") as f:
                self.assertEqual(f.read(), "dummy input content\n")
            self.assertTrue(os.path.isdir(work_sub_path))

    def test_setup_existing_workdir(self):
        sim = PdfSim()
        param = PdfParam()
        param.get_inputfile_contents = lambda: "updated input content\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = os.path.join(tmpdir, "existing_work")
            os.mkdir(workdir)
            os.mkdir(os.path.join(workdir, "fl_Work"))

            # 既存ディレクトリが存在していてもエラーにならずに実行できる
            sim.setup(param, workdir)

            input_path = os.path.join(workdir, "fl_Userinput")
            self.assertTrue(os.path.exists(input_path))
            with open(input_path, "r") as f:
                self.assertEqual(f.read(), "updated input content\n")

    def test_setup_invalid_pdfparam_raises_assertion_error(self):
        sim = PdfSim()
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(AssertionError):
                sim.setup("not_pdfparam", tmpdir)
            with self.assertRaises(AssertionError):
                sim.setup(12345, tmpdir)
            with self.assertRaises(AssertionError):
                sim.setup(object(), tmpdir)

    # ------------------------------------------------------------------
    # 2. sp
    # ------------------------------------------------------------------
    @patch("proteindf_tools.pdfsim.run_pdf")
    def test_sp_dry_run(self, mock_run_pdf):
        sim = PdfSim()
        param = PdfParam()
        param.get_inputfile_contents = lambda: "content"

        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [os.path.join(tmpdir, "dummy_script.py")]
            os.chdir(tmpdir)
            cwd_before = os.getcwd()

            ret = sim.sp(param, workdir="calc_dry", dry_run=True)

            self.assertEqual(ret, (None, None))
            mock_run_pdf.assert_not_called()
            self.assertEqual(os.getcwd(), cwd_before)

    @patch("proteindf_tools.pdfsim.run_pdf")
    def test_sp_run_pdf_calls_and_cwd_restored(self, mock_run_pdf):
        sim = PdfSim()
        param = PdfParam()
        param.get_inputfile_contents = lambda: "content"

        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [os.path.join(tmpdir, "dummy_script.py")]
            os.chdir(tmpdir)
            cwd_before = os.getcwd()

            ret = sim.sp(param, workdir="calc_sp", dry_run=False)

            mock_run_pdf.assert_has_calls([call(["serial", "-o", "pdf.log"]), call("archive")])
            self.assertEqual(mock_run_pdf.call_count, 2)
            self.assertEqual(ret, (None, None))
            self.assertEqual(os.getcwd(), cwd_before)

    @patch("proteindf_tools.pdfsim.PdfParam_H5")
    @patch("proteindf_tools.pdfsim.run_pdf")
    def test_sp_archive_h5_mock_and_cwd_restored(self, mock_run_pdf, mock_h5_cls):
        sim = PdfSim()
        param = PdfParam()
        param.get_inputfile_contents = lambda: "content"

        mock_h5 = mock_h5_cls.return_value
        mock_h5.iterations = 20
        mock_h5.get_total_energy.return_value = -123.456

        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [os.path.join(tmpdir, "dummy_script.py")]
            os.chdir(tmpdir)
            cwd_before = os.getcwd()

            ret = sim.sp(
                param,
                workdir="calc_h5",
                cmd_pdf="parallel",
                cmd_archive="archive-h5",
                dry_run=False,
            )

            mock_run_pdf.assert_has_calls(
                [call(["parallel", "-o", "pdf.log"]), call("archive-h5")]
            )
            mock_h5.open.assert_called_once_with("pdfresults.h5")
            mock_h5.get_total_energy.assert_called_once_with(20)
            self.assertEqual(ret, (20, -123.456))
            self.assertEqual(os.getcwd(), cwd_before)

    # ------------------------------------------------------------------
    # 3. pop
    # ------------------------------------------------------------------
    @patch("proteindf_tools.pdfsim.run_pdf")
    def test_pop_dry_run(self, mock_run_pdf):
        sim = PdfSim()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            cwd_before = os.getcwd()

            sim.pop(iteration=5, workdir="pop_dry", dry_run=True)

            mock_run_pdf.assert_not_called()
            self.assertEqual(os.getcwd(), cwd_before)

    @patch("proteindf_tools.pdfsim.run_pdf")
    def test_pop_without_iteration(self, mock_run_pdf):
        sim = PdfSim()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            cwd_before = os.getcwd()

            sim.pop(workdir="pop_default", dry_run=False)

            mock_run_pdf.assert_called_once_with(["pop-mulliken"])
            self.assertEqual(os.getcwd(), cwd_before)

    @patch("proteindf_tools.pdfsim.run_pdf")
    def test_pop_with_iteration(self, mock_run_pdf):
        sim = PdfSim()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            cwd_before = os.getcwd()

            sim.pop(iteration=7, workdir="pop_iter", dry_run=False)

            mock_run_pdf.assert_called_once_with(["pop-mulliken", "-i", 7])
            self.assertEqual(os.getcwd(), cwd_before)

    # ------------------------------------------------------------------
    # 4. _move
    # ------------------------------------------------------------------
    def test_move_directions(self):
        sim = PdfSim()
        atom = bridge.Atom()
        atom.xyz.x = 1.0
        atom.xyz.y = 2.0
        atom.xyz.z = 3.0

        moved_x = sim._move(atom, 0, 0.5)
        self.assertAlmostEqual(moved_x.xyz.x, 1.5)
        self.assertAlmostEqual(moved_x.xyz.y, 2.0)
        self.assertAlmostEqual(moved_x.xyz.z, 3.0)

        moved_y = sim._move(atom, 1, -0.25)
        self.assertAlmostEqual(moved_y.xyz.x, 1.0)
        self.assertAlmostEqual(moved_y.xyz.y, 1.75)
        self.assertAlmostEqual(moved_y.xyz.z, 3.0)

        moved_z = sim._move(atom, 2, 0.75)
        self.assertAlmostEqual(moved_z.xyz.x, 1.0)
        self.assertAlmostEqual(moved_z.xyz.y, 2.0)
        self.assertAlmostEqual(moved_z.xyz.z, 3.75)

        # 元の原子は変更されない
        self.assertAlmostEqual(atom.xyz.x, 1.0)
        self.assertAlmostEqual(atom.xyz.y, 2.0)
        self.assertAlmostEqual(atom.xyz.z, 3.0)

    def test_move_invalid_direction_raises_assertion_error(self):
        sim = PdfSim()
        atom = bridge.Atom()
        with self.assertRaises(AssertionError):
            sim._move(atom, -1, 0.1)
        with self.assertRaises(AssertionError):
            sim._move(atom, 3, 0.1)
        with self.assertRaises(AssertionError):
            sim._move(atom, 4, 0.1)

    # ------------------------------------------------------------------
    # 5. numerical_grad
    # ------------------------------------------------------------------
    def test_numerical_grad_with_mock_sp(self):
        sim = PdfSim()
        param = PdfParam()
        mol = bridge.AtomGroup()
        atom = bridge.Atom()
        atom.symbol = "H"
        atom.label = "1"
        atom.xyz.x = 1.0
        atom.xyz.y = 2.0
        atom.xyz.z = 3.0
        mol.set_atom("1", atom)
        param.molecule = mol

        # 解析的エネルギー関数: E(x, y, z) = x^2 + 2*y^2 + 3*z^2
        # 解析勾配: dE/dx = 2*x = 2.0, dE/dy = 4*y = 8.0, dE/dz = 6*z = 18.0
        # rms = sqrt((2^2 + 8^2 + 18^2) / 3) = sqrt(392 / 3)
        def mock_sp(p, workdir):
            a = p.molecule["1"]
            x = float(a.xyz.x)
            y = float(a.xyz.y)
            z = float(a.xyz.z)
            energy = x**2 + 2.0 * y**2 + 3.0 * z**2
            return (1, energy)

        sim.sp = mock_sp

        with tempfile.TemporaryDirectory() as tmpdir:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                grad_mat, rms = sim.numerical_grad(
                    param, workdir=tmpdir, accuracy=1.0e-3, delta=0.001
                )

            self.assertEqual(len(grad_mat), 1)
            self.assertAlmostEqual(grad_mat[0][0], 2.0, places=4)
            self.assertAlmostEqual(grad_mat[0][1], 8.0, places=4)
            self.assertAlmostEqual(grad_mat[0][2], 18.0, places=4)
            self.assertAlmostEqual(rms, math.sqrt(392.0 / 3.0), places=4)

            # numerical_grad の呼び出し前後で全原子の座標が保持されること
            self.assertAlmostEqual(param.molecule["1"].xyz.x, 1.0, delta=1e-12)
            self.assertAlmostEqual(param.molecule["1"].xyz.y, 2.0, delta=1e-12)
            self.assertAlmostEqual(param.molecule["1"].xyz.z, 3.0, delta=1e-12)

    def test_numerical_grad_preserves_molecule_coordinates(self):
        """sp をモックした解析的エネルギーで numerical_grad を呼び、
        呼び出し前後で pdfparam.molecule の全原子の座標が許容誤差 1e-12 で変わらないこと。"""
        sim = PdfSim()
        param = PdfParam()
        mol = bridge.AtomGroup()

        atom1 = bridge.Atom()
        atom1.symbol = "H"
        atom1.label = "1"
        atom1.xyz.x = 1.0
        atom1.xyz.y = 2.0
        atom1.xyz.z = 3.0
        mol.set_atom("1", atom1)

        atom2 = bridge.Atom()
        atom2.symbol = "O"
        atom2.label = "2"
        atom2.xyz.x = -1.5
        atom2.xyz.y = 0.5
        atom2.xyz.z = 2.5
        mol.set_atom("2", atom2)

        param.molecule = mol

        def mock_sp(p, workdir):
            energy = 0.0
            for atom_id, atom in p.molecule.atoms():
                energy += float(atom.xyz.x) ** 2 + float(atom.xyz.y) ** 2 + float(atom.xyz.z) ** 2
            return (1, energy)

        sim.sp = mock_sp

        coords_before = {
            atom_id: (float(atom.xyz.x), float(atom.xyz.y), float(atom.xyz.z))
            for atom_id, atom in param.molecule.atoms()
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                sim.numerical_grad(param, workdir=tmpdir, accuracy=1.0e-3, delta=0.001)

        coords_after = {
            atom_id: (float(atom.xyz.x), float(atom.xyz.y), float(atom.xyz.z))
            for atom_id, atom in param.molecule.atoms()
        }

        self.assertEqual(set(coords_before.keys()), set(coords_after.keys()))
        for atom_id in coords_before:
            for c_before, c_after in zip(coords_before[atom_id], coords_after[atom_id]):
                self.assertAlmostEqual(c_before, c_after, delta=1e-12)

    # ------------------------------------------------------------------
    # 6. opt と total_energy
    # ------------------------------------------------------------------
    def test_opt_existing_workdir(self):
        """opt で既存の作業ディレクトリが存在する場合でも正常に opt が完了すること。"""
        sim = PdfSim()
        param = PdfParam()
        mol = bridge.AtomGroup()
        atom = bridge.Atom()
        atom.symbol = "H"
        atom.label = "1"
        atom.xyz.x = 0.0
        atom.xyz.y = 0.0
        atom.xyz.z = 0.0
        mol.set_atom("1", atom)
        param.molecule = mol

        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [os.path.join(tmpdir, "main.py")]
            os.mkdir(os.path.join(tmpdir, "opt_cycle_1"))

            sim.numerical_grad = MagicMock(return_value=([[0.0, 0.0, 0.0]], 0.0))
            # 期待値: 既存ディレクトリが存在していてもエラーにならず opt が完了すること
            sim.opt(param, max_cycle=1)

    def test_opt_multi_atom_index_advances_and_moves_opposite_to_gradient(self):
        """opt の更新ループで原子ごとの添字 i が進み、各原子が自身の勾配の逆向き (-alpha * grad) に更新されること。"""
        sim = PdfSim()
        param = PdfParam()
        mol = bridge.AtomGroup()

        atom1 = bridge.Atom()
        atom1.symbol = "H"
        atom1.label = "1"
        atom1.xyz.x = 0.0
        atom1.xyz.y = 0.0
        atom1.xyz.z = 0.0
        mol.set_atom("1", atom1)

        atom2 = bridge.Atom()
        atom2.symbol = "H"
        atom2.label = "2"
        atom2.xyz.x = 0.0
        atom2.xyz.y = 0.0
        atom2.xyz.z = 0.0
        mol.set_atom("2", atom2)

        param.molecule = mol

        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [os.path.join(tmpdir, "main.py")]

            # atom1 には [1.0, 0.0, 0.0], atom2 には [0.0, 2.0, 0.0] の勾配を返す
            sim.numerical_grad = MagicMock(
                return_value=([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]], 1.0)
            )

            # 期待値: alpha=0.5、勾配の逆向きのため atom1 は x-0.5, atom2 は y-1.0 に移動すること
            sim.opt(param, max_cycle=1)

            self.assertAlmostEqual(mol["1"].xyz.x, -0.5)
            self.assertAlmostEqual(mol["1"].xyz.y, 0.0)
            self.assertAlmostEqual(mol["2"].xyz.x, 0.0)
            self.assertAlmostEqual(mol["2"].xyz.y, -1.0)

    def test_opt_energy_decreases_with_mock_sp(self):
        """sp をモックした解析的エネルギーで opt を実行し、1ステップ後にエネルギーが下がること、
        pdfparam.molecule が更新され複数原子で各原子が自身の勾配で動くことを確認する。"""
        sim = PdfSim()
        param = PdfParam()
        mol = bridge.AtomGroup()

        atom1 = bridge.Atom()
        atom1.symbol = "H"
        atom1.label = "1"
        atom1.xyz.x = 1.0
        atom1.xyz.y = 0.0
        atom1.xyz.z = 0.0
        mol.set_atom("1", atom1)

        atom2 = bridge.Atom()
        atom2.symbol = "H"
        atom2.label = "2"
        atom2.xyz.x = 0.0
        atom2.xyz.y = 2.0
        atom2.xyz.z = 0.0
        mol.set_atom("2", atom2)

        param.molecule = mol

        def calc_energy(m):
            energy = 0.0
            for atom_id, atom in m.atoms():
                energy += float(atom.xyz.x) ** 2 + float(atom.xyz.y) ** 2 + float(atom.xyz.z) ** 2
            return energy

        initial_energy = calc_energy(param.molecule)
        self.assertAlmostEqual(initial_energy, 5.0)

        def mock_sp(p, workdir):
            return (1, calc_energy(p.molecule))

        sim.sp = mock_sp

        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [os.path.join(tmpdir, "main.py")]
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                sim.opt(param, max_cycle=1)

            final_energy = calc_energy(param.molecule)
            # 1ステップ後にエネルギーが下がること
            self.assertLess(final_energy, initial_energy)

            # 各原子が各自の勾配に応じて更新されていること
            self.assertLess(param.molecule["1"].xyz.x, 0.5)
            self.assertAlmostEqual(param.molecule["1"].xyz.y, 0.0, places=2)
            self.assertAlmostEqual(param.molecule["1"].xyz.z, 0.0, places=2)

            self.assertAlmostEqual(param.molecule["2"].xyz.x, 0.0, places=2)
            self.assertLess(param.molecule["2"].xyz.y, 0.5)
            self.assertAlmostEqual(param.molecule["2"].xyz.z, 0.0, places=2)

            # 解析解に収束 (x=0, y=0)
            self.assertAlmostEqual(param.molecule["1"].xyz.x, 0.0, places=2)
            self.assertAlmostEqual(param.molecule["2"].xyz.y, 0.0, places=2)
            self.assertAlmostEqual(final_energy, 0.0, places=2)

    def test_total_energy_calls_sp(self):
        """total_energy は self.sp を呼び出して (iterations, total_energy) を返すこと。"""
        sim = PdfSim()
        param = PdfParam()
        sim.sp = MagicMock(return_value=(10, -50.0))
        ret = sim.total_energy(param, workdir="dummy_dir")
        sim.sp.assert_called_once_with(param, workdir="dummy_dir")
        self.assertEqual(ret, (10, -50.0))


if __name__ == "__main__":
    unittest.main()
