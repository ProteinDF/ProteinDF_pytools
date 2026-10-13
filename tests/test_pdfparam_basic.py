#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from unittest.mock import patch

import proteindf_bridge as bridge
from proteindf_tools.pdfparam import PdfParam


class TestPdfParamBasic(unittest.TestCase):
    # ------------------------------------------------------------------
    # 1. 既定値とsetterの型変換
    # ------------------------------------------------------------------
    def test_default_values(self):
        param = PdfParam()
        self.assertEqual(param.work_path, "./fl_Work")
        self.assertEqual(param.method, "rks")
        self.assertAlmostEqual(param.cut_value, 1.0e-10)
        self.assertAlmostEqual(param.CDAM_tau, 1.0e-10)
        self.assertAlmostEqual(param.CD_epsilon, 1.0e-4)
        self.assertEqual(param.comment, "")

    def test_setter_type_conversion_and_none(self):
        param = PdfParam()

        # work_path
        param.work_path = 123
        self.assertEqual(param.work_path, "123")

        # comment: _data["comment"] が None のときは "" を返す
        param.comment = "initial comment"
        self.assertEqual(param.comment, "initial comment")
        param._data["comment"] = None
        self.assertEqual(param.comment, "")

        # cut_value: float conversion, None is ignored
        param.cut_value = "1.0e-5"
        self.assertAlmostEqual(param.cut_value, 1.0e-5)
        param.cut_value = None
        self.assertAlmostEqual(param.cut_value, 1.0e-5)

        # CDAM_tau: float conversion, None is ignored
        param.CDAM_tau = "2.0e-8"
        self.assertAlmostEqual(param.CDAM_tau, 2.0e-8)
        param.CDAM_tau = None
        self.assertAlmostEqual(param.CDAM_tau, 2.0e-8)

        # CD_epsilon: float conversion, None is ignored
        param.CD_epsilon = "3.0e-3"
        self.assertAlmostEqual(param.CD_epsilon, 3.0e-3)
        param.CD_epsilon = None
        self.assertAlmostEqual(param.CD_epsilon, 3.0e-3)

    def test_method_validation_and_case_conversion(self):
        param = PdfParam()

        # 大文字は小文字に変換
        param.method = "RKS"
        self.assertEqual(param.method, "rks")
        param.method = "UKS"
        self.assertEqual(param.method, "uks")
        param.method = "ROKS"
        self.assertEqual(param.method, "roks")

        # rks/uks/roks 以外で AssertionError
        with self.assertRaises(AssertionError):
            param.method = "hf"
        with self.assertRaises(AssertionError):
            param.method = "dft"

    # ------------------------------------------------------------------
    # 2. run_types, SCF/加速/基底/エンジン等のプロパティ取得・設定
    # ------------------------------------------------------------------
    def test_run_types(self):
        param = PdfParam()
        param.method = "rks"
        self.assertEqual(param.run_types(), ["rks"])

        param.method = "uks"
        self.assertEqual(param.run_types(), ["uks_alpha", "uks_beta"])

        param.method = "roks"
        self.assertEqual(param.run_types(), ["roks_close", "roks_open"])

    def test_iterations_and_max_iterations(self):
        param = PdfParam()
        # iterations既定値はNone
        self.assertIsNone(param.iterations)
        param.iterations = "15"
        self.assertEqual(param.iterations, 15)
        param.iterations = None
        self.assertEqual(param.iterations, 0)

        # max_iterations既定値は100
        self.assertEqual(param.max_iterations, 100)
        param.max_iterations = "200"
        self.assertEqual(param.max_iterations, 200)

    def test_convergence_properties(self):
        param = PdfParam()
        self.assertAlmostEqual(param.convergence_threshold, 1.0e-3)
        self.assertEqual(param.convergence_type, "density")
        self.assertAlmostEqual(param.convergence_threshold_energy, 1.0e-4)

        param.convergence_threshold = "1.0e-5"
        self.assertAlmostEqual(param.convergence_threshold, 1.0e-5)
        param.convergence_threshold = None
        self.assertAlmostEqual(param.convergence_threshold, 1.0e-5)

        param.convergence_type = "energy"
        self.assertEqual(param.convergence_type, "energy")
        param.convergence_type = None
        self.assertEqual(param.convergence_type, "energy")

        param.convergence_threshold_energy = "1.0e-6"
        self.assertAlmostEqual(param.convergence_threshold_energy, 1.0e-6)
        param.convergence_threshold_energy = None
        self.assertAlmostEqual(param.convergence_threshold_energy, 1.0e-6)

    def test_xc_functional_and_engines(self):
        param = PdfParam()
        # 既定値
        self.assertEqual(param.xc_functional, "SVWN")
        self.assertEqual(param.j_engine, "RI_J")
        self.assertEqual(param.k_engine, "conventional")
        self.assertEqual(param.xc_engine, "grid")

        # 設定
        param.xc_functional = "B3LYP"
        self.assertEqual(param.xc_functional, "B3LYP")
        param.xc_functional = None
        self.assertEqual(param.xc_functional, "SVWN")

        param.j_engine = "conventional"
        self.assertEqual(param.j_engine, "conventional")
        param.j_engine = None
        self.assertEqual(param.j_engine, "RI_J")

        param.k_engine = "RI_K"
        self.assertEqual(param.k_engine, "RI_K")
        param.k_engine = None
        self.assertEqual(param.k_engine, "conventional")

        param.xc_engine = "gridfree"
        self.assertEqual(param.xc_engine, "gridfree")
        param.xc_engine = None
        self.assertEqual(param.xc_engine, "grid")

    def test_scf_acceleration_properties(self):
        param = PdfParam()
        self.assertEqual(param.scf_acceleration, "damping")
        self.assertAlmostEqual(param.scf_acceleration_damping_damping_factor, 0.85)
        self.assertEqual(param.scf_acceleration_damping_damping_type, "density_matrix")
        self.assertIsNone(param.scf_acceleration_anderson_start)

        param.scf_acceleration = "anderson"
        self.assertEqual(param.scf_acceleration, "anderson")
        param.scf_acceleration_damping_damping_factor = "0.5"
        self.assertAlmostEqual(param.scf_acceleration_damping_damping_factor, 0.5)
        param.scf_acceleration_damping_damping_type = "density"
        self.assertEqual(param.scf_acceleration_damping_damping_type, "density")
        param.scf_acceleration_anderson_start = "5"
        self.assertEqual(param.scf_acceleration_anderson_start, 5)

    def test_level_shift_virtual_mo(self):
        param = PdfParam()
        self.assertAlmostEqual(param.level_shift_virtual_mo, 0.0)
        param.level_shift_virtual_mo = "0.2"
        self.assertAlmostEqual(param.level_shift_virtual_mo, 0.2)

    # プロパティ・setter
    def test_scf_converged_setter(self):
        """scf_converged に値を設定でき、読み戻せる。"""
        param = PdfParam()
        self.assertFalse(param.scf_converged)
        param.scf_converged = True
        self.assertTrue(param.scf_converged)

    def test_level_shift_setter(self):
        """level_shift に文字列を設定でき、真偽値として読み戻せる。"""
        param = PdfParam()
        self.assertFalse(param.level_shift)
        param.level_shift = "yes"
        self.assertTrue(param.level_shift)

    def test_level_shift_start_iteration_setter(self):
        """level_shift_start_iteration に値を設定でき、読み戻せる。"""
        param = PdfParam()
        self.assertEqual(param.level_shift_start_iteration, 1)
        param.level_shift_start_iteration = 5
        self.assertEqual(param.level_shift_start_iteration, 5)

    def test_scf_acceleration_anderson_damping_factor_default(self):
        """scf_acceleration_anderson_damping_factor の初期値として None が取得できる。"""
        param = PdfParam()
        self.assertIsNone(param.scf_acceleration_anderson_damping_factor)

    # ------------------------------------------------------------------
    # 3. get_raw_data ↔ set_by_raw_data、コピー独立性、__eq__
    # ------------------------------------------------------------------
    def test_copy_independence_from_pdfparam(self):
        p1 = PdfParam()
        p1.method = "uks"
        p1.iterations = 10
        p1.work_path = "/tmp/calc1"
        p1.cut_value = 1.0e-8

        p2 = PdfParam(p1)
        self.assertEqual(p2.method, "uks")
        self.assertEqual(p2.iterations, 10)
        self.assertEqual(p2.work_path, "/tmp/calc1")

        p2.method = "roks"
        p2.iterations = 20
        p2.work_path = "/tmp/calc2"

        self.assertEqual(p1.method, "uks")
        self.assertEqual(p1.iterations, 10)
        self.assertEqual(p1.work_path, "/tmp/calc1")
        self.assertEqual(p2.method, "roks")
        self.assertEqual(p2.iterations, 20)
        self.assertEqual(p2.work_path, "/tmp/calc2")

    def test_copy_independence_from_raw_dict(self):
        raw = {
            "method": "roks",
            "work_path": "/tmp/raw_test",
            "control": {},
            "coordinates": {"groups": {}, "atoms": []},
            "guess": "none",
            "num_of_AOs": 12,
            "num_of_MOs": 6,
            "num_of_atoms": 2,
            "num_of_iterations": 7,
        }
        p = PdfParam(raw)
        self.assertEqual(p.method, "roks")
        self.assertEqual(p.work_path, "/tmp/raw_test")
        self.assertEqual(p.iterations, 7)

        raw["method"] = "rks"
        raw["work_path"] = "/tmp/changed"
        raw["num_of_iterations"] = 99
        self.assertEqual(p.method, "roks")
        self.assertEqual(p.work_path, "/tmp/raw_test")
        self.assertEqual(p.iterations, 7)

    def test_raw_data_roundtrip_with_required_fields(self):
        p1 = PdfParam()
        p1.molecule = bridge.AtomGroup()
        p1.work_path = "/tmp/roundtrip"
        p1.method = "uks"
        p1.comment = "roundtrip test"
        p1.cut_value = 1.0e-8
        _ = p1.scf_converged
        p1.num_of_AOs = 10
        p1.num_of_MOs = 5
        p1.num_of_atoms = 2
        p1.iterations = 3

        raw = p1.get_raw_data()
        raw["guess"] = "core"

        p2 = PdfParam(raw)
        self.assertEqual(p2.work_path, "/tmp/roundtrip")
        self.assertEqual(p2.method, "uks")
        self.assertEqual(p2.comment, "roundtrip test")
        self.assertAlmostEqual(p2.cut_value, 1.0e-8)
        self.assertEqual(p2.num_of_AOs, 10)
        self.assertEqual(p2.num_of_MOs, 5)
        self.assertEqual(p2.num_of_atoms, 2)
        self.assertEqual(p2.iterations, 3)

    def test_equality_operator(self):
        p1 = PdfParam()
        p1.molecule = bridge.AtomGroup()
        p1.num_of_AOs = 10
        p1.num_of_MOs = 5
        p1.num_of_atoms = 2
        p1.iterations = 0

        p2 = PdfParam(p1)
        self.assertEqual(p1, p2)

        # 異なる型との比較
        self.assertNotEqual(p1, "not_a_pdfparam")

        # 属性が異なると非等価
        p3 = PdfParam(p1)
        p3.num_of_AOs = 20
        self.assertNotEqual(p1, p3)

    def test_get_raw_data_uninitialized_molecule(self):
        """molecule 未設定の空 PdfParam から get_raw_data() で辞書を取得できる。"""
        param = PdfParam()
        raw = param.get_raw_data()
        self.assertIsInstance(raw, dict)

    def test_equality_default_objects(self):
        """デフォルト状態の PdfParam 同士を等値比較できる。"""
        p1 = PdfParam()
        p2 = PdfParam()
        self.assertEqual(p1, p2)

    # ------------------------------------------------------------------
    # 4. パス系メソッド
    # ------------------------------------------------------------------
    def test_path_methods_with_file_base_names(self):
        param = PdfParam()
        param.work_path = "/work/calc"
        param._data.setdefault("control", {})["file_base_name"] = {
            "eigenvalues": "eigenvalues.%s.vtr",
            "occupation_vtr": "occ.%s.vtr",
            "Spq_matrix": "s.mat",
            "X_matrix": "x.mat",
            "Xinv_matrix": "xinv.mat",
            "Hpq_matrix": "h.mat",
            "Hpq2_matrix": "h2.mat",
            "Fpq_matrix": "f.%s.mat",
            "C_matrix": "c.%s.mat",
            "Ppq_matrix": "p.%s.mat",
            "Clo_matrix": "clo.%s.mat",
            "gradient": "grad.mat",
            "pop/mulliken": "pop.mulliken.%s.vtr",
        }

        # energy level
        self.assertEqual(
            param.get_energy_level_path("rks", 5),
            "/work/calc/eigenvalues.rks5.vtr",
        )
        self.assertEqual(
            param.get_energy_level_path("uks_alpha", 3),
            "/work/calc/eigenvalues.uks_alpha3.vtr",
        )

        # occupation
        self.assertEqual(param.get_occ_path("rks"), "/work/calc/occ.rks.vtr")
        self.assertEqual(param.get_occ_path("uks_alpha"), "/work/calc/occ.uks_alpha.vtr")

        # overlap and related matrices
        self.assertEqual(param.get_s_mat_path(), "/work/calc/s.mat")
        self.assertEqual(param.get_x_mat_path(), "/work/calc/x.mat")
        self.assertEqual(param.get_xinv_mat_path(), "/work/calc/xinv.mat")
        self.assertEqual(param.get_h_mat_path(), "/work/calc/h.mat")
        self.assertEqual(param.get_h2_mat_path(), "/work/calc/h2.mat")
        self.assertEqual(param.get_gradient_matrix_path(), "/work/calc/grad.mat")

        # F, C, Density matrix (iteration specified)
        self.assertEqual(param.get_f_mat_path("rks", 5), "/work/calc/f.rks5.mat")
        self.assertEqual(param.get_c_mat_path("rks", 5), "/work/calc/c.rks5.mat")
        self.assertEqual(param.get_density_matrix_path("rks", 5), "/work/calc/p.rks5.mat")

        # F, C, Density matrix (default itr=-1 uses param.iterations)
        param.iterations = 10
        self.assertEqual(param.get_f_mat_path("rks"), "/work/calc/f.rks10.mat")
        self.assertEqual(param.get_c_mat_path("rks"), "/work/calc/c.rks10.mat")
        self.assertEqual(param.get_density_matrix_path("rks"), "/work/calc/p.rks10.mat")

        # pop mulliken
        self.assertEqual(
            param.get_pop_mulliken_path("rks", 5),
            "/work/calc/pop.mulliken.rks5.vtr",
        )
        self.assertEqual(
            param.get_pop_mulliken_path("rks"),
            "/work/calc/pop.mulliken.rks10.vtr",
        )

    def test_pop_mulliken_path_default_fallback(self):
        param = PdfParam()
        param.work_path = "/work/default"
        # file_base_name["pop/mulliken"] が未設定の場合は pop.mulliken.%s.vtr にフォールバック
        self.assertEqual(
            param.get_pop_mulliken_path("rks", 1),
            "/work/default/pop.mulliken.rks1.vtr",
        )

    @patch("proteindf_tools.pdfparam_object.Basis2")
    def test_set_basisset_common_unknown_key(self, mock_basis2):
        param = PdfParam()
        with self.assertRaises(ValueError) as ctx:
            param._set_basisset_common("H", "O-DZVP.H", "unknown_key")
        self.assertEqual(str(ctx.exception), "unknown key: unknown_key")


if __name__ == "__main__":
    unittest.main()
