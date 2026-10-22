#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import unittest
import warnings

import proteindf_bridge as bridge
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO
from proteindf_tools.pdfparam import PdfParam
from proteindf_tools.pdfparam_object import PdfParamObject


class TestPdfParamObject(unittest.TestCase):
    def _create_h2_param(self, label1="", label2=""):
        """Create minimal PdfParam representing an H2 molecule."""
        param = PdfParam()

        mol = bridge.AtomGroup()
        atom1 = bridge.Atom()
        atom1.symbol = "H"
        atom1.label = label1
        atom1.xyz = [0.0, 0.0, 0.0]
        mol.set_atom("1", atom1)

        atom2 = bridge.Atom()
        atom2.symbol = "H"
        atom2.label = label2
        atom2.xyz = [0.0, 0.0, 0.74]
        mol.set_atom("2", atom2)

        param.molecule = mol

        bs = BasisSet()
        bs.name = "sto-3g-H"
        cgto = ContractedGTO("s", 1)
        cgto[0] = PrimitiveGTO(3.42525, 0.15433)
        bs.append(cgto)

        key1 = "H@{}".format(label1) if label1 else "H"
        key2 = "H@{}".format(label2) if label2 else "H"
        param.set_basisset(key1, bs)
        param.set_basisset(key2, bs)

        return param

    # --------------------------------------------------------------------------
    # 1. 初期化・複製・ダイジェスト (__init__, _copy_from_pdfparam, digest)
    # --------------------------------------------------------------------------

    def test_init_default(self):
        """PdfParamObject() initializes with empty _data."""
        param = PdfParamObject()
        self.assertEqual(param._data, {})

    def test_init_with_dict(self):
        """PdfParamObject(dict) sets attributes via set_by_raw_data."""
        data = {
            "comment": "sample comment",
            "method": "rks",
            "control": {},
            "coordinates": {"groups": {}, "atoms": []},
            "guess": "none",
            "num_of_AOs": 0,
            "num_of_MOs": 0,
            "num_of_atoms": 0,
            "num_of_iterations": 0,
        }
        param = PdfParamObject(data)
        self.assertEqual(param.comment, "sample comment")
        self.assertEqual(param.method, "rks")

    def test_init_with_pdfparam_and_copy(self):
        """PdfParamObject(PdfParamObject) and _copy_from_pdfparam create deep copy."""
        orig = PdfParamObject()
        orig.comment = "original comment"
        orig._data["extra_list"] = [1, 2, 3]

        copied = PdfParamObject(orig)
        self.assertEqual(copied.comment, "original comment")
        self.assertEqual(copied._data["extra_list"], [1, 2, 3])

        # Verify deep copy isolation
        orig.comment = "modified comment"
        orig._data["extra_list"].append(4)
        self.assertEqual(copied.comment, "original comment")
        self.assertEqual(copied._data["extra_list"], [1, 2, 3])

        # Test _copy_from_pdfparam directly
        dest = PdfParamObject()
        dest._copy_from_pdfparam(copied)
        self.assertEqual(dest.comment, "original comment")
        copied.comment = "changed again"
        self.assertEqual(dest.comment, "original comment")

    def test_digest(self):
        """digest() returns consistent hash and changes when content changes."""
        p1 = PdfParamObject()
        p2 = PdfParamObject()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            d1 = p1.digest()
            d2 = p2.digest()

        self.assertEqual(d1, d2)
        self.assertIsInstance(d1, str)
        self.assertEqual(len(d1), 32)

        p2.comment = "different content"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            d3 = p2.digest()
        self.assertNotEqual(d1, d3)

    # --------------------------------------------------------------------------
    # 2. 設定値のプロパティ
    # --------------------------------------------------------------------------

    def test_step_control(self):
        """step_control property getter, setter, default, and type casting."""
        param = PdfParamObject()
        self.assertEqual(param.step_control, "")
        self.assertEqual(param._get_step_control(), "")

        param.step_control = "diis"
        self.assertEqual(param.step_control, "diis")

        param._set_step_control(12345)
        self.assertEqual(param.step_control, "12345")

    def test_orbital_independence_threshold(self):
        """orbital_independence_threshold property: getter, setter, float casting, None preservation."""
        param = PdfParamObject()
        self.assertIsNone(param.orbital_independence_threshold)
        self.assertIsNone(param._get_orbital_independence_threshold())

        param.orbital_independence_threshold = 1.0e-5
        self.assertEqual(param.orbital_independence_threshold, 1.0e-5)

        param._set_orbital_independence_threshold("1.0e-4")
        self.assertEqual(param.orbital_independence_threshold, 1.0e-4)

        # None does not overwrite existing value
        param.orbital_independence_threshold = None
        self.assertEqual(param.orbital_independence_threshold, 1.0e-4)

    def test_orbital_independence_threshold_canonical(self):
        """orbital_independence_threshold_canonical: getter, setter, float casting, None preservation."""
        param = PdfParamObject()
        self.assertIsNone(param.orbital_independence_threshold_canonical)
        self.assertIsNone(param._get_orbital_independence_threshold_canonical())

        param.orbital_independence_threshold_canonical = 1.0e-6
        self.assertEqual(param.orbital_independence_threshold_canonical, 1.0e-6)

        param._set_orbital_independence_threshold_canonical("1.0e-3")
        self.assertEqual(param.orbital_independence_threshold_canonical, 1.0e-3)

        param.orbital_independence_threshold_canonical = None
        self.assertEqual(param.orbital_independence_threshold_canonical, 1.0e-3)

    def test_orbital_independence_threshold_lowdin(self):
        """orbital_independence_threshold_lowdin: getter, setter, float casting, None preservation."""
        param = PdfParamObject()
        self.assertIsNone(param.orbital_independence_threshold_lowdin)
        self.assertIsNone(param._get_orbital_independence_threshold_lowdin())

        param.orbital_independence_threshold_lowdin = 1.0e-7
        self.assertEqual(param.orbital_independence_threshold_lowdin, 1.0e-7)

        param._set_orbital_independence_threshold_lowdin("1.0e-2")
        self.assertEqual(param.orbital_independence_threshold_lowdin, 1.0e-2)

        param.orbital_independence_threshold_lowdin = None
        self.assertEqual(param.orbital_independence_threshold_lowdin, 1.0e-2)

    def test_gridfree_dedicated_basis(self):
        """gridfree_dedicated_basis property: getter, setter, boolean conversion."""
        param = PdfParamObject()
        self.assertFalse(param.gridfree_dedicated_basis)
        self.assertFalse(param._get_gridfree_dedicated_basis())

        for val in ["YES", "TRUE", "1", "yes", "true", True]:
            param.gridfree_dedicated_basis = val
            self.assertTrue(param.gridfree_dedicated_basis)

        for val in ["NO", "FALSE", "0", "no", "false", False]:
            param._set_gridfree_dedicated_basis(val)
            self.assertFalse(param.gridfree_dedicated_basis)

    def test_gridfree_orthogonalize_method(self):
        """gridfree_orthogonalize_method property: default 'canonical', string setter, None preservation."""
        param = PdfParamObject()
        self.assertEqual(param.gridfree_orthogonalize_method, "canonical")
        self.assertEqual(param._get_gridfree_orthogonalize_method(), "canonical")

        param.gridfree_orthogonalize_method = "lowdin"
        self.assertEqual(param.gridfree_orthogonalize_method, "lowdin")

        param._set_gridfree_orthogonalize_method(None)
        self.assertEqual(param.gridfree_orthogonalize_method, "lowdin")

    def test_gridfree_CDAM_tau(self):
        """gridfree_CDAM_tau property: default 1.0e-10, float setter, None preservation."""
        param = PdfParamObject()
        self.assertEqual(param.gridfree_CDAM_tau, 1.0e-10)
        self.assertEqual(param._get_gridfree_CDAM_tau(), 1.0e-10)

        param.gridfree_CDAM_tau = 1.0e-8
        self.assertEqual(param.gridfree_CDAM_tau, 1.0e-8)

        param._set_gridfree_CDAM_tau("1.0e-6")
        self.assertEqual(param.gridfree_CDAM_tau, 1.0e-6)

        param.gridfree_CDAM_tau = None
        self.assertEqual(param.gridfree_CDAM_tau, 1.0e-6)

    def test_gridfree_CD_epsilon(self):
        """gridfree_CD_epsilon property: default 1.0e-4, float setter, None preservation."""
        param = PdfParamObject()
        self.assertEqual(param.gridfree_CD_epsilon, 1.0e-4)
        self.assertEqual(param._get_gridfree_CD_epsilon(), 1.0e-4)

        param.gridfree_CD_epsilon = 1.0e-3
        self.assertEqual(param.gridfree_CD_epsilon, 1.0e-3)

        param._set_gridfree_CD_epsilon("1.0e-2")
        self.assertEqual(param.gridfree_CD_epsilon, 1.0e-2)

        param.gridfree_CD_epsilon = None
        self.assertEqual(param.gridfree_CD_epsilon, 1.0e-2)

    def test_extra_keywords(self):
        """extra_keywords property: default empty dict, update on dict, delete on None."""
        param = PdfParamObject()
        self.assertEqual(param.extra_keywords, {})
        self.assertEqual(param._get_extra_keywords(), {})

        param.extra_keywords = {"key1": "val1"}
        self.assertEqual(param.extra_keywords, {"key1": "val1"})

        # Update with another dict
        param._set_extra_keywords({"key2": 42})
        self.assertEqual(param.extra_keywords, {"key1": "val1", "key2": 42})

        # Non-dict is ignored
        param.extra_keywords = "not-a-dict"
        self.assertEqual(param.extra_keywords, {"key1": "val1", "key2": 42})

        # None removes extra_keywords
        param.extra_keywords = None
        self.assertEqual(param.extra_keywords, {})

    def test_occupation_level_auto(self):
        """occupation_level automatically computed from num_of_electrons and cached."""
        param = self._create_h2_param()
        # H2 has 2 electrons -> occupation_level should be "[1 - 1]"
        self.assertEqual(param.num_of_electrons, 2)
        self.assertNotIn("occupation_level", param._data)

        level = param.occupation_level
        self.assertEqual(level, "[1 - 1]")
        self.assertEqual(param._data.get("occupation_level"), "[1 - 1]")

    def test_occupation_level_manual(self):
        """occupation_level manual setting and retrieval."""
        param = PdfParamObject()
        param.occupation_level = "[1 - 10]"
        self.assertEqual(param.occupation_level, "[1 - 10]")
        self.assertEqual(param._get_occupation_level(), "[1 - 10]")

        param._set_occupation_level("[1 - 5]")
        self.assertEqual(param.occupation_level, "[1 - 5]")

    # --------------------------------------------------------------------------
    # 3. 分子・基底関数
    # --------------------------------------------------------------------------

    def test_set_moleucule_type_check(self):
        """_set_moleucule accepts AtomGroup and rejects other types with AssertionError."""
        param = PdfParamObject()
        mol = bridge.AtomGroup()
        param.molecule = mol
        self.assertIs(param.molecule, mol)

        with self.assertRaises(AssertionError):
            param.molecule = "invalid_molecule"

        with self.assertRaises(AssertionError):
            param._set_moleucule(12345)

    def test_basis_set_atomlabels_common(self):
        """_get_basis_set_atomlabels_common returns keys for valid groups, asserts on invalid group."""
        param = PdfParamObject()
        with self.assertRaises(AssertionError):
            param._get_basis_set_atomlabels_common("invalid_group")

        bs = BasisSet()
        param.set_basisset("H", bs)
        param.set_basisset_j("H", bs)
        param.set_basisset_xc("C", bs)
        param.set_basisset_gridfree("N", bs)

        self.assertIn("H", param._get_basis_set_atomlabels_common("basis_set"))
        self.assertIn("H", param.get_basisset_atomlabels())

        self.assertIn("H", param._get_basis_set_atomlabels_common("basis_set_j"))
        self.assertIn("H", param.get_basisset_j_atomlabels())

        self.assertIn("C", param._get_basis_set_atomlabels_common("basis_set_xc"))
        self.assertIn("C", param.get_basisset_xc_atomlabels())

        self.assertIn("N", param._get_basis_set_atomlabels_common("basis_set_gridfree"))
        self.assertIn("N", param.get_basisset_gridfree_atomlabels())

    def test_basisset_common(self):
        """_get_basisset_common returns BasisSet for atom label or empty BasisSet if not set."""
        param = PdfParamObject()
        empty_bs = param._get_basisset_common("H", "basis_set_xc")
        self.assertIsInstance(empty_bs, BasisSet)
        self.assertEqual(len(empty_bs), 0)

        bs = BasisSet()
        bs.name = "custom-xc"
        param.set_basisset_xc("H", bs)
        retrieved_xc = param.get_basisset_xc("H")
        self.assertEqual(retrieved_xc.name, "custom-xc")

        bs_gf = BasisSet()
        bs_gf.name = "custom-gf"
        param.set_basisset_gridfree("H", bs_gf)
        retrieved_gf = param.get_basisset_gridfree("H")
        self.assertEqual(retrieved_gf.name, "custom-gf")

    def test_pickup_atom_kinds_without_label(self):
        """_pickup_atom_kinds collects atom kinds when atom labels are empty."""
        param = self._create_h2_param(label1="", label2="")
        kinds = param._pickup_atom_kinds()
        self.assertEqual(kinds, {"H"})

        # Add dummy atom 'X'
        dummy = bridge.Atom()
        dummy.symbol = "X"
        dummy.label = ""
        param.molecule.set_atom("3", dummy)

        # When counterpoise is False, 'X' is ignored
        self.assertEqual(param._pickup_atom_kinds(), {"H"})

        # When counterpoise is True, 'X' is included
        param._data["counterpoise"] = True
        self.assertEqual(param._pickup_atom_kinds(), {"H", "X"})

    @unittest.expectedFailure
    def test_pickup_atom_kinds_with_label_expected_failure(self):
        """_pickup_atom_kinds collects atom kinds with labels formatted as symbol@label.

        proteindf_tools/pdfparam_object.py:1009: atom.label が非空のとき
        atom_kind = "%s@%s" % (atom_symbol, label) で未定義変数 label
        (正しくは atom_label) を参照しているため NameError が発生する。
        修正されたら外す。
        """
        param = self._create_h2_param(label1="1", label2="2")
        kinds = param._pickup_atom_kinds()
        self.assertEqual(kinds, {"H@1", "H@2"})

    def test_get_input_basisset_methods(self):
        """_get_input_basisset, _j, _xc, _gridfree return formatted input strings."""
        param = self._create_h2_param(label1="", label2="")
        param.set_basisset_name("H", "sto-3g")
        param.set_basisset_j_name("H", "sto-3g-j")
        param.set_basisset_xc_name("H", "sto-3g-xc")
        param.set_basisset_gridfree_name("H", "sto-3g-gf")

        self.assertEqual(param._get_input_basisset(), 'H = "sto-3g"\n')
        self.assertEqual(param._get_input_basisset_j(), 'H = "sto-3g-j"\n')
        self.assertEqual(param._get_input_basisset_xc(), 'H = "sto-3g-xc"\n')
        self.assertEqual(param._get_input_basisset_gridfree(), 'H = "sto-3g-gf"\n')

    def test_get_input_geometry_with_atom_group(self):
        """_get_input_geometry with explicit atom_group returns formatted geometry string."""
        param = self._create_h2_param(label1="", label2="")
        geom_str = param._get_input_geometry(param.molecule)
        self.assertIn("H  0.000000  0.000000  0.000000", geom_str)
        self.assertIn("H  0.000000  0.000000  0.740000", geom_str)

    @unittest.expectedFailure
    def test_get_input_geometry_default_expected_failure(self):
        """_get_input_geometry() with default atom_group=None returns formatted geometry.

        proteindf_tools/pdfparam_object.py:912: atom_group=None のとき
        self._get_geometry_input(self.molecule) を呼び出すが、
        _get_geometry_input メソッドが未定義 (正しくは _get_input_geometry) のため
        AttributeError が発生する。
        修正されたら外す。
        """
        param = self._create_h2_param(label1="", label2="")
        geom_str = param._get_input_geometry()
        self.assertIn("H", geom_str)

    @unittest.expectedFailure
    def test_get_input_geometry_with_label_expected_failure(self):
        """_get_input_geometry formats atom labels when label is non-empty.

        proteindf_tools/pdfparam_object.py:926: atom.label が非空のとき
        re.sub(r"\\s+", "_", label) を呼び出すが、モジュールで re が import
        されていないため NameError が発生する。
        修正されたら外す。
        """
        param = self._create_h2_param(label1="1", label2="2")
        geom_str = param._get_input_geometry(param.molecule)
        self.assertIn("H@1", geom_str)

    # --------------------------------------------------------------------------
    # 4. 計算結果の取得
    # --------------------------------------------------------------------------

    def test_counterpoise(self):
        """counterpoise property returns default False and reflects set value."""
        param = PdfParamObject()
        self.assertFalse(param.counterpoise)
        self.assertFalse(param._get_counterpoise())

        param._data["counterpoise"] = True
        self.assertTrue(param.counterpoise)
        self.assertTrue(param._get_counterpoise())

    def test_lo_num_of_iterations(self):
        """lo_num_of_iterations returns default None and reflects set value."""
        param = PdfParamObject()
        self.assertIsNone(param.lo_num_of_iterations)
        self.assertIsNone(param._get_lo_num_of_iterations())

        param._data["lo/num_of_iterations"] = 15
        self.assertEqual(param.lo_num_of_iterations, 15)
        self.assertEqual(param._get_lo_num_of_iterations(), 15)

    def test_lo_satisfied_default(self):
        """lo_satisfied returns False when unset."""
        param = PdfParamObject()
        self.assertFalse(param.lo_satisfied)
        self.assertFalse(param._get_lo_satisfied())

    @unittest.expectedFailure
    def test_lo_satisfied_true_expected_failure(self):
        """lo_satisfied returns True when set to True.

        proteindf_tools/pdfparam_object.py:756: isinstance(answer, bool) のとき
        else 節で answer = answer == "YES" (True == "YES") が評価されて
        False になってしまう不具合があるため失敗する。
        修正されたら外す。
        """
        param = PdfParamObject()
        param._data["lo/satisfied"] = True
        self.assertTrue(param.lo_satisfied)

    def test_get_clo_mat_path(self):
        """get_clo_mat_path returns localized matrix file path."""
        param = PdfParamObject()
        param.work_path = "/work/calc"
        param._data.setdefault("control", {})["file_base_name"] = {
            "Clo_matrix": "clo.%s.mat"
        }

        # Explicit iteration
        path1 = param.get_clo_mat_path("rks", 5)
        self.assertEqual(path1, "/work/calc/clo.rks5.mat")

        # Default iteration (-1 uses lo_num_of_iterations)
        param._data["lo/num_of_iterations"] = 8
        path2 = param.get_clo_mat_path("uks_alpha")
        self.assertEqual(path2, "/work/calc/clo.uks_alpha8.mat")

    # --------------------------------------------------------------------------
    # 5. その他 (_is_fitting_xc, _alias_conversion)
    # --------------------------------------------------------------------------

    def test_is_fitting_xc(self):
        """_is_fitting_xc returns True if xc_functional ends with ~, False otherwise."""
        param = PdfParamObject()
        param.xc_functional = "SVWN"
        self.assertFalse(param._is_fitting_xc())

        param.xc_functional = "B3LYP~"
        self.assertTrue(param._is_fitting_xc())

        param.xc_functional = "PBE"
        self.assertFalse(param._is_fitting_xc())

    def test_alias_conversion(self):
        """_alias_conversion translates alias keys and values to canonical forms."""
        param = PdfParamObject()
        aliases = {
            "method": "nsp",
            "scf-start-guess": "core",
            "xc-potential": "SVWN",
            "max-iteration": 100,
            "custom_key": "custom_val",
        }
        converted = param._alias_conversion(aliases)
        self.assertEqual(converted["method"], "rks")
        self.assertEqual(converted["guess"], "core")
        self.assertEqual(converted["xc_functional"], "SVWN")
        self.assertEqual(converted["max_iteration"], 100)
        self.assertEqual(converted["custom_key"], "custom_val")

        # Test "sp" -> "uks"
        converted_sp = param._alias_conversion({"method": "sp"})
        self.assertEqual(converted_sp["method"], "uks")


if __name__ == "__main__":
    unittest.main()
