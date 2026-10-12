#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO


class TestPrimitiveGTO(unittest.TestCase):
    def test_raw_data_roundtrip(self):
        pgto = PrimitiveGTO(1.234, 0.567)
        raw = pgto.get_raw_data()
        self.assertEqual(raw, {"exp": 1.234, "coef": 0.567})

        pgto_restored = PrimitiveGTO()
        pgto_restored.set_by_raw_data(raw)
        self.assertEqual(pgto, pgto_restored)


class TestContractedGTO(unittest.TestCase):
    def test_shell_type_names_and_ids(self):
        expected_types = ["s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m"]
        self.assertEqual(ContractedGTO.get_supported_shell_types(), expected_types)

        for expected_id, st in enumerate(expected_types):
            self.assertEqual(ContractedGTO.get_shell_type_id(st), expected_id)
            self.assertEqual(ContractedGTO.get_shell_type(expected_id), st)

        with self.assertRaises(ValueError):
            ContractedGTO.get_shell_type_id("invalid")

    def test_shell_type_properties(self):
        cgto = ContractedGTO("p", 1)
        self.assertEqual(cgto.shell_type, "p")
        self.assertEqual(cgto.shell_type_id, 1)

        cgto.shell_type = "d"
        self.assertEqual(cgto.shell_type, "d")
        self.assertEqual(cgto.shell_type_id, 2)

        cgto.shell_type_id = 0
        self.assertEqual(cgto.shell_type, "s")
        self.assertEqual(cgto.shell_type_id, 0)

    def test_get_basis_type(self):
        # s
        self.assertEqual(ContractedGTO.get_basis_type(0, 0), "s")
        # p
        self.assertEqual(ContractedGTO.get_basis_type(1, 0), "px")
        self.assertEqual(ContractedGTO.get_basis_type(1, 1), "py")
        self.assertEqual(ContractedGTO.get_basis_type(1, 2), "pz")
        # d
        self.assertEqual(ContractedGTO.get_basis_type(2, 0), "dxy")
        self.assertEqual(ContractedGTO.get_basis_type(2, 1), "dyz")
        self.assertEqual(ContractedGTO.get_basis_type(2, 2), "dzx")
        self.assertEqual(ContractedGTO.get_basis_type(2, 3), "dxx-yy")
        self.assertEqual(ContractedGTO.get_basis_type(2, 4), "dzz")

    def test_get_basis_type_f_orbital(self):
        self.assertEqual(ContractedGTO.get_basis_type(3, 0), "z3")
        with self.assertRaises(ValueError):
            ContractedGTO.get_basis_type(4, 0)

    def test_normalize(self):
        # Normalized single primitive s GTO
        cgto_s = ContractedGTO("s", 1)
        cgto_s[0] = PrimitiveGTO(1.0, 1.0)
        self.assertAlmostEqual(cgto_s.normalize(), 1.0, places=7)

        # Doubled coefficient should yield factor 0.5
        cgto_s2 = ContractedGTO("s", 1)
        cgto_s2[0] = PrimitiveGTO(1.0, 2.0)
        self.assertAlmostEqual(cgto_s2.normalize(), 0.5, places=7)

        # p and d normalization
        cgto_p = ContractedGTO("p", 1)
        cgto_p[0] = PrimitiveGTO(1.0, 1.0)
        self.assertAlmostEqual(cgto_p.normalize(), 1.0, places=7)

        cgto_d = ContractedGTO("d", 1)
        cgto_d[0] = PrimitiveGTO(1.0, 1.0)
        self.assertAlmostEqual(cgto_d.normalize(), 1.0, places=7)

    def test_expand_standard_shell(self):
        cgto = ContractedGTO("s", 1)
        cgto[0] = PrimitiveGTO(1.0, 1.0)
        expanded = cgto.expand()
        self.assertEqual(len(expanded), 1)
        self.assertIs(expanded[0], cgto)

    def test_raw_data_roundtrip(self):
        cgto = ContractedGTO("p", 2)
        cgto.scale_factor = 1.5
        cgto[0] = PrimitiveGTO(2.0, 0.4)
        cgto[1] = PrimitiveGTO(0.5, 0.6)

        raw = cgto.get_raw_data()
        self.assertEqual(raw["shell_type"], "p")
        self.assertEqual(raw["scale_factor"], 1.5)
        self.assertEqual(len(raw["pGTOs"]), 2)

        restored = ContractedGTO()
        restored.set_by_raw_data(raw)
        self.assertEqual(cgto, restored)


class TestBasisSet(unittest.TestCase):
    def test_get_number_of_aos(self):
        bs = BasisSet("test")
        cgto_s = ContractedGTO("s", 1)
        cgto_s[0] = PrimitiveGTO(1.0, 1.0)
        cgto_p = ContractedGTO("p", 1)
        cgto_p[0] = PrimitiveGTO(1.0, 1.0)
        cgto_d = ContractedGTO("d", 1)
        cgto_d[0] = PrimitiveGTO(1.0, 1.0)

        bs.append(cgto_s)
        self.assertEqual(bs.get_number_of_AOs(), 1)  # 0*2+1 = 1

        bs.append(cgto_p)
        self.assertEqual(bs.get_number_of_AOs(), 4)  # 1 + 3 = 4

        bs.append(cgto_d)
        self.assertEqual(bs.get_number_of_AOs(), 9)  # 4 + 5 = 9

    def test_get_num_of_cgtos(self):
        bs = BasisSet("test")
        bs.append(ContractedGTO("s", 1))
        bs.append(ContractedGTO("s", 1))
        bs.append(ContractedGTO("p", 1))

        self.assertEqual(bs.get_num_of_CGTOs("s"), 2)
        self.assertEqual(bs.get_num_of_CGTOs("p"), 1)
        self.assertEqual(bs.get_num_of_CGTOs("d"), 0)

    def test_sort(self):
        bs = BasisSet("test")
        cgto_d = ContractedGTO("d", 1)
        cgto_s = ContractedGTO("s", 1)
        cgto_p = ContractedGTO("p", 1)

        bs.append(cgto_d)
        bs.append(cgto_s)
        bs.append(cgto_p)

        self.assertEqual([c.shell_type for c in bs], ["d", "s", "p"])
        bs.sort()
        self.assertEqual([c.shell_type for c in bs], ["s", "p", "d"])

    def test_expand(self):
        bs = BasisSet("test_bs")
        cgto_s = ContractedGTO("s", 1)
        cgto_s[0] = PrimitiveGTO(1.0, 1.0)
        cgto_p = ContractedGTO("p", 1)
        cgto_p[0] = PrimitiveGTO(2.0, 1.0)
        bs.append(cgto_s)
        bs.append(cgto_p)

        expanded = bs.expand()
        self.assertEqual(expanded.name, "test_bs")
        self.assertEqual(len(expanded), 2)
        self.assertEqual([c.shell_type for c in expanded], ["s", "p"])
        # Original BasisSet should be emptied by expand()
        self.assertEqual(len(bs), 0)

    def test_raw_data_roundtrip(self):
        bs = BasisSet("water_basis")
        cgto_s = ContractedGTO("s", 1)
        cgto_s[0] = PrimitiveGTO(10.0, 0.2)
        cgto_p = ContractedGTO("p", 1)
        cgto_p[0] = PrimitiveGTO(1.0, 0.8)
        bs.append(cgto_s)
        bs.append(cgto_p)

        raw = bs.get_raw_data()
        self.assertEqual(raw["name"], "water_basis")
        self.assertEqual(len(raw["cGTOs"]), 2)

        restored = BasisSet()
        restored.set_by_raw_data(raw)
        self.assertEqual(bs, restored)


if __name__ == "__main__":
    unittest.main()
