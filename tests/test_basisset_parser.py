#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from proteindf_tools.basisset_parser import BasisSetParser


class TestBasisSetParser(unittest.TestCase):
    def setUp(self):
        self.parser = BasisSetParser()

    def test_parse_s_and_p_shells(self):
        # Gaussian94 format with single-character shells (s, p)
        sample = """
H 0
S 2 1.00
      33.8700000D+00           0.0254938
       5.0950000D+00           0.1903730
P 1 1.00
       1.1590000E+00           1.0000000
****
"""
        basissets, basissets_density, basissets_xc = self.parser.parse(sample)
        self.assertIn("H", basissets)
        self.assertEqual(len(basissets["H"]), 2)

        cgto_s = basissets["H"][0]
        self.assertEqual(cgto_s.shell_type, "s")
        self.assertEqual(len(cgto_s), 2)
        self.assertAlmostEqual(cgto_s[0].exp, 33.87)
        self.assertAlmostEqual(cgto_s[0].coef, 0.0254938)
        self.assertAlmostEqual(cgto_s[1].exp, 5.095)
        self.assertAlmostEqual(cgto_s[1].coef, 0.1903730)

        cgto_p = basissets["H"][1]
        self.assertEqual(cgto_p.shell_type, "p")
        self.assertEqual(len(cgto_p), 1)
        self.assertAlmostEqual(cgto_p[0].exp, 1.159)
        self.assertAlmostEqual(cgto_p[0].coef, 1.0)

    def test_parse_sp_shell(self):
        # SP shell should be expanded into separate s and p CGTOs
        sample = """
C 0
SP 2 1.00
       2.9412494D+00     -0.09996723D+00     0.15591627D+00
       0.6834831D+00      0.39951283D+00     0.60768372D+00
****
"""
        basissets, _, _ = self.parser.parse(sample)
        self.assertIn("C", basissets)
        self.assertEqual(len(basissets["C"]), 2)

        cgto_s = basissets["C"][0]
        cgto_p = basissets["C"][1]

        self.assertEqual(cgto_s.shell_type, "s")
        self.assertEqual(cgto_p.shell_type, "p")

        self.assertEqual(len(cgto_s), 2)
        self.assertEqual(len(cgto_p), 2)

        # exponents are shared
        self.assertAlmostEqual(cgto_s[0].exp, 2.9412494)
        self.assertAlmostEqual(cgto_p[0].exp, 2.9412494)
        self.assertAlmostEqual(cgto_s[1].exp, 0.6834831)
        self.assertAlmostEqual(cgto_p[1].exp, 0.6834831)

        # coefficients differ
        self.assertAlmostEqual(cgto_s[0].coef, -0.09996723)
        self.assertAlmostEqual(cgto_p[0].coef, 0.15591627)
        self.assertAlmostEqual(cgto_s[1].coef, 0.39951283)
        self.assertAlmostEqual(cgto_p[1].coef, 0.60768372)

    def test_parse_multiple_blocks(self):
        sample = """
H 0
S 1 1.00
  1.0 1.0
****
C 0
S 1 1.00
  2.0 1.0
****
H 0
S 1 1.00
  3.0 1.0
****
H 0
S 1 1.00
  4.0 1.0
****
"""
        basissets, basissets_density, basissets_xc = self.parser.parse(sample)
        self.assertIn("H", basissets)
        self.assertIn("C", basissets)
        self.assertIn("H", basissets_density)
        self.assertIn("H", basissets_xc)

        self.assertAlmostEqual(basissets["H"][0][0].exp, 1.0)
        self.assertAlmostEqual(basissets["C"][0][0].exp, 2.0)
        self.assertAlmostEqual(basissets_density["H"][0][0].exp, 3.0)
        self.assertAlmostEqual(basissets_xc["H"][0][0].exp, 4.0)

    def test_invalid_atom_line_ignored(self):
        sample = """
INVALID_LINE_WITHOUT_ATOM_0
S 1 1.00
  1.0 1.0
****
"""
        basissets, basissets_density, basissets_xc = self.parser.parse(sample)
        self.assertEqual(basissets, {})
        self.assertEqual(basissets_density, {})
        self.assertEqual(basissets_xc, {})

    def test_unsupported_shell_type_raises(self):
        sample = """
C 0
SPD 1 1.00
  1.0 1.0
****
"""
        with self.assertRaises(RuntimeError):
            self.parser.parse(sample)

    def test_too_many_basis_blocks_same_atom(self):
        # 4th block for the same atom attempts sys.stderr.write, but sys is not imported
        # in basisset_parser.py, raising NameError: name 'sys' is not defined.
        sample = """
H 0
S 1 1.00
  1.0 1.0
****
H 0
S 1 1.00
  1.0 1.0
****
H 0
S 1 1.00
  1.0 1.0
****
H 0
S 1 1.00
  1.0 1.0
****
"""
        # Expected behavior: parses without unhandled NameError (warns or ignores)
        self.parser.parse(sample)


if __name__ == "__main__":
    unittest.main()
