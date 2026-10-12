#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import textwrap
import unittest

from proteindf_tools.basis2 import Basis2
from proteindf_tools.basisset import BasisSet, ContractedGTO, PrimitiveGTO


class TestBasis2(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_crete(self):
        bs2 = Basis2()

    def test_DZVP(self):
        expect = """
        O-DZVP.H
         2 0 0 0 0 0 0 0 0 0 0
         4
             5.099918e+01  9.660500e-03
             7.483218e+00  7.372890e-02
             1.777468e+00  2.958581e-01
             5.193295e-01  7.159053e-01
         1
             1.541100e-01  1.000000e+00
        """
        expect = textwrap.dedent(expect).strip()

        bs2 = Basis2()
        dzvp = bs2.get_basisset("O-DZVP.H")
        self.assertEqual(str(dzvp).strip(), expect)

    def test_DZVP2(self):
        expect = """
        O-DZVP2.C
         3 2 1 0 0 0 0 0 0 0 0
         7
             5.784157e+03  8.190000e-04
             8.693035e+02  6.293500e-03
             1.985116e+02  3.178120e-02
             5.642990e+01  1.172734e-01
             1.828546e+01  3.034763e-01
             6.448715e+00  4.535214e-01
             2.341860e+00  2.430591e-01
         2
             5.459533e+00 -7.780440e-02
             4.781968e-01  5.714947e-01
         1
             1.457301e-01  1.000000e+00
         5
             3.425856e+01  5.804300e-03
             7.863895e+00  4.064030e-02
             2.344519e+00  1.550219e-01
             7.961715e-01  3.531444e-01
             2.726804e-01  4.550062e-01
         1
             8.926050e-02  1.000000e+00
         1
             6.000000e-01  1.000000e+00
        """
        expect = textwrap.dedent(expect).strip()

        bs2 = Basis2()
        dzvp2 = bs2.get_basisset("O-DZVP2.C")
        self.assertEqual(str(dzvp2).strip(), expect)


class TestBasis2Additional(unittest.TestCase):
    def test_basis_set_get_basis2(self):
        bs = BasisSet("O-TEST.H")
        cgto_s = ContractedGTO("s", 2)
        cgto_s[0] = PrimitiveGTO(50.0, 0.1)
        cgto_s[1] = PrimitiveGTO(0.5, 0.9)
        bs.append(cgto_s)

        cgto_p = ContractedGTO("p", 1)
        cgto_p[0] = PrimitiveGTO(1.5, 1.0)
        bs.append(cgto_p)

        basis2_str = bs.get_basis2()
        lines = basis2_str.strip().splitlines()
        self.assertEqual(lines[0], "O-TEST.H")
        self.assertEqual(
            lines[1].split(),
            ["1", "1", "0", "0", "0", "0", "0", "0", "0", "0", "0"],
        )

    def test_basis2_roundtrip_orbital(self):
        sample = textwrap.dedent("""\
        O-SAMPLE.H
         1 1 0 0 0 0 0 0 0 0 0
         2
             5.000000e+01  1.000000e-01
             5.000000e-01  9.000000e-01
         1
             1.500000e-01  1.000000e+00
        """)
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(sample)
            path1 = f.name
        try:
            bs2_first = Basis2(path1)
            output_str = bs2_first.get_basis2()

            with tempfile.NamedTemporaryFile("w", delete=False) as f2:
                f2.write(output_str)
                path2 = f2.name
            try:
                bs2_second = Basis2(path2)
                self.assertEqual(
                    bs2_first.get_basisset("O-SAMPLE.H"),
                    bs2_second.get_basisset("O-SAMPLE.H"),
                )
            finally:
                os.remove(path2)
        finally:
            os.remove(path1)

    def test_basis2_roundtrip_with_auxiliary(self):
        # Basis2.get_basis2() outputs the basis set name for both basis_j and basis_xc,
        # but Basis2._load() only reads the name line once before expecting J then XC data,
        # causing ValueError when re-parsing the output.
        sample = textwrap.dedent("""\
        A-SAMPLE.H
         1 0 0 0
         1
             1.000000e+00
         1 0 0 0
         1
             2.000000e+00
        """)
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(sample)
            path1 = f.name
        try:
            bs2_first = Basis2(path1)
            output_str = bs2_first.get_basis2()

            with tempfile.NamedTemporaryFile("w", delete=False) as f2:
                f2.write(output_str)
                path2 = f2.name
            try:
                bs2_second = Basis2(path2)
                self.assertEqual(
                    bs2_first.get_basisset_j("A-SAMPLE.H"),
                    bs2_second.get_basisset_j("A-SAMPLE.H"),
                )
            finally:
                os.remove(path2)
        finally:
            os.remove(path1)

    @unittest.skipUnless(
        os.environ.get("PDF_HOME")
        and os.path.exists(os.path.join(os.environ.get("PDF_HOME", ""), "data/basis2")),
        "PDF_HOME data/basis2 is not available",
    )
    def test_pdf_home_basis2_load(self):
        bs2 = Basis2()
        dzvp = bs2.get_basisset("O-DZVP.H")
        self.assertTrue(len(dzvp) > 0)


if __name__ == "__main__":
    unittest.main()
