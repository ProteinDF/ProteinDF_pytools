#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from unittest.mock import patch

import proteindf_bridge as bridge
from proteindf_tools.pdfcommon import load_pdfparam, pdf_home, save_pdfparam, set_basisset
from proteindf_tools.pdfparam import PdfParam


class TestPdfCommon(unittest.TestCase):
    def test_pdf_home_set(self):
        with patch.dict(os.environ, {"PDF_HOME": "/custom/pdf/home"}):
            self.assertEqual(pdf_home(), "/custom/pdf/home")

    def test_pdf_home_unset(self):
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(pdf_home(), "")

    def test_save_load_pdfparam_roundtrip(self):
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
        with tempfile.NamedTemporaryFile(delete=False) as f1:
            path1 = f1.name
        with tempfile.NamedTemporaryFile(delete=False) as f2:
            path2 = f2.name
        try:
            bridge.save_msgpack(raw, path1)
            p1 = load_pdfparam(path1)
            self.assertIsInstance(p1, PdfParam)
            self.assertEqual(p1.method, "roks")
            self.assertEqual(p1.work_path, "/tmp/raw_test")
            self.assertEqual(p1.iterations, 7)

            save_pdfparam(p1, path2)
            p2 = load_pdfparam(path2)
            self.assertIsInstance(p2, PdfParam)
            self.assertEqual(p2.method, "roks")
            self.assertEqual(p2.work_path, "/tmp/raw_test")
            self.assertEqual(p2.iterations, 7)
        finally:
            for p in (path1, path2):
                if os.path.exists(p):
                    os.remove(p)

    def test_save_pdfparam_non_string_path(self):
        p = PdfParam()
        with self.assertRaises(AssertionError):
            save_pdfparam(p, 12345)

    @unittest.skipUnless(
        os.environ.get("PDF_HOME")
        and os.path.exists(os.path.join(os.environ.get("PDF_HOME", ""), "data/basis2")),
        "PDF_HOME data/basis2 is not available",
    )
    def test_set_basisset(self):
        """set_basisset(pdfparam) は基底関数を設定した PdfParam を返すべき。

        実際には pdfcommon.py 内で Basis2 がインポートされておらず、
        NameError: name 'Basis2' is not defined が発生する。
        """
        p = PdfParam()
        result = set_basisset(p)
        self.assertIsInstance(result, PdfParam)


if __name__ == "__main__":
    unittest.main()
