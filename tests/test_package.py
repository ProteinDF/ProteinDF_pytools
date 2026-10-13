#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

import proteindf_tools


class TestPackage(unittest.TestCase):
    def test_import_package(self):
        self.assertIsNotNone(proteindf_tools)

    def test_version(self):
        self.assertTrue(hasattr(proteindf_tools, "__version__"))
        self.assertIsInstance(proteindf_tools.__version__, str)
        self.assertTrue(len(proteindf_tools.__version__) > 0)

    def test_all_symbols_exist_except_pdfarchive(self):
        self.assertTrue(hasattr(proteindf_tools, "__all__"))
        for name in proteindf_tools.__all__:
            if name == "PdfArchive":
                continue
            with self.subTest(symbol=name):
                self.assertTrue(
                    hasattr(proteindf_tools, name),
                    f"'{name}' in __all__ is not defined in proteindf_tools",
                )

    def test_pdfarchive_in_all_defined(self):
        """__all__ に含まれる 'PdfArchive' がパッケージ自身で定義されていること。"""
        self.assertTrue(
            hasattr(proteindf_tools, "PdfArchive"),
            "'PdfArchive' in __all__ should be defined in proteindf_tools",
        )
        self.assertIs(
            proteindf_tools.PdfArchive,
            proteindf_tools.PdfArchive_Sqlite3,
        )

    def test_wildcard_import(self):
        local_scope = {}
        exec("from proteindf_tools import *", {}, local_scope)
        self.assertIn("PdfParam", local_scope)
        self.assertIn("QmSim", local_scope)


if __name__ == "__main__":
    unittest.main()
