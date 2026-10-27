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

    def test_all_symbols_exist(self):
        self.assertTrue(hasattr(proteindf_tools, "__all__"))
        for name in proteindf_tools.__all__:
            with self.subTest(symbol=name):
                self.assertTrue(
                    hasattr(proteindf_tools, name),
                    f"'{name}' in __all__ is not defined in proteindf_tools",
                )

    def test_wildcard_import(self):
        local_scope = {}
        exec("from proteindf_tools import *", {}, local_scope)
        self.assertIn("PdfParam", local_scope)
        self.assertIn("QmSim", local_scope)


if __name__ == "__main__":
    unittest.main()
