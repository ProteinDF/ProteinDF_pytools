#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import os
from pathlib import Path
import unittest
from unittest import mock

# scripts/pdf-test-eri.py contains hyphens in the filename, so load via importlib
SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "pdf-test-eri.py"
spec = importlib.util.spec_from_file_location("pdf_test_eri", SCRIPT_PATH)
pdf_test_eri = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pdf_test_eri)


class TestPdfTestEri(unittest.TestCase):
    def test_gau2pdf_corresponding_index(self):
        # Basis type 4 adds +4
        self.assertEqual(pdf_test_eri.gau2pdf_corresponding_index(10, 4), 14)
        # Basis type 8 subtracts -4
        self.assertEqual(pdf_test_eri.gau2pdf_corresponding_index(10, 8), 6)
        # Other basis types remain unchanged
        for bt in [0, 1, 2, 3, 5, 6, 7, 9]:
            self.assertEqual(pdf_test_eri.gau2pdf_corresponding_index(10, bt), 10)

    def test_translate_5d_index(self):
        orbinfo = {
            4: {"basis_type_id": 4},
            5: {"basis_type_id": 5},
            6: {"basis_type_id": 6},
            7: {"basis_type_id": 7},
            8: {"basis_type_id": 8},
        }
        # In 5d basis, orbital 4 (type 4) becomes 8 (+4), orbital 8 (type 8) becomes 4 (-4),
        # and orbitals 5-7 (types 5-7) remain unchanged.
        res = pdf_test_eri.translate_5d_index(orbinfo, 4, 5, 6, 8)
        self.assertEqual(res, (8, 5, 6, 4))

    @mock.patch("sys.stderr")
    def test_check_eri(self, mock_stderr):
        # Clearly within threshold -> True
        self.assertTrue(pdf_test_eri.check_eri(1.0, 1.0001, 1e-3, 0, 0, 0, 0))
        # Clearly exceeds threshold -> False
        self.assertFalse(pdf_test_eri.check_eri(1.0, 1.01, 1e-3, 0, 0, 0, 0))

    @mock.patch.dict(os.environ, {"PDF_HOME": "/dummy/pdf/home"})
    @mock.patch("builtins.print")
    @mock.patch.object(pdf_test_eri.bridge, "load_msgpack")
    @mock.patch("subprocess.Popen")
    def test_make_orbinfo_table(self, mock_popen, mock_load, mock_print):
        mock_proc = mock.MagicMock()
        mock_proc.stdout.readline.side_effect = ["header line\n", ""]
        mock_proc.wait.return_value = 0
        mock_popen.return_value = mock_proc
        mock_load.return_value = {0: {"basis_type_id": 0}}

        res = pdf_test_eri.make_orbinfo_table("param.in")

        mock_popen.assert_called_once()
        args, kwargs = mock_popen.call_args
        self.assertEqual(
            args[0],
            ["/dummy/pdf/home/bin/pdf-info-orb", "-p", "param.in", "-w", "orbinfo.mpac"],
        )
        self.assertTrue(kwargs.get("universal_newlines"))
        mock_load.assert_called_once_with("orbinfo.mpac")
        self.assertEqual(res, {0: {"basis_type_id": 0}})

    @mock.patch.dict(os.environ, {"PDF_HOME": "/dummy/pdf/home"})
    @mock.patch("subprocess.Popen")
    def test_get_pdf_eri_list(self, mock_popen):
        mock_proc = mock.MagicMock()
        mock_proc.stdout.readline.side_effect = [
            "(0 0|0 0)=  1.5\n",
            "dummy non-matching line\n",
            "(1 2|3 4)= -0.12345\n",
            "",
        ]
        mock_proc.wait.return_value = 0
        mock_popen.return_value = mock_proc

        eri_list = pdf_test_eri.get_pdf_eri_list("param.in", "indices.mpac")
        self.assertEqual(eri_list, [1.5, -0.12345])

        mock_popen.assert_called_once()
        args, kwargs = mock_popen.call_args
        self.assertEqual(
            args[0],
            ["/dummy/pdf/home/bin/pdf-eri", "-p", "param.in", "indices.mpac"],
        )
        self.assertTrue(kwargs.get("universal_newlines"))

    @mock.patch("sys.stderr")
    @mock.patch.object(pdf_test_eri, "get_pdf_eri_list")
    @mock.patch.object(pdf_test_eri, "make_orbinfo_table")
    @mock.patch.object(pdf_test_eri.bridge, "save_msgpack")
    @mock.patch.object(pdf_test_eri.bridge, "load_msgpack")
    def test_main_success(
        self, mock_load, mock_save, mock_make_orb, mock_get_eri, mock_stderr
    ):
        mock_load.return_value = [[4, 5, 6, 8, 1.5]]
        mock_make_orb.return_value = {
            4: {"basis_type_id": 4},
            5: {"basis_type_id": 5},
            6: {"basis_type_id": 6},
            8: {"basis_type_id": 8},
        }
        mock_get_eri.return_value = [1.5001]

        test_args = ["pdf-test-eri.py", "dummy_eri.mpac", "dummy_param.in", "-t", "1e-3"]
        with mock.patch("sys.argv", test_args):
            pdf_test_eri.main()

        mock_save.assert_called_once_with([(8, 5, 6, 4)], pdf_test_eri.PDF_INDEX_FILE)

    @mock.patch("sys.stderr")
    @mock.patch.object(pdf_test_eri, "get_pdf_eri_list")
    @mock.patch.object(pdf_test_eri, "make_orbinfo_table")
    @mock.patch.object(pdf_test_eri.bridge, "save_msgpack")
    @mock.patch.object(pdf_test_eri.bridge, "load_msgpack")
    def test_main_failure(
        self, mock_load, mock_save, mock_make_orb, mock_get_eri, mock_stderr
    ):
        mock_load.return_value = [[0, 0, 0, 0, 1.5]]
        mock_make_orb.return_value = {0: {"basis_type_id": 0}}
        mock_get_eri.return_value = [1.6]

        test_args = ["pdf-test-eri.py", "dummy_eri.mpac", "dummy_param.in", "-t", "1e-3"]
        with mock.patch("sys.argv", test_args):
            with self.assertRaises(SystemExit) as cm:
                pdf_test_eri.main()
        self.assertIn("1/1 ERROR FOUND", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
