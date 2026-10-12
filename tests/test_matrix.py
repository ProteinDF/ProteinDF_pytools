#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import tempfile
import unittest

from proteindf_tools.matrix import Matrix, SymmetricMatrix


class TestMatrix(unittest.TestCase):
    def setUp(self):
        self.tmpfiles = []

    def tearDown(self):
        for f in self.tmpfiles:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    def _create_tempfile(self, content=b""):
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            if content:
                fp.write(content)
            path = fp.name
        self.tmpfiles.append(path)
        return path

    def test_matrix_roundtrip_row_major(self):
        # matrix_type = 0 (row-major), non-square (2 rows, 3 cols)
        m = Matrix(2, 3)
        values = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        for r in range(2):
            for c in range(3):
                m.set(r, c, values[r][c])

        tmpfile = self._create_tempfile()
        m.save(tmpfile, matrix_type=0)

        self.assertTrue(Matrix.is_loadable(tmpfile))
        self.assertEqual(Matrix.get_size(tmpfile), (2, 3))

        loaded = Matrix()
        ret = loaded.load(tmpfile)
        self.assertTrue(ret)
        self.assertEqual(loaded.rows, 2)
        self.assertEqual(loaded.cols, 3)
        for r in range(2):
            for c in range(3):
                self.assertAlmostEqual(loaded.get(r, c), values[r][c])

    def test_matrix_roundtrip_col_major(self):
        # matrix_type = 1 (col-major), non-square (3 rows, 2 cols)
        m = Matrix(3, 2)
        values = [[1.5, 2.5], [3.5, 4.5], [5.5, 6.5]]
        for r in range(3):
            for c in range(2):
                m.set(r, c, values[r][c])

        tmpfile = self._create_tempfile()
        m.save(tmpfile, matrix_type=1)

        self.assertTrue(Matrix.is_loadable(tmpfile))
        self.assertEqual(Matrix.get_size(tmpfile), (3, 2))

        loaded = Matrix()
        ret = loaded.load(tmpfile)
        self.assertTrue(ret)
        self.assertEqual(loaded.rows, 3)
        self.assertEqual(loaded.cols, 2)
        for r in range(3):
            for c in range(2):
                self.assertAlmostEqual(loaded.get(r, c), values[r][c])

    def test_matrix_find_header_struct(self):
        m = Matrix(2, 4)
        tmpfile = self._create_tempfile()
        m.save(tmpfile, matrix_type=0)

        answer, row, col, endian, header_struct, matrix_type = Matrix.find_header_struct(tmpfile)
        self.assertTrue(answer)
        self.assertEqual(row, 2)
        self.assertEqual(col, 4)
        self.assertEqual(matrix_type, 0)
        self.assertIn(endian, ["=", "<"])
        self.assertIn(header_struct, ["iii", "bii", "ill", "bll"])

    def test_symmetric_matrix_roundtrip(self):
        # SymmetricMatrix 3x3: save and load roundtrip
        sm = SymmetricMatrix(3)
        # set lower triangle elements
        sm.set(0, 0, 1.0)
        sm.set(1, 0, 2.0)
        sm.set(1, 1, 3.0)
        sm.set(2, 0, 4.0)
        sm.set(2, 1, 5.0)
        sm.set(2, 2, 6.0)

        tmpfile = self._create_tempfile()
        sm.save(tmpfile)

        # Verify file size: header (bii = 1+4+4=9 bytes) + lower triangle 3*(3+1)/2 = 6 doubles (48 bytes) = 57 bytes
        expected_size = struct.calcsize("=bii") + struct.calcsize("=d") * (3 * 4 // 2)
        self.assertEqual(os.path.getsize(tmpfile), expected_size)

        self.assertTrue(SymmetricMatrix.is_loadable(tmpfile))
        self.assertEqual(SymmetricMatrix.get_size(tmpfile), (3, 3))

        answer, row, col, endian, header_struct, matrix_type = SymmetricMatrix.find_header_struct(
            tmpfile
        )
        self.assertTrue(answer)
        self.assertEqual((row, col), (3, 3))
        self.assertEqual(matrix_type, 2)

        loaded = SymmetricMatrix()
        loaded.load(tmpfile)
        self.assertEqual(loaded.rows, 3)
        self.assertEqual(loaded.cols, 3)
        self.assertAlmostEqual(loaded.get(0, 0), 1.0)
        self.assertAlmostEqual(loaded.get(1, 0), 2.0)
        self.assertAlmostEqual(loaded.get(1, 1), 3.0)
        self.assertAlmostEqual(loaded.get(2, 0), 4.0)
        self.assertAlmostEqual(loaded.get(2, 1), 5.0)
        self.assertAlmostEqual(loaded.get(2, 2), 6.0)

    def test_matrix_load_nonexistent_file(self):
        m = Matrix()
        ret = m.load("/nonexistent_file_path_12345.bin")
        self.assertFalse(ret)

    def test_symmetric_matrix_load_nonexistent_file(self):
        # SymmetricMatrix.load logs error but does not return False when file does not exist
        # (returns None). Expected behavior according to docstring/Matrix convention is False.
        sm = SymmetricMatrix()
        ret = sm.load("/nonexistent_file_path_12345.bin")
        self.assertIs(ret, False)

    def test_matrix_load_unknown_matrix_type(self):
        # Matrix.load currently executes a bare 'raise' on unknown matrix_type,
        # which raises 'RuntimeError: No active exception to reraise'.
        # Expected behavior is raising ValueError (or returning False).
        hdr = struct.pack("<iii", 99, 2, 2)
        body = struct.pack("<4d", 1.0, 2.0, 3.0, 4.0)
        tmpfile = self._create_tempfile(content=hdr + body)
        m = Matrix()
        with self.assertRaises(ValueError):
            m.load(tmpfile)

    def test_symmetric_matrix_find_header_struct_non_square(self):
        # Non-square header (row != col) where file size happens to match row*(row+1)/2.
        # A symmetric matrix must be square (row == col). Currently, find_header_struct
        # does not check row == col and returns answer=True.
        # Expected behavior is returning answer=False.
        hdr = struct.pack("<iii", 2, 2, 3)
        body = struct.pack("<3d", 1.0, 2.0, 3.0)
        tmpfile = self._create_tempfile(content=hdr + body)
        answer, row, col, endian, header_struct, matrix_type = (
            SymmetricMatrix.find_header_struct(tmpfile)
        )
        self.assertFalse(answer)

    def test_symmetric_matrix_is_loadable_non_square_matrix_file(self):
        # A standard non-square Matrix file (e.g. 2x3, 6 doubles)
        # is not loadable as SymmetricMatrix
        m = Matrix(2, 3)
        tmpfile = self._create_tempfile()
        m.save(tmpfile)
        self.assertFalse(SymmetricMatrix.is_loadable(tmpfile))

    def test_is_loadable_wrong_size(self):
        tmpfile = self._create_tempfile(content=b"\x00" * 100)
        self.assertFalse(Matrix.is_loadable(tmpfile))
        self.assertFalse(SymmetricMatrix.is_loadable(tmpfile))

    def test_matrix_crafted_endian(self):
        # Craft binary files with little-endian and big-endian headers and data
        for endian_name, prefix in [("little", "<"), ("big", ">")]:
            for hfmt in ["iii", "bii", "ill", "bll"]:
                for mtype in [0, 1]:
                    row, col = 2, 3
                    if mtype == 0:
                        # row-major order: (0,0), (0,1), (0,2), (1,0), (1,1), (1,2)
                        raw_data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
                    else:
                        # col-major order: (0,0), (1,0), (0,1), (1,1), (0,2), (1,2)
                        raw_data = [1.0, 4.0, 2.0, 5.0, 3.0, 6.0]

                    hdr = struct.pack(prefix + hfmt, mtype, row, col)
                    body = struct.pack(prefix + f"{len(raw_data)}d", *raw_data)
                    tmpfile = self._create_tempfile(content=hdr + body)

                    with self.subTest(endian=endian_name, hfmt=hfmt, mtype=mtype):
                        self.assertTrue(Matrix.is_loadable(tmpfile))
                        self.assertEqual(Matrix.get_size(tmpfile), (2, 3))
                        m = Matrix()
                        self.assertTrue(m.load(tmpfile))
                        self.assertEqual(m.rows, 2)
                        self.assertEqual(m.cols, 3)
                        expected = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
                        for r in range(2):
                            for c in range(3):
                                self.assertAlmostEqual(m.get(r, c), expected[r][c])

    def test_symmetric_matrix_crafted_endian(self):
        # Craft binary files for SymmetricMatrix with little-endian and big-endian
        for endian_name, prefix in [("little", "<"), ("big", ">")]:
            for hfmt in ["iii", "bii", "ill", "bll"]:
                dim = 3
                # lower triangle: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2)
                raw_data = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
                hdr = struct.pack(prefix + hfmt, 2, dim, dim)
                body = struct.pack(prefix + f"{len(raw_data)}d", *raw_data)
                tmpfile = self._create_tempfile(content=hdr + body)

                with self.subTest(endian=endian_name, hfmt=hfmt):
                    self.assertTrue(SymmetricMatrix.is_loadable(tmpfile))
                    self.assertEqual(SymmetricMatrix.get_size(tmpfile), (3, 3))
                    sm = SymmetricMatrix()
                    sm.load(tmpfile)
                    self.assertEqual(sm.rows, 3)
                    self.assertEqual(sm.cols, 3)
                    expected = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
                    idx = 0
                    for r in range(dim):
                        for c in range(r + 1):
                            self.assertAlmostEqual(sm.get(r, c), expected[idx])
                            idx += 1


if __name__ == "__main__":
    unittest.main()
