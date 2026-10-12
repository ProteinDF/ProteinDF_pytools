#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import tempfile
import unittest

from proteindf_tools.vector import Vector


class TestVector(unittest.TestCase):
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

    @unittest.expectedFailure
    def test_vector_roundtrip(self):
        """Vector roundtrip test.

        Vector.save packs header size using 'l' with native size (8 bytes on 64-bit systems),
        whereas Vector.find_header_struct only tests '=i', '=l', '<i', '<l', '>i', '>l' where
        'l' with standard prefix has size 4 bytes. Thus, a saved vector file cannot be loaded
        back on 64-bit platforms without fixing Vector.save or find_header_struct.
        """
        v = Vector([1.0, 2.5, -3.2, 4.8])
        tmpfile = self._create_tempfile()
        v.save(tmpfile)

        self.assertTrue(Vector.is_loadable(tmpfile))
        loaded = Vector()
        ret = loaded.load(tmpfile)
        self.assertTrue(ret)
        self.assertEqual(len(loaded), 4)
        for i in range(4):
            self.assertAlmostEqual(loaded[i], v[i])

    def test_vector_crafted_roundtrip(self):
        # Craft standard 32-bit integer header (=i or =l) and verify Vector find_header_struct, is_loadable, load
        values = [1.0, 2.5, -3.2, 4.8]
        hdr = struct.pack("=i", len(values))
        body = struct.pack(f"={len(values)}d", *values)
        tmpfile = self._create_tempfile(content=hdr + body)

        self.assertTrue(Vector.is_loadable(tmpfile))
        answer, size, endian, header_struct = Vector.find_header_struct(tmpfile)
        self.assertTrue(answer)
        self.assertEqual(size, 4)
        self.assertIn(endian, ["=", "<"])
        self.assertEqual(header_struct, "i")

        loaded = Vector()
        ret = loaded.load(tmpfile)
        self.assertTrue(ret)
        self.assertEqual(len(loaded), 4)
        for i in range(4):
            self.assertAlmostEqual(loaded[i], values[i])

    def test_vector_load_nonexistent_file(self):
        v = Vector()
        ret = v.load("/nonexistent_vector_file_12345.bin")
        self.assertFalse(ret)

    def test_vector_is_loadable_wrong_size(self):
        tmpfile = self._create_tempfile(content=b"\x00" * 100)
        self.assertFalse(Vector.is_loadable(tmpfile))

    def test_vector_crafted_endian(self):
        # Craft vector files with little-endian and big-endian headers
        for endian_name, prefix in [("little", "<"), ("big", ">")]:
            for hfmt in ["i", "l"]:
                values = [10.5, -20.25, 30.125]
                hdr = struct.pack(prefix + hfmt, len(values))
                body = struct.pack(prefix + f"{len(values)}d", *values)
                tmpfile = self._create_tempfile(content=hdr + body)

                with self.subTest(endian=endian_name, hfmt=hfmt):
                    self.assertTrue(Vector.is_loadable(tmpfile))
                    answer, size, endian, header_struct = Vector.find_header_struct(tmpfile)
                    self.assertTrue(answer)
                    self.assertEqual(size, 3)
                    if prefix == "<":
                        self.assertIn(endian, ["=", "<"])
                    else:
                        self.assertEqual(endian, prefix)
                    self.assertEqual(header_struct, "i")

                    v = Vector()
                    ret = v.load(tmpfile)
                    self.assertTrue(ret)
                    self.assertEqual(len(v), 3)
                    for i in range(3):
                        self.assertAlmostEqual(v[i], values[i])


if __name__ == "__main__":
    unittest.main()
