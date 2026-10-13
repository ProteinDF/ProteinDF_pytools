#!/usr/bin/env python
# -*- coding: utf-8 -*-

import contextlib
import io
import os
import tempfile
import unittest

import proteindf_bridge as bridge
from proteindf_tools.poputils import PopEspUtils, PopUtils


class TestPopUtils(unittest.TestCase):
    def test_get_atomlist(self):
        mol = bridge.AtomGroup()
        h1 = bridge.Atom()
        h1.symbol = "H"
        mol.set_atom("1", h1)

        h2 = bridge.Atom()
        h2.symbol = "H"
        mol.set_atom("2", h2)

        x = bridge.Atom()
        x.symbol = "X"
        mol.set_atom("3", x)

        atomlist, charge = PopUtils._get_atomlist(mol)
        symbols = [a.symbol for a in atomlist]
        self.assertEqual(symbols, ["H", "H"])
        self.assertAlmostEqual(charge, 0.0)

    def test_get_atomlist_by_param(self):
        raw = {
            "control": {},
            "guess": "none",
            "coordinates": {"groups": {}, "atoms": [{"symbol": "H", "xyz": [0.0, 0.0, 0.0]}]},
            "num_of_AOs": 1,
            "num_of_MOs": 1,
            "num_of_atoms": 1,
            "num_of_iterations": 1,
        }
        with tempfile.NamedTemporaryFile(delete=False) as f:
            path = f.name
        try:
            bridge.save_msgpack(raw, path)
            atomlist, charge = PopUtils.get_atomlist_by_param(path)
            self.assertEqual(len(atomlist), 1)
            self.assertEqual(atomlist[0].symbol, "H")
            self.assertAlmostEqual(charge, 0.0)
        finally:
            if os.path.exists(path):
                os.remove(path)

    def test_set_atomlist(self):
        a1 = bridge.Atom()
        a2 = bridge.Atom()
        stdout_buf = io.StringIO()
        with contextlib.redirect_stdout(stdout_buf):
            PopUtils.set_atomlist([0.3, -0.3], [a1, a2])
        self.assertAlmostEqual(a1.charge, 0.3)
        self.assertAlmostEqual(a2.charge, -0.3)
        self.assertIn("total charges:  0.000", stdout_buf.getvalue())


class TestPopEspUtils(unittest.TestCase):
    def setUp(self):
        self.esp_util = PopEspUtils()

    def test_calc_esp(self):
        a1 = bridge.Atom()
        a1.xyz = bridge.Position(0.0, 0.0, 0.0)
        a1.charge = 2.0

        pos = bridge.Position(0.0, 0.0, 2.0)
        esp = self.esp_util._calc_esp(pos, [a1])
        # esp = 2.0 / 2.0 = 1.0
        self.assertAlmostEqual(esp, 1.0)

        a2 = bridge.Atom()
        a2.xyz = bridge.Position(0.0, 0.0, 4.0)
        a2.charge = -1.0
        esp2 = self.esp_util._calc_esp(pos, [a1, a2])
        # esp = 2.0 / 2.0 + (-1.0) / 2.0 = 0.5
        self.assertAlmostEqual(esp2, 0.5)

    def test_calc_rrms(self):
        a1 = bridge.Atom()
        a1.xyz = bridge.Position(0.0, 0.0, 0.0)
        a1.charge = 2.0

        grids = [bridge.Position(0.0, 0.0, 2.0)]
        esps = [2.0]

        # estimate = 1.0, exact = 2.0 -> delta = -1.0
        # sum_delta2 = 1.0, sum_v2 = 4.0 -> rrms = sqrt(1.0 / 4.0) = 0.5
        rrms = self.esp_util._calcRRMS(grids, esps, [a1])
        self.assertAlmostEqual(rrms, 0.5)

    def test_calc_rrms_exact_match(self):
        a1 = bridge.Atom()
        a1.xyz = bridge.Position(0.0, 0.0, 0.0)
        a1.charge = 2.0

        grids = [bridge.Position(0.0, 0.0, 2.0)]
        esps = [1.0]

        # estimate = 1.0, exact = 1.0 -> delta = 0.0 -> rrms = 0.0
        rrms = self.esp_util._calcRRMS(grids, esps, [a1])
        self.assertAlmostEqual(rrms, 0.0)


if __name__ == "__main__":
    unittest.main()
