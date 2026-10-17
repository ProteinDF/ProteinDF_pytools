#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import io
import os
import pickle
import re
import sys
import tempfile
import unittest
import unittest.mock as mock
from contextlib import redirect_stdout

import proteindf_bridge as bridge
import proteindf_tools as pdf

# Setup module aliases for legacy scripts importing pdfbridge / pdfpytools
if "pdfbridge" not in sys.modules:
    sys.modules["pdfbridge"] = bridge
if "pdfpytools" not in sys.modules:
    sys.modules["pdfpytools"] = pdf

# Mock or provide minimal jinja2 if not installed
try:
    import jinja2
except ImportError:

    def _render_simple_template(tpl, **kwargs):
        # Match loop: {%-? for item in items %}(body){%-? endfor %}
        for_pattern = re.compile(
            r"\{%-?\s*for\s+(\w+)\s+in\s+(\w+)\s*-?%\}(.*?)\{%-?\s*endfor\s*-?%\}",
            re.DOTALL,
        )

        def repl_for(m):
            item_var = m.group(1)
            items_var = m.group(2)
            body = m.group(3)
            items = kwargs.get(items_var, [])
            out = []
            for it in items:
                cur = body

                def repl_var(vm):
                    expr = vm.group(1).strip()
                    m_fmt = re.match(r'\"([^\"]+)\"\|format\((\w+)\.(\w+)\)', expr)
                    if m_fmt:
                        fmt, obj, attr = m_fmt.group(1), m_fmt.group(2), m_fmt.group(3)
                        if obj == item_var:
                            val = it[attr] if isinstance(it, dict) else getattr(it, attr)
                            return fmt % val
                    m_dot = re.match(r"(\w+)\.(\w+)", expr)
                    if m_dot:
                        obj, attr = m_dot.group(1), m_dot.group(2)
                        if obj == item_var:
                            val = it[attr] if isinstance(it, dict) else getattr(it, attr)
                            return str(val)
                    return vm.group(0)

                cur = re.sub(r"\{\{\s*(.*?)\s*\}\}", repl_var, cur)
                out.append(cur)
            return "".join(out)

        res = for_pattern.sub(repl_for, tpl)

        def repl_top_var(vm):
            expr = vm.group(1).strip()
            m_fmt = re.match(r'\"([^\"]+)\"\|format\((.*?)\)', expr)
            if m_fmt:
                fmt, subexpr = m_fmt.group(1), m_fmt.group(2).strip()
                try:
                    val = eval(subexpr, {}, kwargs)
                    return fmt % val
                except Exception:
                    return ""
            if expr in kwargs:
                val = kwargs[expr]
                return "" if val is None else str(val)
            return ""

        res = re.sub(r"\{\{\s*(.*?)\s*\}\}", repl_top_var, res)
        return res

    class FakeTemplate(object):
        def __init__(self, templ_str):
            self.templ_str = templ_str

        def render(self, **kwargs):
            return _render_simple_template(self.templ_str, **kwargs)

    class FakeFileSystemLoader(object):
        def __init__(self, searchpath="./", encoding="utf8"):
            self.searchpath = searchpath

    class FakeEnvironment(object):
        def __init__(self, loader=None):
            self.loader = loader

        def from_string(self, source):
            return FakeTemplate(source)

        def get_template(self, name):
            with open(name, "r", encoding="utf-8") as f:
                return FakeTemplate(f.read())

    jinja2_mock = mock.MagicMock()
    jinja2_mock.Environment = FakeEnvironment
    jinja2_mock.FileSystemLoader = FakeFileSystemLoader
    sys.modules["jinja2"] = jinja2_mock

# Mock requests and bs4 if not available
try:
    import requests
except ImportError:
    requests_mock = mock.MagicMock()
    sys.modules["requests"] = requests_mock

try:
    import bs4
except ImportError:
    bs4_mock = mock.MagicMock()
    sys.modules["bs4"] = bs4_mock


SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


def _load_script(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    module_name = filename.replace("-", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestXyz2Pdf(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_xyz2pdf_normal_atom_order_and_coords(self):
        """xyz2pdf preserves atom ordering and Angstrom coordinates for CH3CHO."""
        # 7-atom system: CH3CHO
        xyz_content = """7
Acetaldehyde (CH3CHO)
C  -0.38400   0.40700   0.00000
C   1.04700  -0.08900   0.00000
O   1.38500  -1.25000   0.00000
H  -0.45700   1.49300   0.00000
H  -0.90200   0.03800   0.88700
H  -0.90200   0.03800  -0.88700
H   1.81000   0.71800   0.00000
"""
        xyz_path = os.path.join(self.tmpdir.name, "input.xyz")
        with open(xyz_path, "w") as f:
            f.write(xyz_content)

        mod = _load_script("xyz2pdf.py")
        # Run with default output ("fl_Userinput" in current working directory)
        sys.argv = ["xyz2pdf.py", xyz_path]
        mod.main()

        out_path = "fl_Userinput"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()

        # Check unit specification (Angstrom, coordinates unchanged)
        self.assertIn("geometry/cartesian/unit = angstrom", content)

        # Check number of electrons (C:6*2 + O:8 + H:1*4 = 24 electrons)
        self.assertIn("method/rks/electron-number      = 24", content)
        self.assertIn("method/rks/occlevel     = [ 1 - 12 ]", content)

        # Extract atoms block and verify exact ordering and coordinates
        m = re.search(r"geometry/cartesian/input\s*=\s*\{(.*?)\}end", content, re.DOTALL)
        self.assertIsNotNone(m)
        lines = [line.strip() for line in m.group(1).strip().splitlines() if line.strip()]
        self.assertEqual(len(lines), 7)

        expected_atoms = [
            ("C", -0.384, 0.407, 0.0),
            ("C", 1.047, -0.089, 0.0),
            ("O", 1.385, -1.250, 0.0),
            ("H", -0.457, 1.493, 0.0),
            ("H", -0.902, 0.038, 0.887),
            ("H", -0.902, 0.038, -0.887),
            ("H", 1.810, 0.718, 0.0),
        ]
        for line, (sym, x, y, z) in zip(lines, expected_atoms):
            tokens = line.split()
            self.assertEqual(tokens[0], sym)
            self.assertAlmostEqual(float(tokens[1]), x, places=4)
            self.assertAlmostEqual(float(tokens[2]), y, places=4)
            self.assertAlmostEqual(float(tokens[3]), z, places=4)

    def test_xyz2pdf_water_normal(self):
        """xyz2pdf preserves ordering and calculates electrons for water."""
        xyz_content = """3
Water molecule
O   0.00000   0.00000   0.11700
H   0.00000   0.75700  -0.46900
H   0.00000  -0.75700  -0.46900
"""
        xyz_path = os.path.join(self.tmpdir.name, "water.xyz")
        with open(xyz_path, "w") as f:
            f.write(xyz_content)

        mod = _load_script("xyz2pdf.py")
        sys.argv = ["xyz2pdf.py", xyz_path]
        mod.main()

        out_path = "fl_Userinput"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()

        # O(8) + 2*H(1) = 10 electrons
        self.assertIn("method/rks/electron-number      = 10", content)
        self.assertIn("method/rks/occlevel     = [ 1 - 5 ]", content)

        m = re.search(r"geometry/cartesian/input\s*=\s*\{(.*?)\}end", content, re.DOTALL)
        lines = [l.strip() for l in m.group(1).strip().splitlines() if l.strip()]
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("O"))
        self.assertTrue(lines[1].startswith("H"))
        self.assertTrue(lines[2].startswith("H"))

    @unittest.expectedFailure
    def test_xyz2pdf_custom_output_option(self):
        """xyz2pdf writes output to the path specified by -o option.

        scripts/xyz2pdf.py:101: args.output[0] で文字列の1文字目だけが
        ファイル名として扱われてしまう不具合があるため失敗する。
        修正されたら外す。
        """
        xyz_content = """1
H
H 0.0 0.0 0.0
"""
        xyz_path = os.path.join(self.tmpdir.name, "h.xyz")
        with open(xyz_path, "w") as f:
            f.write(xyz_content)

        custom_out = os.path.join(self.tmpdir.name, "custom_output.inp")
        mod = _load_script("xyz2pdf.py")
        sys.argv = ["xyz2pdf.py", xyz_path, "-o", custom_out]
        mod.main()

        self.assertTrue(os.path.exists(custom_out))
        with open(custom_out, "r") as f:
            content = f.read()
        self.assertIn("geometry/cartesian/unit = angstrom", content)

    def test_xyz2pdf_empty_file(self):
        """xyz2pdf raises ValueError when loading an empty xyz file."""
        empty_path = os.path.join(self.tmpdir.name, "empty.xyz")
        with open(empty_path, "w") as f:
            pass

        mod = _load_script("xyz2pdf.py")
        sys.argv = ["xyz2pdf.py", empty_path]
        with self.assertRaises(ValueError):
            mod.main()

    def test_xyz2pdf_atom_count_mismatch(self):
        """xyz2pdf raises IndexError when atom count header does not match content."""
        mismatch_path = os.path.join(self.tmpdir.name, "mismatch.xyz")
        with open(mismatch_path, "w") as f:
            f.write("5\ntitle\nH 0.0 0.0 0.0\n")

        mod = _load_script("xyz2pdf.py")
        sys.argv = ["xyz2pdf.py", mismatch_path]
        with self.assertRaises(IndexError):
            mod.main()

    def test_xyz2pdf_nonexistent_file(self):
        """bridge.Xyz().load() does not raise error on missing file and creates 0 atoms."""
        nonexistent_path = os.path.join(self.tmpdir.name, "nonexistent.xyz")

        mod = _load_script("xyz2pdf.py")
        sys.argv = ["xyz2pdf.py", nonexistent_path]
        mod.main()

        out_path = "fl_Userinput"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()
        self.assertIn("method/rks/electron-number      = 0", content)


class TestXyz2Gau(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    @unittest.expectedFailure
    def test_xyz2gau_default_template(self):
        """xyz2gau creates Gaussian input using default template.

        scripts/xyz2gau.py:77: g09_templ_str が未定義(正しくは templ_str)のため
        NameError が発生する。
        修正されたら外す。
        """
        xyz_content = """1
Hydrogen
H   0.00000 0.00000 0.00000
"""
        xyz_path = os.path.join(self.tmpdir.name, "h.xyz")
        with open(xyz_path, "w") as f:
            f.write(xyz_content)

        mod = _load_script("xyz2gau.py")
        sys.argv = ["xyz2gau.py", xyz_path]
        mod.main()

        out_path = "sample.gau"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()
        self.assertIn("H    0.00000  0.00000  0.00000", content)

    @unittest.expectedFailure
    def test_xyz2gau_custom_template(self):
        """xyz2gau creates Gaussian input using custom template specified by -t.

        scripts/xyz2gau.py:77: -t を指定した場合でも77行目の g09_templ_str 未定義により
        NameError が発生する。
        修正されたら外す。
        """
        xyz_content = """1
Hydrogen
H   0.00000 0.00000 0.00000
"""
        xyz_path = os.path.join(self.tmpdir.name, "h.xyz")
        with open(xyz_path, "w") as f:
            f.write(xyz_content)

        tpl_content = """# custom gaussian template
{%- for atom in atomlist %}
{{ atom.symbol }} {{ atom.x }} {{ atom.y }} {{ atom.z }}
{%- endfor %}
"""
        tpl_path = os.path.join(self.tmpdir.name, "custom.tpl")
        with open(tpl_path, "w") as f:
            f.write(tpl_content)

        mod = _load_script("xyz2gau.py")
        sys.argv = ["xyz2gau.py", xyz_path, "-t", tpl_path]
        mod.main()

        out_path = "sample.gau"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()
        self.assertIn("# custom gaussian template", content)
        self.assertIn("H 0.0 0.0 0.0", content)

    def test_xyz2gau_empty_file(self):
        """xyz2gau raises ValueError on empty xyz file."""
        empty_path = os.path.join(self.tmpdir.name, "empty.xyz")
        with open(empty_path, "w") as f:
            pass

        mod = _load_script("xyz2gau.py")
        sys.argv = ["xyz2gau.py", empty_path]
        with self.assertRaises(ValueError):
            mod.main()

    def test_xyz2gau_atom_count_mismatch(self):
        """xyz2gau raises IndexError on mismatched xyz file."""
        mismatch_path = os.path.join(self.tmpdir.name, "mismatch.xyz")
        with open(mismatch_path, "w") as f:
            f.write("3\ntitle\nH 0.0 0.0 0.0\n")

        mod = _load_script("xyz2gau.py")
        sys.argv = ["xyz2gau.py", mismatch_path]
        with self.assertRaises(IndexError):
            mod.main()


class TestPdfSetCharges(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def _create_synthetic_param(self, path):
        param = pdf.PdfParam()
        param.comment = "synthetic test system"
        param.num_of_atoms = 2
        param.num_of_AOs = 2
        param.num_of_MOs = 2
        param.method = "rks"
        param.iterations = 1
        param.guess = "none"

        mol = bridge.AtomGroup()
        a1 = bridge.Atom(symbol="H", position=bridge.Position([0.0, 0.0, 0.0]), charge=0.0, name="H1")
        a2 = bridge.Atom(symbol="O", position=bridge.Position([0.0, 0.0, 1.0]), charge=0.0, name="O1")
        mol.set_atom("0", a1)
        mol.set_atom("1", a2)
        param.molecule = mol

        param._data.setdefault("control", {})["file_base_name"] = {}
        pdf.save_pdfparam(param, path)
        return param

    @unittest.expectedFailure
    def test_pdf_set_charges_normal(self):
        """pdf-set-charges.py sets charges from CSV to atoms in pdfparam.

        scripts/pdf-set-charges.py:123: count が初期化・定義されていないため
        NameError が発生する。
        修正されたら外す。
        """
        param_path = os.path.join(self.tmpdir.name, "pdfparam.mpac")
        csv_path = os.path.join(self.tmpdir.name, "charges.csv")
        out_path = os.path.join(self.tmpdir.name, "new_pdfparam.mpac")

        self._create_synthetic_param(param_path)
        with open(csv_path, "w") as f:
            f.write("0.4\n-0.8\n")

        mod = _load_script("pdf-set-charges.py")
        sys.argv = [
            "pdf-set-charges.py",
            "-i",
            param_path,
            "-c",
            csv_path,
            "-o",
            out_path,
        ]
        with redirect_stdout(io.StringIO()):
            mod.main()

        self.assertTrue(os.path.exists(out_path))
        new_param = pdf.load_pdfparam(out_path)
        atoms = new_param.molecule.get_atom_list()
        self.assertAlmostEqual(atoms[0].charge, 0.4)
        self.assertAlmostEqual(atoms[1].charge, -0.8)

    def test_pdf_set_charges_atom_count_mismatch(self):
        """pdf-set-charges.py fails with AssertionError if CSV rows < atom count."""
        param_path = os.path.join(self.tmpdir.name, "pdfparam.mpac")
        csv_path = os.path.join(self.tmpdir.name, "charges.csv")

        self._create_synthetic_param(param_path)
        # 1 charge for 2 atoms
        with open(csv_path, "w") as f:
            f.write("0.5\n")

        mod = _load_script("pdf-set-charges.py")
        sys.argv = ["pdf-set-charges.py", "-i", param_path, "-c", csv_path]
        with redirect_stdout(io.StringIO()):
            with self.assertRaises(AssertionError):
                mod.main()

    def test_pdf_set_charges_invalid_csv(self):
        """pdf-set-charges.py fails with ValueError if CSV contains non-numeric rows."""
        csv_path = os.path.join(self.tmpdir.name, "invalid.csv")
        with open(csv_path, "w") as f:
            f.write("not_a_number\n")

        mod = _load_script("pdf-set-charges.py")
        with redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                mod.load_charges(csv_path)


class TestPdfTransBasis(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    @unittest.expectedFailure
    def test_pdf_trans_basis_gaussian94_normal(self):
        """pdf-trans-basis.py parses Gaussian94 basis set file.

        scripts/pdf-trans-basis.py:73: _parse_block_gaussian94 が
        _parse_block_aussian94 と誤記されており AttributeError が発生する
        (また ContractedGTO 等の未定義もある)。
        修正されたら外す。
        """
        bs_content = """****
H     0
S   1   1.00
      3.42525091           1.00000000
****
"""
        bs_path = os.path.join(self.tmpdir.name, "sto3g.g94")
        with open(bs_path, "w") as f:
            f.write(bs_content)

        mod = _load_script("pdf-trans-basis.py")
        parser = mod.BasisSetParser()
        parser.load(bs_path)

        self.assertIn("H", parser._basissets)

    def test_pdf_trans_basis_nonexistent_file(self):
        """pdf-trans-basis.py does not fail and leaves basissets empty on nonexistent file."""
        nonexistent = os.path.join(self.tmpdir.name, "none.g94")
        mod = _load_script("pdf-trans-basis.py")
        parser = mod.BasisSetParser()
        parser.load(nonexistent)
        self.assertEqual(len(parser._basissets), 0)

    def test_pdf_trans_basis_empty_file(self):
        """pdf-trans-basis.py leaves basissets empty on file without sections."""
        empty_path = os.path.join(self.tmpdir.name, "empty.g94")
        with open(empty_path, "w") as f:
            f.write("no section marker\n")

        mod = _load_script("pdf-trans-basis.py")
        parser = mod.BasisSetParser()
        parser.load(empty_path)
        self.assertEqual(len(parser._basissets), 0)


class TestPdfMakeMat(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def test_pdf_make_mat_generation_and_values(self):
        """pdf-make-mat.py generates A.mat and B.mat with expected row/col ordering."""
        mod = _load_script("pdf-make-mat.py")
        sys.argv = ["pdf-make-mat.py"]
        mod.main()

        self.assertTrue(os.path.exists("A.mat"))
        self.assertTrue(os.path.exists("B.mat"))

        mat_a = pdf.Matrix()
        mat_a.load("A.mat")
        self.assertEqual((mat_a.rows, mat_a.cols), (3, 3))
        # Row-major: A.set(r, c, count) with count += 1 in inner c loop
        count = 0.0
        for r in range(3):
            for c in range(3):
                self.assertAlmostEqual(mat_a.get(r, c), count)
                count += 1.0

        mat_b = pdf.Matrix()
        mat_b.load("B.mat")
        self.assertEqual((mat_b.rows, mat_b.cols), (3, 3))
        # Col-major: B.set(r, c, count) with count += 1 in inner r loop
        count = 0.0
        for c in range(3):
            for r in range(3):
                self.assertAlmostEqual(mat_b.get(r, c), count)
                count += 1.0


class TestPdfMakeBasis2(unittest.TestCase):
    def setUp(self):
        self._orig_argv = sys.argv
        self._orig_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)

    def tearDown(self):
        os.chdir(self._orig_cwd)
        sys.argv = self._orig_argv
        self.tmpdir.cleanup()

    def _create_minimal_cache(self, cache_path):
        url_index = (
            "https://bse.pnl.gov/bse/portal/user/anon/js_peid/11535052407933/panel/Main/template/content"
        )
        sample_index_html = (
            'new basisSet("bs_sto3g", "sto-3g", "orbital", "[H]", "published", "F", "F", "2014-01-01")'
        )
        sample_g94 = """H 0
S 3 1.00
      3.42525091           0.15432897
      0.62391373           0.53532814
      0.16885540           0.44463454
****
"""
        cache_data = {
            (url_index, ()): sample_index_html,
            "Gaussian94": {"sto-3g": sample_g94},
            "ProteinDF": {},
            "ProteinDF-DFT": {},
            "ProteinDF-MRD": {},
        }
        with open(cache_path, "wb") as f:
            pickle.dump(cache_data, f)

    def test_pdf_make_basis2_with_cache(self):
        """MakeBasis2 generates ProteinDF basis format from minimal cached data without network."""
        cache_path = os.path.join(self.tmpdir.name, "basis2.cache")
        out_path = os.path.join(self.tmpdir.name, "basis2.out")
        self._create_minimal_cache(cache_path)

        mod = _load_script("pdf-make-basis2.py")
        bs2 = mod.MakeBasis2(cache_path=cache_path, output_path=out_path, uncontract_shells="s")
        bs2.run(basisset_names=["sto-3g"])
        bs2._stdout.close()

        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()

        # Check ProteinDF orbital basis format
        self.assertIn("O-sto-3g.H", content)
        self.assertIn("3.425251e+00  1.543290e-01", content)
        # Check uncontracted shell (MRD) basis format
        self.assertIn("O-sto-3g-s.H", content)
        self.assertIn("3.425251e+00  1.000000e+00", content)

    def test_pdf_make_basis2_show_list(self):
        """MakeBasis2.show_list prints available basis set names to stdout."""
        cache_path = os.path.join(self.tmpdir.name, "basis2.cache")
        self._create_minimal_cache(cache_path)

        mod = _load_script("pdf-make-basis2.py")
        bs2 = mod.MakeBasis2(cache_path=cache_path)
        buf = io.StringIO()
        bs2._stdout = buf
        bs2.show_list()
        self.assertIn("sto-3g", buf.getvalue())

    def test_pdf_make_basis2_missing_pdf_home(self):
        """pdf-make-basis2.py main() requires PDF_HOME environment variable for default cache."""
        mod = _load_script("pdf-make-basis2.py")
        orig_pdf_home = os.environ.pop("PDF_HOME", None)
        try:
            sys.argv = ["pdf-make-basis2.py", "--help"]
            # When PDF_HOME is not set, referencing os.environ['PDF_HOME'] raises KeyError
            with self.assertRaises(KeyError):
                mod.main()
        finally:
            if orig_pdf_home is not None:
                os.environ["PDF_HOME"] = orig_pdf_home

    def test_pdf_make_basis2_unknown_basisset(self):
        """MakeBasis2 logs error and produces no basis definition for unknown basis set."""
        cache_path = os.path.join(self.tmpdir.name, "basis2.cache")
        out_path = os.path.join(self.tmpdir.name, "unknown.out")
        self._create_minimal_cache(cache_path)

        mod = _load_script("pdf-make-basis2.py")
        bs2 = mod.MakeBasis2(cache_path=cache_path, output_path=out_path)
        bs2.run(basisset_names=["unknown-basis-xyz"])
        bs2._stdout.close()

        with open(out_path, "r") as f:
            content = f.read()
        self.assertNotIn("O-unknown-basis-xyz", content)


if __name__ == "__main__":
    unittest.main()
