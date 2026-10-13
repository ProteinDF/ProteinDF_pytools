#!/usr/bin/env python
# -*- coding: utf-8 -*-

import contextlib
import io
import os
import tempfile
import unittest

from proteindf_tools.process import Process


class TestProcess(unittest.TestCase):
    def test_cmd_echo(self):
        p = Process()
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            ret = p.cmd("echo hello").commit()
        self.assertEqual(ret, 0)
        self.assertIn("hello", stdout_buf.getvalue())

    def test_cmd_true(self):
        p = Process()
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            ret = p.cmd("true").commit()
        self.assertEqual(ret, 0)

    def test_cmd_failure(self):
        p = Process()
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            ret = p.cmd("false").commit()
        self.assertNotEqual(ret, 0)

    def test_pipe(self):
        p = Process()
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            ret = p.cmd('echo "foo bar"').pipe("grep foo").commit()
        self.assertEqual(ret, 0)
        self.assertIn("foo bar", stdout_buf.getvalue())

    def test_pipe_failure(self):
        p = Process()
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            ret = p.cmd('echo "foo bar"').pipe("grep nonexistent").commit()
        self.assertNotEqual(ret, 0)

    def test_invalid_command(self):
        p = Process()
        stderr_buf = io.StringIO()
        with contextlib.redirect_stderr(stderr_buf):
            with self.assertRaises(OSError):
                p.cmd("nonexistent_command_12345")

    def test_is_shell_script(self):
        p = Process()
        self.assertTrue(p._is_shell_script("run.sh arg1"))
        self.assertFalse(p._is_shell_script("echo hello"))
        self.assertFalse(p._is_shell_script(""))

    def test_commit_with_logfile(self):
        """logfile_path 指定時、Process.commit はログファイルに出力を書き込んで終了するべき。"""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            log_path = f.name
        try:
            p = Process(logfile_path=log_path)
            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()
            with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                ret = p.cmd("echo logged_output").commit()
            self.assertEqual(ret, 0)
            with open(log_path, "r", encoding="utf-8") as f:
                content = f.read()
            self.assertIn("logged_output", content)

            p2 = Process(logfile_path=log_path)
            with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                ret2 = p2.cmd("sh -c 'echo logged_err 1>&2'").commit()
            self.assertEqual(ret2, 0)
            with open(log_path, "r", encoding="utf-8") as f:
                content = f.read()
            self.assertIn("logged_output", content)
            self.assertIn("logged_err", content)
        finally:
            if os.path.exists(log_path):
                os.remove(log_path)


if __name__ == "__main__":
    unittest.main()
