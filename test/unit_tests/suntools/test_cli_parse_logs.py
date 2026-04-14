#!/usr/bin/env python3
# ---------------------------------------------------------------
# Programmer(s): SUNDIALS contributors
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Test suite for `python -m suntools.cli parse_logs ...`
# ---------------------------------------------------------------

import os
import subprocess
import sys
import unittest


class TestCliParseLogs(unittest.TestCase):
    def _run(self, stdin_text, args):
        env = os.environ.copy()
        # Make sure we import the repo-local package.
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
        env["PYTHONPATH"] = os.path.join(repo_root, "suntools", "src") + os.pathsep + env.get(
            "PYTHONPATH", ""
        )
        proc = subprocess.run(
            [sys.executable, "-m", "suntools.cli", "parse_logs", *args],
            input=stdin_text,
            text=True,
            capture_output=True,
            env=env,
            check=False,
        )
        return proc

    def test_filters_integrator_only(self):
        log = (
            "[INFO][rank 0][CVode][begin-step-attempt] step = 1, tn = 0.0, h = 1.0\n"
            "[INFO][rank 0][CVode][begin-nonlinear-solve] tol = 1e-4\n"
            "[INFO][rank 0][CVode][end-nonlinear-solve] status = success, iters = 3\n"
            "[INFO][rank 0][CVode][end-step-attempt] status = success\n"
        )
        proc = self._run(log, ['--filter=integrator'])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("begin-step-attempt", proc.stdout)
        self.assertIn("end-step-attempt", proc.stdout)
        self.assertNotIn("begin-nonlinear-solve", proc.stdout)
        self.assertNotIn("end-nonlinear-solve", proc.stdout)

    def test_filters_nonlinear_and_linear(self):
        log = (
            "[INFO][rank 0][CVode][begin-linear-solve] iterative = 1\n"
            "[INFO][rank 0][CVode][end-linear-solve] status = success, iters = 5\n"
            "[INFO][rank 0][CVode][begin-nonlinear-solve] tol = 1e-4\n"
            "[INFO][rank 0][CVode][end-nonlinear-solve] status = success, iters = 3\n"
        )
        proc = self._run(log, ['--filter=nonlinear,linear'])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("begin-linear-solve", proc.stdout)
        self.assertIn("end-linear-solve", proc.stdout)
        self.assertIn("begin-nonlinear-solve", proc.stdout)
        self.assertIn("end-nonlinear-solve", proc.stdout)

    def test_preserves_array_continuations(self):
        log = (
            "[DEBUG][rank 0][Scope][begin-step-attempt] step = 1, tn = 0.0, h = 1.0\n"
            "[DEBUG][rank 0][Scope][label] u_0(:) =\n"
            " 1\n"
            " 2\n"
            "\n"
            "[DEBUG][rank 0][Scope][end-step-attempt] status = success\n"
        )
        proc = self._run(log, ['--filter=integrator'])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("u_0(:) =", proc.stdout)
        self.assertIn("\n 1\n", proc.stdout)
        self.assertIn("\n 2\n", proc.stdout)

    def test_nonlinear_excludes_nested_linear_regions(self):
        log = (
            "[INFO][rank 0][CVode][begin-nonlinear-solve] tol = 1e-4\n"
            "[INFO][rank 0][CVode][newton-iter] iter = 0\n"
            "[INFO][rank 0][CVode][begin-linear-solve] iterative = 1\n"
            "[INFO][rank 0][CVode][some-linear-detail] x = 1\n"
            "[INFO][rank 0][CVode][end-linear-solve] status = success, iters = 2\n"
            "[INFO][rank 0][CVode][newton-iter] iter = 1\n"
            "[INFO][rank 0][CVode][end-nonlinear-solve] status = success, iters = 2\n"
        )
        proc = self._run(log, ['--filter=nonlinear'])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("begin-nonlinear-solve", proc.stdout)
        self.assertIn("end-nonlinear-solve", proc.stdout)
        self.assertIn("newton-iter", proc.stdout)
        self.assertNotIn("begin-linear-solve", proc.stdout)
        self.assertNotIn("end-linear-solve", proc.stdout)
        self.assertNotIn("some-linear-detail", proc.stdout)

    def test_nonlinear_includes_nonlog_lines_in_region(self):
        log = (
            "[INFO][rank 0][CVode][begin-nonlinear-solve] tol = 1e-4\n"
            "this is not a logger line but is inside nonlinear\n"
            "[INFO][rank 0][CVode][begin-linear-solve] iterative = 1\n"
            "this is not a logger line but is inside linear\n"
            "[INFO][rank 0][CVode][end-linear-solve] status = success, iters = 2\n"
            "[INFO][rank 0][CVode][end-nonlinear-solve] status = success, iters = 2\n"
        )
        proc = self._run(log, ["--filter=nonlinear"])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("inside nonlinear", proc.stdout)
        self.assertNotIn("inside linear", proc.stdout)


def run_tests():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestCliParseLogs))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
