#!/usr/bin/env python3
# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# Test suite for SUNDIALS log parser
# ---------------------------------------------------------------

import unittest
import tempfile
import os
import sys

# Import the logs module
import logs


class TestGetHistory(unittest.TestCase):
    """Test the get_history function for extracting time series data."""

    def setUp(self):
        """Create a multi-step log for history extraction."""
        self.test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        self.test_log.write(
            # ---- Step 1 ----
            "[INFO][rank 0][TestScope][begin-step-attempt] step = 1, tn = 0.0, h = 0.1\n"
            # Stages 1 - explicit
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 1\n"
            "[INFO][rank 0][TestScope][end-stages-list] status = stage success\n"
            # Stage 2 - implicit
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 2\n"
            # Nonlinear solve
            "[INFO][rank 0][TestScope][begin-nonlinear-solve] solver = Newton\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            # Iteration 1
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            # Linear solve
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            # Iteration 2
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 2\n"
            # Linear solve
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            # End nonlinear solve
            "[INFO][rank 0][TestScope][end-nonlinear-solve] status = nls success\n"
            "[INFO][rank 0][TestScope][end-stages-list]\n"
            # Stage 3
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 3, status = stage success\n"
            "[INFO][rank 0][TestScope][begin-nonlinear-solve] solver = Newton\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 2\n"
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-nonlinear-solve] status = nls success\n"
            "[INFO][rank 0][TestScope][end-stages-list]\n"
            "[INFO][rank 0][TestScope][end-step-attempt] status = step success\n"
            # ---- Step 2 ----
            "[INFO][rank 0][TestScope][begin-step-attempt] step = 2, tn = 0.1, h = 0.2\n"
            # Stages 1 - explicit
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 1\n"
            "[INFO][rank 0][TestScope][end-stages-list] status = stage success\n"
            # Stage 2 - implicit
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 2\n"
            # Nonlinear solve
            "[INFO][rank 0][TestScope][begin-nonlinear-solve] solver = Newton\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            # Iteration 1
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            # Linear solve
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            # Iteration 2
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 2\n"
            # Linear solve
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            # End nonlinear solve
            "[INFO][rank 0][TestScope][end-nonlinear-solve] status = nls success\n"
            "[INFO][rank 0][TestScope][end-stages-list]\n"
            # Stage 3
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 3, status = stage success\n"
            "[INFO][rank 0][TestScope][begin-nonlinear-solve] solver = Newton\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 2\n"
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-nonlinear-solve] status = nls success\n"
            "[INFO][rank 0][TestScope][end-stages-list]\n"
            "[INFO][rank 0][TestScope][end-step-attempt] status = step success\n"
            # ---- Step 2, Attempt 2 ----
            "[INFO][rank 0][TestScope][begin-step-attempt] step = 2, tn = 0.1, h = 0.15\n"
            # Stages 1 - explicit
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 1\n"
            "[INFO][rank 0][TestScope][end-stages-list] status = stage success\n"
            # Stage 2 - implicit
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 2\n"
            # Nonlinear solve
            "[INFO][rank 0][TestScope][begin-nonlinear-solve] solver = Newton\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            # Iteration 1
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            # Linear solve
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            # Iteration 2
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 2\n"
            # Linear solve
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            # End nonlinear solve
            "[INFO][rank 0][TestScope][end-nonlinear-solve] status = nls success\n"
            "[INFO][rank 0][TestScope][end-stages-list]\n"
            # Stage 3
            "[INFO][rank 0][TestScope][begin-stages-list]\n"
            "[INFO][rank 0][TestScope][stage] stage = 3, status = stage success\n"
            "[INFO][rank 0][TestScope][begin-nonlinear-solve] solver = Newton\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 2\n"
            "[INFO][rank 0][TestScope][begin-linear-solve] solver = SPGMR\n"
            "[INFO][rank 0][TestScope][begin-iterations-list]\n"
            "[INFO][rank 0][TestScope][iteration] iteration = 1\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-linear-solve] status = ls success\n"
            "[INFO][rank 0][TestScope][end-iterations-list]\n"
            "[INFO][rank 0][TestScope][end-nonlinear-solve] status = nls success\n"
            "[INFO][rank 0][TestScope][end-stages-list]\n"
            "[INFO][rank 0][TestScope][end-step-attempt] status = step success\n"
        )
        self.test_log.close()
        self.data = logs.log_file_to_list(self.test_log.name)
        logs.print_log(self.data)

    def tearDown(self):
        """Clean up test files."""
        try:
            os.unlink(self.test_log.name)
        except:
            pass

    def test_extract_nested_history(self):
        """Test extracting history of a key."""

        result = logs.extract(self.data, "step")
        print()
        print(70*"=")
        print("result = ", result)
        print(70*"=")

        # result = logs.extract(self.data, ["foo"])
        # print()
        # print(f"result is dict {isinstance(result, dict)}")
        # print("result = ", result)

        # result = logs.extract(self.data, ["stages.stage"])
        # print()
        # print(f"result is dict {isinstance(result, dict)}")
        # print("result = ", result)

        #result = logs.extract(self.data, ["step", "tn", "h"])
        #result = logs.extract(self.data, "stages.stage")
        #result = logs.extract(self.data, ["step", "stages.stage.nonlinear-solve"]) # Should raise error, stage is a key not a dict or list, but does not because it breaks the next case e.g., explicit first stage
        result = logs.extract(self.data, ["step", "stages.stage", "stages.nonlinear-solve.status"]) # RIGHT WAY
        #result = logs.extract(self.data, ["stages.nonlinear-solve.iterations.iteration", "stages.nonlinear-solve.iterations.linear-solve.status"]) # RIGHT WAY

        print()
        print(70*"=")
        print("result = ", result)
        print(70*"=")

        # result = logs.extract(self.data, ["step", "stage"]) # WOULD like this to work will be complicated

        # result = logs.extract2(self.data, "step")

        # print()
        # print(f"result is dict {isinstance(result, dict)}")
        # print("result = ", result)

        # Would like to do subsubkey1 but this would cause problems for sub step sizes
        # result = logs.extract(self.data, ["step", "subsection.subkey1", "subsection.subsubsection.subsubkey1"])

        # print()
        # print(f"result is dict {isinstance(result, dict)}")
        # print("result = ", result)

    # To only pass stages -- need to look recursively for that key
    # checking dicts and lists within the current dict (or list of dicts)
    # THIS IS COMPLICATED AND PROBABLY A BAD IDEA -- take stage need to find the stages
    # list then loop over all the dicts in that list to get all the stage values

# def key_exists():





    # def test_extract_nested_2_history(self):
    #     """Test extracting history of a key."""
    #     steps, times, values = logs.get_history(self.data, "subsubsection.subsubkey2")

    #     # Should return three lists
    #     self.assertEqual(len(steps), len(times))
    #     self.assertEqual(len(times), len(values))

    #     # Should have all steps
    #     self.assertEqual(len(steps), 1)

    #     print()
    #     print(f"steps  = {steps}")
    #     print(f"times  = {times}")
    #     print(f"values = {values}")


    # def test_extract_basic_history(self):
    #     """Test extracting history of a key."""
    #     steps, times, values = logs.get_history(self.data, "h")

    #     # Should return three lists
    #     self.assertEqual(len(steps), len(times))
    #     self.assertEqual(len(times), len(values))

    #     # Should have all steps
    #     self.assertEqual(len(steps), 4)

    # def test_filter_by_status(self):
    #     """Test filtering history by step status."""
    #     steps, times, values = logs.get_history(self.data, "metric", step_status="success")

    #     # Should exclude the failed step
    #     self.assertEqual(len(steps), 3)

    #     # All included steps should have 'success' status
    #     for step_data in self.data:
    #         if "failed" in step_data.get("status", ""):
    #             self.assertNotIn(step_data["metric"], values)

    # def test_filter_by_time_range(self):
    #     """Test filtering history by time range."""
    #     steps, times, values = logs.get_history(self.data, "h", time_range=[0.01, 0.03])

    #     # All times should be within range
    #     for t in times:
    #         self.assertGreaterEqual(t, 0.01)
    #         self.assertLessEqual(t, 0.03)

    # def test_filter_by_step_range(self):
    #     """Test filtering history by step number range."""
    #     steps, times, values = logs.get_history(self.data, "h", step_range=[2, 3])

    #     # All steps should be within range
    #     for s in steps:
    #         self.assertGreaterEqual(s, 2)
    #         self.assertLessEqual(s, 3)

    # def test_nonexistent_key(self):
    #     """Test extracting history for a nonexistent key."""
    #     steps, times, values = logs.get_history(self.data, "nonexistent_key")

    #     # Should return empty lists
    #     self.assertEqual(len(steps), 0)
    #     self.assertEqual(len(times), 0)
    #     self.assertEqual(len(values), 0)

    # def test_group_by_level(self):
    #     """Test grouping history by time level."""
    #     steps_by_level, times_by_level, values_by_level = logs.get_history(
    #         self.data, "h", group_by_level=True
    #     )

    #     # Should return dictionaries keyed by level
    #     self.assertIsInstance(steps_by_level, dict)
    #     self.assertIsInstance(times_by_level, dict)
    #     self.assertIsInstance(values_by_level, dict)

    #     # Should have level 0 (top level)
    #     self.assertIn(0, steps_by_level)


def run_tests():
    """Run all tests and print summary."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestGetHistory))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print("=" * 70)

    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
