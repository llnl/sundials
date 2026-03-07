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
# Test suite for _parse_logfile_line function
# ---------------------------------------------------------------

import unittest
import sys

# Import the logs module
import logs


class TestParseLogfileLine(unittest.TestCase):
    """Test the _parse_logfile_line function."""

    def test_parse_basic_log_line(self):
        """Test parsing a basic log line with all components."""
        line = "[INFO][rank 0][TestScope][test-label] key1 = value1, key2 = value2"
        result = logs._parse_logfile_line(line, 0, [])

        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 5)
        expected_keys = {"loglvl", "rank", "scope", "label", "payload"}
        self.assertEqual(set(result.keys()), expected_keys)
        self.assertIsInstance(result["loglvl"], str)
        self.assertEqual(result["loglvl"], "INFO")
        self.assertIsInstance(result["rank"], int)
        self.assertEqual(result["rank"], 0)
        self.assertIsInstance(result["scope"], str)
        self.assertEqual(result["scope"], "TestScope")
        self.assertIsInstance(result["label"], str)
        self.assertEqual(result["label"], "test-label")
        self.assertIsInstance(result["payload"], dict)
        self.assertEqual(len(result["payload"]), 2)
        self.assertEqual(result["payload"]["key1"], "value1")
        self.assertEqual(result["payload"]["key2"], "value2")

    def test_parse_empty_payload(self):
        """Test parsing line with empty payload."""
        line = "[INFO][rank 0][TestScope][label]"
        result = logs._parse_logfile_line(line, 0, [])

        self.assertIn("payload", result)
        self.assertIsInstance(result["payload"], dict)
        self.assertEqual(len(result["payload"]), 0)

    def test_parse_whitespace_only_payload(self):
        """Test parsing line with whitespace-only payload."""
        line = "[INFO][rank 0][TestScope][label]    "
        result = logs._parse_logfile_line(line, 0, [])

        self.assertIn("payload", result)
        self.assertIsInstance(result["payload"], dict)
        # Whitespace should be trimmed, resulting in empty payload
        self.assertEqual(len(result["payload"]), 0)

    def test_parse_payload_with_array(self):
        """Test parsing line with array in payload."""
        all_lines = [
            "[INFO][rank 0][TestScope][label] array(:) =",
            " 1.0",
            " 2.0",
            " 3.0",
            "[INFO][rank 0][TestScope][next-label]",
        ]

        result = logs._parse_logfile_line(all_lines[0], 0, all_lines)

        self.assertIn("payload", result)
        self.assertIn("array(:)", result["payload"])
        self.assertIsInstance(result["payload"]["array(:)"], list)
        self.assertEqual(len(result["payload"]["array(:)"]), 3)

    def test_parse_line_with_spaces_in_scope_and_label(self):
        """Test parsing line with spaces in scope and label."""
        line = "[INFO][rank 0][Test Scope][test label here] key = val"
        result = logs._parse_logfile_line(line, 0, [])

        # The regex (.*) for scope should capture anything, including spaces
        self.assertEqual(result["scope"], "Test Scope")
        self.assertEqual(result["label"], "test label here")

    def test_parse_line_with_nested_brackets_in_payload(self):
        """Test parsing line with nested brackets in the payload."""
        line = "[INFO][rank 0][TestScope][TestLabel] data[i] = 42"
        result = logs._parse_logfile_line(line, 0, [])

        self.assertIn("data[i]", result["payload"])
        self.assertEqual(result["payload"]["data[i]"], 42)

    def test_empty_result_on_no_match(self):
        """Test that empty dict is returned when line doesn't match pattern."""
        bad_lines = [
            "",  # empty line
            "  ",  # only white space
            "Not a log line at all",
            "# This is a comment",
            "[INFO] missing rank",
            "[INFO][0] missing rank keyword",
            "[INFO][rank 0] missing scope and label",
            "[INFO][rank 0][Scope] missing label",
            "INFO][rank 0][Scope][Label] missing opening bracket"
            "[INFO DEBUG][rank 0][Scope][Label] bad log type",
        ]

        for line in bad_lines:
            result = logs._parse_logfile_line(line, 0, [])
            self.assertIsInstance(result, dict)
            self.assertEqual(len(result), 0, f"Failed for line: {line}")


def run_tests():
    """Run all tests and print summary."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestParseLogfileLine))

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

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
