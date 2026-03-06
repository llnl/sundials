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
# Test suite for _parse_logfile_payload function
# ---------------------------------------------------------------

import unittest
import sys

# Import the logs module
import logs


class TestParseLogfilePayload(unittest.TestCase):
    """Test the _parse_logfile_payload function."""

    def test_simple_payload(self):
        """Test parsing simple payload of key=value pairs."""
        payload = "  int key = 1,   float_key=0.123, sci-key=1.23e+2,strkey = foo, equation = x=y+z,empty_key  "
        result = logs._parse_logfile_payload(payload, 0, [])

        self.assertEqual(len(result), 6)
        self.assertIsInstance(result["int key"], int)
        self.assertEqual(result["int key"], 1)
        self.assertIsInstance(result["float_key"], float)
        self.assertAlmostEqual(result["float_key"], 0.123)
        self.assertIsInstance(result["sci-key"], float)
        self.assertAlmostEqual(result["sci-key"], 1.23e2)
        self.assertIsInstance(result["strkey"], str)
        self.assertEqual(result["strkey"], "foo")
        self.assertIsInstance(result["equation"], str)
        self.assertEqual(result["equation"], "x=y+z")
        self.assertEqual(result["empty_key"], "")

    def test_empty_payload(self):
        """Test parsing empty payload."""
        payload = ""
        result = logs._parse_logfile_payload(payload, 0, [])

        self.assertEqual(len(result), 0)

    def test_whitespace_only_payload(self):
        """Test parsing whitespace-only payload."""
        payload = "   "
        result = logs._parse_logfile_payload(payload, 0, [])

        self.assertEqual(len(result), 0)

    def test_array_indicator_default(self):
        """Test parsing array with default (:) indicator."""
        all_lines = [
            "array(:) =",
            " 1",
            " 0.123",
            " 1.23e+2",
            " foo",
            "[INFO][rank 0][Scope][next-label] key = val",
            " 3.0",  # This should not be included
        ]
        payload = all_lines[0]
        result = logs._parse_logfile_payload(payload, 0, all_lines)

        self.assertIn("array(:)", result)
        self.assertIsInstance(result["array(:)"], list)
        # Should only have 4 values, stopping at the next log line
        self.assertEqual(len(result["array(:)"]), 4)
        self.assertIsInstance(result["array(:)"][0], int)
        self.assertEqual(result["array(:)"][0], 1)
        self.assertIsInstance(result["array(:)"][1], float)
        self.assertAlmostEqual(result["array(:)"][1], 0.123)
        self.assertIsInstance(result["array(:)"][2], float)
        self.assertAlmostEqual(result["array(:)"][2], 1.23e2)
        self.assertIsInstance(result["array(:)"][3], str)
        self.assertAlmostEqual(result["array(:)"][3], "foo")

    def test_array_stops_at_empty_line(self):
        """Test that array parsing stops at empty line."""
        all_lines = [
            "[INFO][rank 0][Scope][label] array(:) =",
            " 1.0",
            " 2.0",
            "",  # Empty line
            " 3.0",  # This should not be included
        ]
        payload = "array(:) ="
        result = logs._parse_logfile_payload(payload, 0, all_lines)

        # Should only have 2 values, stopping at the empty line
        self.assertEqual(len(result["array(:)"]), 2)

    def test_empty_array(self):
        """Test parsing array with no values."""
        all_lines = [
            "[INFO][rank 0][Scope][label] empty(:) =",
            "[INFO][rank 0][Scope][next-label]",
        ]
        payload = "empty(:) ="
        result = logs._parse_logfile_payload(payload, 0, all_lines)

        # Empty array should still be created
        self.assertIn("empty(:)", result)
        self.assertIsInstance(result["empty(:)"], list)
        self.assertEqual(len(result["empty(:)"]), 0)

    def test_custom_array_indicator(self):
        """Test parsing with custom array indicator."""
        all_lines = [
            "[INFO][rank 0][Scope][label] data[*] =",
            " 1.0",
            " 2.0",
            "[INFO][rank 0][Scope][next]",
        ]
        payload = "data[*] ="
        result = logs._parse_logfile_payload(payload, 0, all_lines, array_indicator="[*]")

        self.assertIn("data[*]", result)
        self.assertIsInstance(result["data[*]"], list)
        self.assertEqual(len(result["data[*]"]), 2)


def run_tests():
    """Run all tests and print summary."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestParseLogfilePayload))

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
