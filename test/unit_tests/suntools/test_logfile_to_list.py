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
import json
import os
import sys

# Import the logs module
import logs


class TestLogParsing(unittest.TestCase):
    """Test parsing of individual log line components."""

    def test_simple_logfile(self):
        """Test that a simple log file of key=value pairs is parsed correctly."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            # Step 0 - only begin and end labels
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key2 = val2\n"
            # Step 1 - key-value pairs within the step
            "[INFO][rank 0][TestScope][begin-step-attempt] key3 = val3\n"
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][label] key5 = val5\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key6 = val6\n"
            # Step 2 - empty payloads
            "[INFO][rank 0][TestScope][begin-step-attempt]\n"
            "[INFO][rank 0][TestScope][label]\n"
            "[INFO][rank 0][TestScope][end-step-attempt]\n"
        )
        test_log.close()

        try:
            data = logs.log_file_to_list(test_log.name)
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 3)

            data0 = data[0]
            self.assertIsInstance(data0, dict)
            self.assertEqual(len(data0), 3)  # parser adds time level to dict

            self.assertIn("level", data0)
            self.assertEqual(data0["level"], 0)

            for i in range(1, 3):
                self.assertIn(f"key{i}", data0)
                self.assertEqual(data0[f"key{i}"], f"val{i}")

            data1 = data[1]
            self.assertIsInstance(data1, dict)
            self.assertEqual(len(data1), 5)  # parser adds 1 key (time level)

            self.assertIn("level", data1)
            self.assertEqual(data1["level"], 0)

            for i in range(3, 7):
                self.assertIn(f"key{i}", data1)
                self.assertEqual(data1[f"key{i}"], f"val{i}")

            data2 = data[2]
            self.assertIsInstance(data2, dict)
            self.assertEqual(len(data2), 1)  # parser adds 1 key (time level)

            self.assertIn("level", data2)
            self.assertEqual(data2["level"], 0)

        finally:
            os.unlink(test_log.name)

    def test_nested_dictionary_structure(self):
        """Test that begin-X / end-X pairs create nested dictionaries."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            # Step 0 - nested dictionary
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-subsection] subkey1 = subval1\n"
            "[INFO][rank 0][TestScope][label] subkey2 = subval2\n"
            "[INFO][rank 0][TestScope][end-subsection] subkey3 = subval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key4 = val4\n"
            # Step 1 - nested dictionary with multi string name
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-subsection-long-name] subkey1 = subval1\n"
            "[INFO][rank 0][TestScope][label] subkey2 = subval2\n"
            "[INFO][rank 0][TestScope][end-subsection-long-name] subkey3 = subval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key4 = val4\n"
            # Step 2 - more than on nested dictionary
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-subsection-1] subkey1 = subval1\n"
            "[INFO][rank 0][TestScope][label] subkey2 = subval2\n"
            "[INFO][rank 0][TestScope][end-subsection-1] subkey3 = subval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][begin-subsection-2] subkey4 = subval4\n"
            "[INFO][rank 0][TestScope][label] subkey5 = subval5\n"
            "[INFO][rank 0][TestScope][end-subsection-2] subkey6 = subval6\n"
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key5 = val5\n"
            # Step 3 - Reopen and update dictionary
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-subsection] subkey1 = subval1\n"
            "[INFO][rank 0][TestScope][label] subkey2 = subval2\n"
            "[INFO][rank 0][TestScope][end-subsection] subkey3 = subval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][begin-subsection] subkey4 = subval4\n"
            "[INFO][rank 0][TestScope][label] subkey5 = subval5\n"
            "[INFO][rank 0][TestScope][end-subsection] subkey6 = subval6\n"
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key5 = val5\n"
            # Step 4 - nested dictionaries with nested dictionaries
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-subsection] subkey1 = subval1\n"
            "[INFO][rank 0][TestScope][label] subkey2 = subval2\n"
            "[INFO][rank 0][TestScope][begin-subsubsection] subsubkey1 = subsubval1\n"
            "[INFO][rank 0][TestScope][label] subsubkey2 = subsubval2\n"
            "[INFO][rank 0][TestScope][begin-subsubsubsection] subsubsubkey1 = subsubsubval1\n"
            "[INFO][rank 0][TestScope][label] subsubsubkey2 = subsubsubval2\n"
            "[INFO][rank 0][TestScope][end-subsubsubsection] subsubsubkey3 = subsubsubval3\n"
            "[INFO][rank 0][TestScope][label] subsubkey3 = subsubval3\n"
            "[INFO][rank 0][TestScope][end-subsubsection] subsubkey4 = subsubval4\n"
            "[INFO][rank 0][TestScope][label] subkey3 = subval3\n"
            "[INFO][rank 0][TestScope][end-subsection] subkey4 = subval4\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key4 = val4\n"
        )
        test_log.close()

        try:
            data = logs.log_file_to_list(test_log.name)
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 5)

            # Step 0 - nested dictionary
            data0 = data[0]
            self.assertIsInstance(data0, dict)
            self.assertEqual(len(data0), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data0)
            self.assertEqual(data0["level"], 0)

            for i in range(1, 5):
                self.assertIn(f"key{i}", data0)
                self.assertEqual(data0[f"key{i}"], f"val{i}")

            self.assertIn("subsection", data0)
            data0_subsection = data0["subsection"]
            self.assertIsInstance(data0_subsection, dict)
            self.assertEqual(len(data0_subsection), 3)

            for i in range(1, 4):
                self.assertIn(f"subkey{i}", data0_subsection)
                self.assertEqual(data0_subsection[f"subkey{i}"], f"subval{i}")

            # Step 1 - nested dictionary with multi string name
            data1 = data[1]
            self.assertIsInstance(data1, dict)
            self.assertEqual(len(data1), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data1)
            self.assertEqual(data1["level"], 0)

            for i in range(1, 5):
                self.assertIn(f"key{i}", data1)
                self.assertEqual(data1[f"key{i}"], f"val{i}")

            self.assertIn("subsection-long-name", data1)
            data1_subsection = data1["subsection-long-name"]
            self.assertIsInstance(data1_subsection, dict)
            self.assertEqual(len(data1_subsection), 3)

            for i in range(1, 4):
                self.assertIn(f"subkey{i}", data1_subsection)
                self.assertEqual(data1_subsection[f"subkey{i}"], f"subval{i}")

            # Step 2 - more than one nested dictionary
            data2 = data[2]
            self.assertIsInstance(data2, dict)
            self.assertEqual(len(data2), 8)  # parser adds 1 key (time level)

            self.assertIn("level", data2)
            self.assertEqual(data2["level"], 0)

            for i in range(1, 6):
                self.assertIn(f"key{i}", data2)
                self.assertEqual(data2[f"key{i}"], f"val{i}")

            self.assertIn("subsection-1", data2)
            data2_subsection1 = data2["subsection-1"]
            self.assertIsInstance(data2_subsection1, dict)
            self.assertEqual(len(data2_subsection1), 3)

            for i in range(1, 4):
                self.assertIn(f"subkey{i}", data2_subsection1)
                self.assertEqual(data2_subsection1[f"subkey{i}"], f"subval{i}")

            self.assertIn("subsection-2", data2)
            data2_subsection2 = data2["subsection-2"]
            self.assertIsInstance(data2_subsection2, dict)
            self.assertEqual(len(data2_subsection2), 3)

            for i in range(4, 7):
                self.assertIn(f"subkey{i}", data2_subsection2)
                self.assertEqual(data2_subsection2[f"subkey{i}"], f"subval{i}")

            # Step 3 - Reopen and update dictionary
            data3 = data[3]
            self.assertIsInstance(data3, dict)
            self.assertEqual(len(data3), 7)  # parser adds 1 key (time level)

            self.assertIn("level", data3)
            self.assertEqual(data3["level"], 0)

            for i in range(1, 6):
                self.assertIn(f"key{i}", data3)
                self.assertEqual(data3[f"key{i}"], f"val{i}")

            self.assertIn("subsection", data3)
            data3_subsection = data3["subsection"]
            self.assertIsInstance(data3_subsection, dict)
            self.assertEqual(len(data3_subsection), 6)

            for i in range(1, 7):
                self.assertIn(f"subkey{i}", data3_subsection)
                self.assertEqual(data3_subsection[f"subkey{i}"], f"subval{i}")

            # Step 4 - nested dictionaries with nested dictionaries
            data4 = data[4]
            self.assertIsInstance(data4, dict)
            self.assertEqual(len(data4), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data4)
            self.assertEqual(data4["level"], 0)

            for i in range(1, 5):
                self.assertIn(f"key{i}", data4)
                self.assertEqual(data4[f"key{i}"], f"val{i}")

            self.assertIn("subsection", data4)
            data4_subsection = data4["subsection"]
            self.assertIsInstance(data4_subsection, dict)
            self.assertEqual(len(data4_subsection), 5)

            for i in range(1, 5):
                self.assertIn(f"subkey{i}", data4_subsection)
                self.assertEqual(data4_subsection[f"subkey{i}"], f"subval{i}")

            self.assertIn("subsubsection", data4_subsection)
            data4_subsubsection = data4_subsection["subsubsection"]
            self.assertIsInstance(data4_subsubsection, dict)
            self.assertEqual(len(data4_subsubsection), 5)

            for i in range(1, 5):
                self.assertIn(f"subsubkey{i}", data4_subsubsection)
                self.assertEqual(data4_subsubsection[f"subsubkey{i}"], f"subsubval{i}")

            self.assertIn("subsubsubsection", data4_subsubsection)
            data4_subsubsubsection = data4_subsubsection["subsubsubsection"]
            self.assertIsInstance(data4_subsubsubsection, dict)
            self.assertEqual(len(data4_subsubsubsection), 3)

            for i in range(1, 4):
                self.assertIn(f"subsubsubkey{i}", data4_subsubsubsection)
                self.assertEqual(data4_subsubsubsection[f"subsubsubkey{i}"], f"subsubsubval{i}")

        finally:
            os.unlink(test_log.name)

    def test_list_structure(self):
        """Test that begin-X-list / end-X-list pairs create lists."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            # Step 0 - add list with one entry
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-items-list] listkey1 = listval1\n"
            "[INFO][rank 0][TestScope][label] listkey2 = listval2\n"
            "[INFO][rank 0][TestScope][end-items-list] listkey3 = listval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key4 = val4\n"
            # Step 1 - add list, reopen and add additional entry
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-items-list] listkey1 = listval1\n"
            "[INFO][rank 0][TestScope][label] listkey2 = listval2\n"
            "[INFO][rank 0][TestScope][end-items-list] listkey3 = listval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][begin-items-list] listkey4 = listval4\n"
            "[INFO][rank 0][TestScope][label] listkey5 = listval5\n"
            "[INFO][rank 0][TestScope][end-items-list] listkey6 = listval6\n"
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key5 = val5\n"
            # Step 2 - add a second list
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-items-first-list] listkey1 = listval1\n"
            "[INFO][rank 0][TestScope][label] listkey2 = listval2\n"
            "[INFO][rank 0][TestScope][end-items-first-list] listkey3 = listval3\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][begin-items-second-list] listkey4 = listval4\n"
            "[INFO][rank 0][TestScope][label] listkey5 = listval5\n"
            "[INFO][rank 0][TestScope][end-items-second-list] listkey6 = listval6\n"
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key5 = val5\n"
            # Step 3 - add a list within a list
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-items-outer-list] outerlistkey1 = outerlistval1\n"
            "[INFO][rank 0][TestScope][label] outerlistkey2 = outerlistval2\n"
            "[INFO][rank 0][TestScope][begin-items-inner-list] innerlistkey1 = innerlistval1\n"
            "[INFO][rank 0][TestScope][label] innerlistkey2 = innerlistval2\n"
            "[INFO][rank 0][TestScope][end-items-inner-list] innerlistkey3 = innerlistval3\n"
            "[INFO][rank 0][TestScope][label] outerlistkey3 = outerlistval3\n"
            "[INFO][rank 0][TestScope][end-items-outer-list] outerlistkey4 = outerlistval4\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key4 = val4\n"
        )
        test_log.close()

        try:
            data = logs.log_file_to_list(test_log.name)
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 4)

            # Step 0 - add list
            data0 = data[0]
            self.assertIsInstance(data0, dict)
            self.assertEqual(len(data0), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data0)
            self.assertEqual(data0["level"], 0)

            for i in range(1, 4):
                self.assertIn(f"key{i}", data0)
                self.assertEqual(data0[f"key{i}"], f"val{i}")

            self.assertIn("items", data0)
            data0_items = data0["items"]
            self.assertIsInstance(data0_items, list)
            self.assertEqual(len(data0_items), 1)  # opening a list appends a new dictionary

            data0_items0 = data0_items[0]
            self.assertIsInstance(data0_items0, dict)
            self.assertEqual(len(data0_items0), 3)

            for i in range(1, 4):
                self.assertIn(f"listkey{i}", data0_items0)
                self.assertEqual(data0_items0[f"listkey{i}"], f"listval{i}")

            # Step 1 - add list, reopen and add additional entry
            data1 = data[1]
            self.assertIsInstance(data1, dict)
            self.assertEqual(len(data1), 7)  # parser adds 1 key (time level)

            self.assertIn("level", data1)
            self.assertEqual(data1["level"], 0)

            for i in range(1, 4):
                self.assertIn(f"key{i}", data1)
                self.assertEqual(data1[f"key{i}"], f"val{i}")

            self.assertIn("items", data1)
            data1_items = data1["items"]
            self.assertIsInstance(data1_items, list)
            self.assertEqual(len(data1_items), 2)  # opening a list appends a new dictionary

            data1_items0 = data1_items[0]
            self.assertIsInstance(data1_items0, dict)
            self.assertEqual(len(data1_items0), 3)

            for i in range(1, 4):
                self.assertIn(f"listkey{i}", data1_items0)
                self.assertEqual(data1_items0[f"listkey{i}"], f"listval{i}")

            data1_items1 = data1_items[1]
            self.assertIsInstance(data1_items1, dict)
            self.assertEqual(len(data1_items1), 3)

            for i in range(4, 7):
                self.assertIn(f"listkey{i}", data1_items1)
                self.assertEqual(data1_items1[f"listkey{i}"], f"listval{i}")

            # Step 2 - add a second list
            data2 = data[2]
            self.assertIsInstance(data2, dict)
            self.assertEqual(len(data2), 8)  # parser adds 1 key (time level)

            self.assertIn("level", data2)
            self.assertEqual(data2["level"], 0)

            for i in range(1, 4):
                self.assertIn(f"key{i}", data2)
                self.assertEqual(data2[f"key{i}"], f"val{i}")

            self.assertIn("items-first", data2)
            data2_items_first = data2["items-first"]
            self.assertIsInstance(data2_items_first, list)
            self.assertEqual(len(data2_items_first), 1)  # opening a list appends a new dictionary

            data2_items_first0 = data2_items_first[0]
            self.assertIsInstance(data2_items_first0, dict)
            self.assertEqual(len(data2_items_first0), 3)

            for i in range(1, 4):
                self.assertIn(f"listkey{i}", data2_items_first0)
                self.assertEqual(data2_items_first0[f"listkey{i}"], f"listval{i}")

            self.assertIn("items-second", data2)
            data2_items_second = data2["items-second"]
            self.assertIsInstance(data2_items_second, list)
            self.assertEqual(len(data2_items_second), 1)  # opening a list appends a new dictionary

            data2_items_second0 = data2_items_second[0]
            self.assertIsInstance(data2_items_second0, dict)
            self.assertEqual(len(data2_items_second0), 3)

            for i in range(4, 7):
                self.assertIn(f"listkey{i}", data2_items_second0)
                self.assertEqual(data2_items_second0[f"listkey{i}"], f"listval{i}")

            # Step 3 - add a list within a list
            data3 = data[3]
            self.assertIsInstance(data3, dict)
            self.assertEqual(len(data3), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data3)
            self.assertEqual(data3["level"], 0)

            for i in range(1, 5):
                self.assertIn(f"key{i}", data3)
                self.assertEqual(data3[f"key{i}"], f"val{i}")

            self.assertIn("items-outer", data3)
            data3_items_outer = data3["items-outer"]
            self.assertIsInstance(data3_items_outer, list)
            self.assertEqual(len(data3_items_outer), 1)  # opening a list appends a new dictionary

            data3_items_outer0 = data3_items_outer[0]
            self.assertIsInstance(data3_items_outer0, dict)
            self.assertEqual(len(data3_items_outer0), 5)

            for i in range(1, 5):
                self.assertIn(f"outerlistkey{i}", data3_items_outer0)
                self.assertEqual(data3_items_outer0[f"outerlistkey{i}"], f"outerlistval{i}")

            self.assertIn("items-inner", data3_items_outer0)
            data3_items_outer0_items_inner = data3_items_outer0["items-inner"]
            self.assertIsInstance(data3_items_outer0_items_inner, list)
            self.assertEqual(
                len(data3_items_outer0_items_inner), 1
            )  # opening a list appends a new dictionary

            data3_items_outer0_items_inner0 = data3_items_outer0_items_inner[0]
            self.assertIsInstance(data3_items_outer0_items_inner0, dict)
            self.assertEqual(len(data3_items_outer0_items_inner0), 3)

            for i in range(1, 4):
                self.assertIn(f"innerlistkey{i}", data3_items_outer0_items_inner0)
                self.assertEqual(
                    data3_items_outer0_items_inner0[f"innerlistkey{i}"], f"innerlistval{i}"
                )

        finally:
            os.unlink(test_log.name)

    def test_mixed_dict_and_list_nesting(self):
        """Test nesting dictionaries within lists and vice versa."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            # Step 0 - list inside dict
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-outer-dict] dictkey1 = dictval1\n"
            "[INFO][rank 0][TestScope][label] dictkey2 = dictval2\n"
            "[INFO][rank 0][TestScope][begin-inner-list-list] listkey1 = listval1\n"
            "[INFO][rank 0][TestScope][label] listkey2 = listval2\n"
            "[INFO][rank 0][TestScope][end-inner-list-list] listkey3 = listval3\n"
            "[INFO][rank 0][TestScope][label] dictkey3 = dictval3\n"
            "[INFO][rank 0][TestScope][begin-inner-list-list] listkey4 = listval4\n"
            "[INFO][rank 0][TestScope][label] listkey5 = listval5\n"
            "[INFO][rank 0][TestScope][end-inner-list-list] listkey6 = listval6\n"
            "[INFO][rank 0][TestScope][label] dictkey4 = dictval4\n"
            "[INFO][rank 0][TestScope][end-outer-dict] dictkey5 = dictval5\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key4 = val4\n"
            # Step 1 - dict inside list
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            "[INFO][rank 0][TestScope][begin-outer-list-list] listkey1 = listval1\n"
            "[INFO][rank 0][TestScope][label] listkey2 = listval2\n"
            "[INFO][rank 0][TestScope][begin-inner-dict] dictkey1 = dictval1\n"
            "[INFO][rank 0][TestScope][label] dictkey2 = dictval2\n"
            "[INFO][rank 0][TestScope][end-inner-dict] dictkey3 = dictval3\n"
            "[INFO][rank 0][TestScope][label] listkey3 = listval3\n"
            "[INFO][rank 0][TestScope][end-outer-list-list] listkey4 = listval4\n"
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            "[INFO][rank 0][TestScope][begin-outer-list-list] listkey5 = listval5\n"
            "[INFO][rank 0][TestScope][label] listkey6 = listval6\n"
            "[INFO][rank 0][TestScope][begin-inner-dict] dictkey4 = dictval4\n"
            "[INFO][rank 0][TestScope][label] dictkey5 = dictval5\n"
            "[INFO][rank 0][TestScope][end-inner-dict] dictkey6 = dictval6\n"
            "[INFO][rank 0][TestScope][label] listkey7 = listval7\n"
            "[INFO][rank 0][TestScope][end-outer-list-list] listkey8 = listval8\n"
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key5 = val5\n"
        )
        test_log.close()

        try:
            data = logs.log_file_to_list(test_log.name)
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 2)

            # Step 0 - list inside dict
            data0 = data[0]
            self.assertIsInstance(data0, dict)
            self.assertEqual(len(data0), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data0)
            self.assertEqual(data0["level"], 0)

            for i in range(1, 5):
                self.assertIn(f"key{i}", data0)
                self.assertEqual(data0[f"key{i}"], f"val{i}")

            self.assertIn("outer-dict", data0)
            data0_outer_dict = data0["outer-dict"]
            self.assertIsInstance(data0_outer_dict, dict)
            self.assertEqual(len(data0_outer_dict), 6)  # opening a list appends a new dictionary

            for i in range(1, 6):
                self.assertIn(f"dictkey{i}", data0_outer_dict)
                self.assertEqual(data0_outer_dict[f"dictkey{i}"], f"dictval{i}")

            self.assertIn("inner-list", data0_outer_dict)
            data0_outer_dict_inner_list = data0_outer_dict["inner-list"]
            self.assertIsInstance(data0_outer_dict_inner_list, list)
            self.assertEqual(
                len(data0_outer_dict_inner_list), 2
            )  # opening a list appends a new dictionary

            data0_outer_dict_inner_list0 = data0_outer_dict_inner_list[0]
            self.assertIsInstance(data0_outer_dict_inner_list0, dict)
            self.assertEqual(len(data0_outer_dict_inner_list0), 3)

            for i in range(1, 4):
                self.assertIn(f"listkey{i}", data0_outer_dict_inner_list0)
                self.assertEqual(data0_outer_dict_inner_list0[f"listkey{i}"], f"listval{i}")

            data0_outer_dict_inner_list1 = data0_outer_dict_inner_list[1]
            self.assertIsInstance(data0_outer_dict_inner_list1, dict)
            self.assertEqual(len(data0_outer_dict_inner_list1), 3)

            for i in range(4, 7):
                self.assertIn(f"listkey{i}", data0_outer_dict_inner_list1)
                self.assertEqual(data0_outer_dict_inner_list1[f"listkey{i}"], f"listval{i}")

            # Step 1 - dict inside list
            data1 = data[1]
            self.assertIsInstance(data1, dict)
            self.assertEqual(len(data1), 7)  # parser adds 1 key (time level)

            self.assertIn("level", data1)
            self.assertEqual(data1["level"], 0)

            for i in range(1, 6):
                self.assertIn(f"key{i}", data1)
                self.assertEqual(data1[f"key{i}"], f"val{i}")

            self.assertIn("outer-list", data1)
            data1_outer_list = data1["outer-list"]
            self.assertIsInstance(data1_outer_list, list)
            self.assertEqual(len(data1_outer_list), 2)  # opening a list appends a new dictionary

            data1_outer_list0 = data1_outer_list[0]
            self.assertIsInstance(data1_outer_list0, dict)
            self.assertEqual(len(data1_outer_list0), 5)

            for i in range(1, 5):
                self.assertIn(f"listkey{i}", data1_outer_list0)
                self.assertEqual(data1_outer_list0[f"listkey{i}"], f"listval{i}")

            self.assertIn("inner-dict", data1_outer_list0)
            data1_outer_list0_dict = data1_outer_list0["inner-dict"]
            self.assertIsInstance(data1_outer_list0_dict, dict)
            self.assertEqual(len(data1_outer_list0_dict), 3)

            for i in range(1, 4):
                self.assertIn(f"dictkey{i}", data1_outer_list0_dict)
                self.assertEqual(data1_outer_list0_dict[f"dictkey{i}"], f"dictval{i}")

            data1_outer_list1 = data1_outer_list[1]
            self.assertIsInstance(data1_outer_list1, dict)
            self.assertEqual(len(data1_outer_list1), 5)

            for i in range(5, 9):
                self.assertIn(f"listkey{i}", data1_outer_list1)
                self.assertEqual(data1_outer_list1[f"listkey{i}"], f"listval{i}")

            self.assertIn("inner-dict", data1_outer_list1)
            data1_outer_list1_dict = data1_outer_list1["inner-dict"]
            self.assertIsInstance(data1_outer_list1_dict, dict)
            self.assertEqual(len(data1_outer_list1_dict), 3)

            for i in range(4, 7):
                self.assertIn(f"dictkey{i}", data1_outer_list1_dict)
                self.assertEqual(data1_outer_list1_dict[f"dictkey{i}"], f"dictval{i}")

        finally:
            os.unlink(test_log.name)

    def test_fast_steps_create_nested_levels(self):
        """Test that begin/end-fast-steps creates nested time levels."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            # Start outer step 1
            "[INFO][rank 0][OuterScope][begin-step-attempt] outerkey1 = outerval1\n"
            "[INFO][rank 0][OuterScope][label] outerkey2 = outerval2\n"
            "[INFO][rank 0][OuterScope][begin-fast-steps]\n"
            # Start middle step 1
            "[INFO][rank 0][MiddleScope][begin-step-attempt] middlekey1 = middleval1\n"
            "[INFO][rank 0][MiddleScope][label] middlekey2 = middleval2\n"
            "[INFO][rank 0][MiddleScope][begin-fast-steps]\n"
            # Inner steps 1
            "[INFO][rank 0][InnerScope][begin-step-attempt] innerkey1 = innerval1\n"
            "[INFO][rank 0][InnerScope][label] innerkey2 = innerval2\n"
            "[INFO][rank 0][InnerScope][end-step-attempt] innerkey3 = innerval3\n"
            "[INFO][rank 0][InnerScope][begin-step-attempt] innerkey4 = innerval4\n"
            "[INFO][rank 0][InnerScope][label] innerkey5 = innerval5\n"
            "[INFO][rank 0][InnerScope][end-step-attempt] innerkey6 = innerval6\n"
            # End middle step 1
            "[INFO][rank 0][MiddleScope][end-fast-steps]\n"
            "[INFO][rank 0][MiddleScope][label] middlekey3 = middleval3\n"
            "[INFO][rank 0][MiddleScope][end-step-attempt] middlekey4 = middleval4\n"
            # Start middle step 2
            "[INFO][rank 0][MiddleScope][begin-step-attempt] middlekey5 = middleval5\n"
            "[INFO][rank 0][MiddleScope][label] middlekey6 = middleval6\n"
            "[INFO][rank 0][MiddleScope][begin-fast-steps]\n"
            # Inner steps 2
            "[INFO][rank 0][InnerScope][begin-step-attempt] innerkey7 = innerval7\n"
            "[INFO][rank 0][InnerScope][label] innerkey8 = innerval8\n"
            "[INFO][rank 0][InnerScope][end-step-attempt] innerkey9 = innerval9\n"
            "[INFO][rank 0][InnerScope][begin-step-attempt] innerkey10 = innerval10\n"
            "[INFO][rank 0][InnerScope][label] innerkey11 = innerval11\n"
            "[INFO][rank 0][InnerScope][end-step-attempt] innerkey12 = innerval12\n"
            # End middle step 2
            "[INFO][rank 0][MiddleScope][end-fast-steps]\n"
            "[INFO][rank 0][MiddleScope][label] middlekey7 = middleval7\n"
            "[INFO][rank 0][MiddleScope][end-step-attempt] middlekey8 = middleval8\n"
            # End outer step 1
            "[INFO][rank 0][OuterScope][end-fast-steps]\n"
            "[INFO][rank 0][OuterScope][label] outerkey3 = outerval3\n"
            "[INFO][rank 0][Outer][end-step-attempt] outerkey4 = outerval4\n"
        )
        test_log.close()

        try:
            data = logs.log_file_to_list(test_log.name)
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 1)

            # Outer step
            data0 = data[0]
            self.assertIsInstance(data0, dict)
            self.assertEqual(len(data0), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data0)
            self.assertEqual(data0["level"], 0)

            for i in range(1, 5):
                self.assertIn(f"outerkey{i}", data0)
                self.assertEqual(data0[f"outerkey{i}"], f"outerval{i}")

            # Middle step 1
            self.assertIn("time-level-1", data0)
            data0_time_level1 = data0["time-level-1"]
            self.assertIsInstance(data0_time_level1, list)
            self.assertEqual(len(data0_time_level1), 2)

            data0_time_level1_step0 = data0_time_level1[0]
            self.assertIsInstance(data0_time_level1_step0, dict)
            self.assertEqual(len(data0_time_level1_step0), 6)  # parser adds 1 key (time level)

            self.assertIn("level", data0_time_level1_step0)
            self.assertEqual(data0_time_level1_step0["level"], 1)

            for i in range(1, 5):
                self.assertIn(f"middlekey{i}", data0_time_level1_step0)
                self.assertEqual(data0_time_level1_step0[f"middlekey{i}"], f"middleval{i}")

            # Inner steps
            self.assertIn("time-level-2", data0_time_level1_step0)
            data0_time_level2 = data0_time_level1_step0["time-level-2"]
            self.assertIsInstance(data0_time_level2, list)
            self.assertEqual(len(data0_time_level2), 2)

            # Inner step 1
            data0_time_level2_step0 = data0_time_level2[0]
            self.assertIsInstance(data0_time_level2_step0, dict)
            self.assertEqual(len(data0_time_level2_step0), 4)  # parser adds 1 key (time level)

            self.assertIn("level", data0_time_level2_step0)
            self.assertEqual(data0_time_level2_step0["level"], 2)

            for i in range(1, 4):
                self.assertIn(f"innerkey{i}", data0_time_level2_step0)
                self.assertEqual(data0_time_level2_step0[f"innerkey{i}"], f"innerval{i}")

            # Inner step 2
            data0_time_level2_step1 = data0_time_level2[1]
            self.assertIsInstance(data0_time_level2_step1, dict)
            self.assertEqual(len(data0_time_level2_step1), 4)  # parser adds 1 key (time level)

            self.assertIn("level", data0_time_level2_step1)
            self.assertEqual(data0_time_level2_step1["level"], 2)

            for i in range(4, 7):
                self.assertIn(f"innerkey{i}", data0_time_level2_step1)
                self.assertEqual(data0_time_level2_step1[f"innerkey{i}"], f"innerval{i}")

            # Middle step 2
            data0_time_level1_step1 = data0_time_level1[1]
            self.assertIsInstance(data0_time_level1_step1, dict)
            self.assertEqual(len(data0_time_level1_step1), 6)

            self.assertIn("level", data0_time_level1_step1)
            self.assertEqual(data0_time_level1_step1["level"], 1)

            for i in range(5, 9):
                self.assertIn(f"middlekey{i}", data0_time_level1_step1)
                self.assertEqual(data0_time_level1_step1[f"middlekey{i}"], f"middleval{i}")

            # Inner steps
            self.assertIn("time-level-2", data0_time_level1_step1)
            data0_time_level2 = data0_time_level1_step1["time-level-2"]
            self.assertIsInstance(data0_time_level2, list)
            self.assertEqual(len(data0_time_level2), 2)

            # Inner step 3
            data0_time_level2_step0 = data0_time_level2[0]
            self.assertIsInstance(data0_time_level2_step0, dict)
            self.assertEqual(len(data0_time_level2_step0), 4)  # parser adds 1 key (time level)

            self.assertIn("level", data0_time_level2_step0)
            self.assertEqual(data0_time_level2_step0["level"], 2)

            for i in range(7, 10):
                self.assertIn(f"innerkey{i}", data0_time_level2_step0)
                self.assertEqual(data0_time_level2_step0[f"innerkey{i}"], f"innerval{i}")

            # Inner step 4
            data0_time_level2_step1 = data0_time_level2[1]
            self.assertIsInstance(data0_time_level2_step1, dict)
            self.assertEqual(len(data0_time_level2_step1), 4)  # parser adds 1 key (time level)

            self.assertIn("level", data0_time_level2_step1)
            self.assertEqual(data0_time_level2_step1["level"], 2)

            for i in range(10, 13):
                self.assertIn(f"innerkey{i}", data0_time_level2_step1)
                self.assertEqual(data0_time_level2_step1[f"innerkey{i}"], f"innerval{i}")

        finally:
            os.unlink(test_log.name)

    def test_partitions_list(self):
        """Test that begin/end-partitions-list creates partition structures."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            # Overall step
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][label] key2 = val2\n"
            # First partition -- creates list
            "[INFO][rank 0][TestScope][begin-partitions-list] parkey1 = parval1\n"
            "[INFO][rank 0][TestScope][label] parkey2 = parval2\n"
            # Steps -- adds evolve entry
            "[INFO][rank 0][TestScope][begin-step-attempt] parstepkey1 = parstepval1\n"
            "[INFO][rank 0][TestScope][label] parstepkey2 = parstepval2\n"
            "[INFO][rank 0][TestScope][end-step-attempt] parstepkey3 = parstepval3\n"
            # First partition -- updates first partition entry
            "[INFO][rank 0][TestScope][label] parkey3 = parval3\n"
            # Steps -- adds new entry to evolve list
            "[INFO][rank 0][TestScope][begin-step-attempt] parstepkey4 = parstepval4\n"
            "[INFO][rank 0][TestScope][label] parstepkey5 = parstepval5\n"
            "[INFO][rank 0][TestScope][end-step-attempt] parstepkey6 = parstepval6\n"
            # First partition -- updates first partition entry
            "[INFO][rank 0][TestScope][label] parkey4 = parval4\n"
            "[INFO][rank 0][TestScope][end-partitions-list] parkey5 = parval5\n"
            # Overall step -- updates overall entry
            "[INFO][rank 0][TestScope][label] key3 = val3\n"
            # Second partition -- adds new entry to partition list
            "[INFO][rank 0][TestScope][begin-partitions-list] parkey6 = parval6\n"
            "[INFO][rank 0][TestScope][label] parkey7 = parval7\n"
            # Steps -- creates evolve entry
            "[INFO][rank 0][TestScope][begin-step-attempt] parstepkey7 = parstepval7\n"
            "[INFO][rank 0][TestScope][label] parstepkey8 = parstepval8\n"
            "[INFO][rank 0][TestScope][end-step-attempt] parstepkey9 = parstepval9\n"
            # Second partition -- updates second partition entry
            "[INFO][rank 0][TestScope][label] parkey8 = parval8\n"
            # Steps -- adds new entry to evolve list
            "[INFO][rank 0][TestScope][begin-step-attempt] parstepkey10 = parstepval10\n"
            "[INFO][rank 0][TestScope][label] parstepkey11 = parstepval11\n"
            "[INFO][rank 0][TestScope][end-step-attempt] parstepkey12 = parstepval12\n"
            # Second partition -- updates second partition entry
            "[INFO][rank 0][TestScope][label] parkey9 = parval9\n"
            "[INFO][rank 0][TestScope][end-partitions-list] parkey10 = parval10\n"
            # Overall step -- updates overall step
            "[INFO][rank 0][TestScope][label] key4 = val4\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key5 = val5\n"
        )
        test_log.close()

        try:
            data = logs.log_file_to_list(test_log.name)
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 1)

            # Overall step
            data0 = data[0]
            self.assertIsInstance(data0, dict)
            self.assertEqual(len(data0), 7)  # parser adds 1 key (time level)

            self.assertIn("level", data0)
            self.assertEqual(data0["level"], 0)

            for i in range(1, 6):
                self.assertIn(f"key{i}", data0)
                self.assertEqual(data0[f"key{i}"], f"val{i}")

            # Partitions
            self.assertIn("partitions", data0)
            data0_partitions = data0["partitions"]
            self.assertIsInstance(data0_partitions, list)
            self.assertEqual(len(data0_partitions), 2)

            # First partition
            data0_partitions0 = data0_partitions[0]
            self.assertIsInstance(data0_partitions0, dict)
            self.assertEqual(len(data0_partitions0), 6)  # time level only added to step attempts

            for i in range(1, 6):
                self.assertIn(f"parkey{i}", data0_partitions0)
                self.assertEqual(data0_partitions0[f"parkey{i}"], f"parval{i}")

            # First partition steps -- begin-step-attempt adds evolve list to dict
            self.assertIn("evolve", data0_partitions0)
            data0_partitions0_evolve = data0_partitions0["evolve"]
            self.assertIsInstance(data0_partitions0_evolve, list)
            self.assertEqual(len(data0_partitions0_evolve), 2)

            # Step 1
            data0_partitions0_evolve0 = data0_partitions0_evolve[0]
            self.assertIsInstance(data0_partitions0_evolve0, dict)
            self.assertEqual(len(data0_partitions0_evolve0), 4)  # parser adds 1 key (time level)

            for i in range(1, 4):
                self.assertIn(f"parstepkey{i}", data0_partitions0_evolve0)
                self.assertEqual(data0_partitions0_evolve0[f"parstepkey{i}"], f"parstepval{i}")

            # Step 2
            data0_partitions0_evolve1 = data0_partitions0_evolve[1]
            self.assertIsInstance(data0_partitions0_evolve1, dict)
            self.assertEqual(len(data0_partitions0_evolve1), 4)  # parser adds 1 key (time level)

            for i in range(4, 7):
                self.assertIn(f"parstepkey{i}", data0_partitions0_evolve1)
                self.assertEqual(data0_partitions0_evolve1[f"parstepkey{i}"], f"parstepval{i}")

            # Second partition
            data0_partitions1 = data0_partitions[1]
            self.assertIsInstance(data0_partitions1, dict)
            self.assertEqual(len(data0_partitions1), 6)  # time level only added to step attempts

            for i in range(6, 11):
                self.assertIn(f"parkey{i}", data0_partitions1)
                self.assertEqual(data0_partitions1[f"parkey{i}"], f"parval{i}")

            # Second partition steps -- begin-step-attempt adds evolve list to dict
            self.assertIn("evolve", data0_partitions1)
            data0_partitions1_evolve = data0_partitions1["evolve"]
            self.assertIsInstance(data0_partitions1_evolve, list)
            self.assertEqual(len(data0_partitions1_evolve), 2)

            # Step 1
            data0_partitions1_evolve0 = data0_partitions1_evolve[0]
            self.assertIsInstance(data0_partitions1_evolve0, dict)
            self.assertEqual(len(data0_partitions1_evolve0), 4)  # parser adds 1 key (time level)

            for i in range(7, 10):
                self.assertIn(f"parstepkey{i}", data0_partitions1_evolve0)
                self.assertEqual(data0_partitions1_evolve0[f"parstepkey{i}"], f"parstepval{i}")

            # Step 2
            data0_partitions1_evolve1 = data0_partitions1_evolve[1]
            self.assertIsInstance(data0_partitions1_evolve1, dict)
            self.assertEqual(len(data0_partitions1_evolve1), 4)  # parser adds 1 key (time level)

            for i in range(10, 13):
                self.assertIn(f"parstepkey{i}", data0_partitions1_evolve1)
                self.assertEqual(data0_partitions1_evolve1[f"parstepkey{i}"], f"parstepval{i}")

        finally:
            os.unlink(test_log.name)

    def test_update_with_duplicate_keys(self):
        """Test that update raises an exception when trying to overwrite an existing key."""
        test_log = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
        test_log.write(
            "[INFO][rank 0][TestScope][begin-step-attempt] key1 = val1\n"
            "[INFO][rank 0][TestScope][end-step-attempt] key1 = val2\n"
        )
        test_log.close()

        try:
            with self.assertRaises(KeyError):
                logs.log_file_to_list(test_log.name)

        finally:
            os.unlink(test_log.name)


def run_tests():
    """Run all tests and print summary."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestLogParsing))

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
