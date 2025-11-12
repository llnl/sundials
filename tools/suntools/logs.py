#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody Balos and David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# -----------------------------------------------------------------------------
# Module of Python functions that may be useful to SUNDIALS developers writing
# scripts to parse logs produced by SUNLogger.
# -----------------------------------------------------------------------------

import re
from collections import ChainMap


def _convert_to_num(s):
    """Try to convert a string to an int or float

    :param str s: The string to convert.
    :returns: If the string is a numerical value, an integer (long long) or floating
              point (double) value, otherwise the input string.
    """
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def _parse_logfile_payload(payload, line_number, all_lines, array_indicator="(:)"):
    """Parse the payload of a log file line into a dictionary.

    :param str payload: The payload of a log file line.
    :param int line_number: The line number of payload in the log file.
    :param str all_lines: All the lines in the log file.
    :param str array_indicator: The string that denotes an array output in the log file.
    :returns: A dictionary of key-value pairs from the payload.
    """
    kvpstrs = payload.split(",")
    kvp_dict = {}
    for kvpstr in kvpstrs:
        kvp = kvpstr.split("=")
        if len(kvp) == 1:
            # Check for empty payload
            if not kvp[0].strip():
                continue
            kvp_dict[kvp[0].strip()] = ""
        else:
            key, value = kvp
            values = []
            if array_indicator in key:
                for line in all_lines[line_number + 1 :]:
                    if line.startswith("[") or not line.strip():
                        break
                    values.append(float(line))
                kvp_dict[key.strip()] = values
            else:
                kvp_dict[key.strip()] = value.strip()
    return kvp_dict


def _parse_logfile_line(line, line_number, all_lines):
    """Parse a line from a log file into a dictionary.

    :param str line: The log file line to parse.
    :param int line_number: The line number of the line in the log file.
    :param str all_lines: All the lines in the log file.
    :returns: A dictionary of key-value pairs from the line payload.

    A log file line begins a preamble containing the logging level (ERROR, WARNING,
    INFO, or DEBUG), the MPI rank that issued the message, the function that issued the
    message (scope), and a label with additional context for the message. For
    informational or debugging logs the preamble is followed by the payload which is
    either a comma-separated list of key-value pairs

    .. code-block:: none

       [loglvl][rank][scope][label] key1 = value1, key2 = value2

    or multiline output with one value per line for keys corresponding to a vector or
    array

    .. code-block:: none

       [loglvl][rank][scope][label] y(:) =
       y_1
       y_2
       ...
    """
    pattern = re.compile(r"\[(\w+)\]\[(rank \d+)\]\[(.*)\]\[(.*)\](.*)")
    matches = pattern.findall(line)
    line_dict = {}
    if matches:
        line_dict["loglvl"] = matches[0][0]
        line_dict["rank"] = _convert_to_num(matches[0][1].split()[1])
        line_dict["scope"] = matches[0][2]
        line_dict["label"] = matches[0][3]
        line_dict["payload"] = _parse_logfile_payload(matches[0][4], line_number, all_lines)
    return line_dict


class StepData:
    """Helper class for parsing a step attempt from a log file into a hierarchical
    dictionary where entries may be lists of dictionaries.
    """

    def __init__(self):
        self.container = [ChainMap()]
        self.parent_keys = ["main"]

    def __repr__(self):
        tmp = "Container:"
        for entry in self.container:
            tmp += f"\n  {entry}"
        tmp += "\nParent Keys:"
        for entry in self.parent_keys:
            tmp += f"\n  {entry}"
        return tmp

    def update(self, data):
        """Update the active dictionary"""
        self.container[-1].update(data)

    def open_dict(self, key):
        """Activate a dictionary"""
        self.container[-1][key] = {}
        self.container[-1] = self.container[-1].new_child(self.container[-1][key])

    def close_dict(self):
        """Deactivate the active dictionary"""
        self.container[-1] = self.container[-1].parents

    def open_list(self, key):
        """Activate a list of dictionaries"""
        if key in self.container[-1]:
            self.container.append(ChainMap())
        else:
            self.container[-1][key] = []
            self.container.append(ChainMap())
        self.parent_keys.append(key)

    def close_list(self):
        """Deactivate the active list"""
        tmp = self.container[-1].maps[0]
        self.container[-2][self.parent_keys[-1]].append(tmp)
        self.parent_keys.pop()
        self.container.pop()

    def get_step(self):
        """Get the step dictionary and reset the container"""
        tmp = self.container.pop().maps[0]
        self.container = [ChainMap()]
        self.parent_keys = ["main"]
        return tmp


def log_file_to_list(filename):
    """Parses a log file and returns a list of dictionaries.

    :param str filename: The name of the log file to parse.
    :returns: A list of dictionaries.

    The list returned for a time integrator log file will contain a dictionary for each
    step attempt e.g.,

    .. code-block:: none

       [ {step : 1, tn : 0.0, h : 0.01, ...}, {step : 2, tn : 0.01, h : 0.10, ...}, ...]
    """
    with open(filename, "r") as logfile:

        # List of step attempts, each entry is a dictionary for one attempt
        step_attempts = []

        # Time level for nested integrators e.g., MRI methods
        level = 0

        # Partition for split integrators e.g., operator splitting methods
        partition = 0

        # Read the log file
        all_lines = logfile.readlines()

        # Create instance of helper class for building attempt dictionary
        s = StepData()

        for line_number, line in enumerate(all_lines):

            line_dict = _parse_logfile_line(line.rstrip(), line_number, all_lines)

            if not line_dict:
                continue

            label = line_dict["label"]
            label_split = label.split("-")
            list_or_dict = label_split[-1]
            if list_or_dict == "list":
                region_name = "-".join(label_split[1:-1])
            else:
                region_name = "-".join(label_split[1:])

            if label == "begin-step-attempt":
                line_dict["payload"]["level"] = level
                if level > 0:
                    s.open_list(f"time-level-{level}")
                if partition > 0:
                    s.open_list(f"evolve")
                s.update(line_dict["payload"])
                continue
            elif label == "end-step-attempt":
                s.update(line_dict["payload"])
                if level > 0 or partition > 0:
                    s.close_list()
                else:
                    step_attempts.append(s.get_step())
                continue

            if label == "begin-fast-steps":
                level += 1
                continue
            elif label == "end-fast-steps":
                level -= 1
                continue

            if label == "begin-partitions-list":
                s.open_list("partitions")
                s.update(line_dict["payload"])
                partition += 1
                continue
            elif label == "end-partitions-list":
                s.update(line_dict["payload"])
                s.close_list()
                partition -= 1
                continue

            if label_split[0] == "begin":
                if list_or_dict == "list":
                    s.open_list(region_name)
                else:
                    s.open_dict(region_name)
                s.update(line_dict["payload"])
                continue
            elif label_split[0] == "end":
                s.update(line_dict["payload"])
                if list_or_dict == "list":
                    s.close_list()
                else:
                    s.close_dict()
                continue

            s.update(line_dict["payload"])

    return step_attempts


def _print_log_list(a_list, indent=0):
    """Print list value from a log entry dictionary"""
    spaces = (indent + 2) * " "
    for entry in a_list:
        if type(entry) is list:
            print(f"{spaces}[")
            _print_log_list(entry, indent + 2)
            print(f"{spaces}]")
        elif type(entry) is dict:
            print(f"{spaces}{{")
            _print_log_dict(entry, indent + 2)
            print(f"{spaces}}}")
        else:
            print(f"{spaces}{entry}")


def _print_log_dict(a_dict, indent=0):
    """Print dictionary value from a log entry dictionary"""
    spaces = (indent + 2) * " "
    for key in a_dict:
        if type(a_dict[key]) is list:
            print(f"{spaces}{key} :")
            print(f"{spaces}[")
            _print_log_list(a_dict[key], indent=indent + 2)
            print(f"{spaces}]")
        elif type(a_dict[key]) is dict:
            print(f"{spaces}{key} :")
            print(f"{spaces}{{")
            _print_log_dict(a_dict[key], indent=indent + 2)
            print(f"{spaces}}}")
        else:
            print(f"{spaces}{key} : {a_dict[key]}")


def print_log(log, indent=0):
    """Print a log file list created by :py:func:`log_file_to_list`.

    :param list log: The log file list to print.
    :param int indent: The number of spaces to indent the output.
    """
    spaces = indent * " "
    subspaces = (indent + 2) * " "
    for entry in log:
        print(f"{spaces}{{")
        for key in entry:
            if type(entry[key]) is list:
                print(f"{subspaces}{key} :")
                print(f"{subspaces}[")
                _print_log_list(entry[key], indent=indent + 2)
                print(f"{subspaces}]")
            elif type(entry[key]) is dict:
                print(f"{subspaces}{key} :")
                print(f"{subspaces}{{")
                _print_log_dict(entry[key], indent=indent + 2)
                print(f"{subspaces}}}")
            else:
                print(f"{subspaces}{key} : {entry[key]}")
        print(f"{spaces}}}")


def get_history(
    log, key, step_status=None, time_range=None, step_range=None, group_by_level=False
):
    """Extract the history of a key from a log file list created by
    :py:func:`log_file_to_list`.

    :param list log: The log file list to extract values from.
    :param str key: The key to extract.
    :param str step_status: Only extract values for steps which match the given status
                            e.g., "success" or "failed".
    :param time_range: Only extract values in the time interval, [low, high].
    :type time_range: [float, float]
    :param step_range: Only extract values in the step number interval, [low, high].
    :type step_range: [int, int]
    :param bool group_by_level: Group outputs by time level.
    :returns: A list of steps, times, and values
    """

    steps, times, values, levels = _get_history(log, key, step_status, time_range, step_range)

    if group_by_level:
        from collections import defaultdict

        steps_by_level = defaultdict(list)
        times_by_level = defaultdict(list)
        values_by_level = defaultdict(list)
        for s, t, v, l in zip(steps, times, values, levels):
            steps_by_level[l].append(s)
            times_by_level[l].append(t)
            values_by_level[l].append(v)
        return steps_by_level, times_by_level, values_by_level
    else:
        return steps, times, values


def _get_history(log, key, step_status, time_range, step_range):
    """Extract the step/time series of the requested value."""

    steps = []
    times = []
    values = []
    levels = []

    for entry in log:

        step = int(entry["step"])
        time = float(entry["tn"])
        level = entry["level"]

        if time_range is not None:
            if time < time_range[0] or time > time_range[1]:
                continue

        if step_range is not None:
            if step < step_range[0] or step > step_range[1]:
                continue

        save_data = True
        if step_status is not None:
            if step_status not in entry["status"]:
                save_data = False

        if key in entry and save_data:
            steps.append(step)
            times.append(time)
            values.append(_convert_to_num(entry[key]))
            levels.append(level)

        if "stages" in entry:
            for s in entry["stages"]:
                next_level_key = f"time-level-{level + 1}"
                if next_level_key in s:
                    sub_steps, sub_times, sub_values, sub_levels = _get_history(
                        s[next_level_key], key, step_status, time_range, None
                    )
                    steps.extend(sub_steps)
                    times.extend(sub_times)
                    values.extend(sub_values)
                    levels.extend(sub_levels)

        if "compute-embedding" in entry:
            next_level_key = f"time-level-{level + 1}"
            if next_level_key in entry["compute-embedding"]:
                sub_steps, sub_times, sub_values, sub_levels = _get_history(
                    entry["compute-embedding"][next_level_key], key, step_status, time_range, None
                )
                steps.extend(sub_steps)
                times.extend(sub_times)
                values.extend(sub_values)
                levels.extend(sub_levels)

    return steps, times, values, levels
