#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody Balos and David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Module of Python functions that may be useful to SUNDIALS developers writing
# scripts to parse logs produced by SUNLogger.
# -----------------------------------------------------------------------------

import re
import json
from typing import Iterable


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
        kvp = kvpstr.split("=", 1)  # only split on first =
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
                    values.append(_convert_to_num(line.strip()))
                kvp_dict[key.strip()] = values
            else:
                kvp_dict[key.strip()] = _convert_to_num(value.strip())
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
    pattern = re.compile(r"\[(\w+)\]\[(rank \d+)\]\[(.*?)\]\[(.*?)\](.*)")
    matches = pattern.findall(line)
    line_dict = {}
    if matches:
        line_dict["loglvl"] = matches[0][0]
        line_dict["rank"] = _convert_to_num(matches[0][1].split()[1])
        line_dict["scope"] = matches[0][2]
        line_dict["label"] = matches[0][3]
        line_dict["payload"] = _parse_logfile_payload(matches[0][4], line_number, all_lines)
    return line_dict


def _label_region_info(label):
    """Return hierarchical region info derived from a SUNLogger label.

    The `log_file_to_list` parser treats labels of the form:

    - `begin-<region>` / `end-<region>` as dictionary regions
    - `begin-<region>-list` / `end-<region>-list` as list-of-dicts regions

    This helper centralizes the label parsing logic so other tools (e.g., the
    `suntools` CLI) can stay consistent with :py:func:`log_file_to_list`.

    :param str label: The SUNLogger label.
    :returns: (event, region_name, is_list) where event is "begin", "end", or None.
    :rtype: tuple[str|None, str|None, bool]
    """
    if not label:
        return None, None, False

    label_split = label.split("-")
    if not label_split:
        return None, None, False

    event = label_split[0]
    if event not in {"begin", "end"}:
        return None, None, False

    is_list = label_split[-1] == "list"
    if is_list:
        region_name = "-".join(label_split[1:-1])
    else:
        region_name = "-".join(label_split[1:])

    return event, region_name, is_list


def _split_filter_list(value: str) -> set[str]:
    parts = [p.strip().lower() for p in value.split(",")]
    return {p for p in parts if p}


def _is_log_line(line: str) -> bool:
    return line.startswith("[")


def _find_continuation_end(lines: list[str], start_index: int) -> int:
    """Return the index of the next log line or blank line (exclusive)."""
    i = start_index + 1
    while i < len(lines):
        if _is_log_line(lines[i]) or not lines[i].strip():
            break
        i += 1
    return i


def _iter_filtered_lines(
    lines: list[str], selected: set[str], invert: bool = False
) -> Iterable[str]:
    step_depth = 0
    nonlinear_depth = 0
    linear_depth = 0

    i = 0
    while i < len(lines):
        line = lines[i]
        if not _is_log_line(line):
            # Attribute non-log lines to the current region context.
            cats: set[str] = set()
            if linear_depth > 0:
                cats.add("linear")
            elif nonlinear_depth > 0:
                cats.add("nonlinear")
            elif step_depth > 0:
                cats.add("integrator")

            keep = bool(cats & selected)
            if invert:
                keep = not keep
            if keep:
                yield line

            i += 1
            continue

        entry = _parse_logfile_line(line.rstrip("\n"), i, lines)
        if not entry:
            i += 1
            continue

        label = entry.get("label", "")
        event, region_name, _ = _label_region_info(label)

        # Hierarchical region tracking:
        # - everything between begin/end-nonlinear-solve is "nonlinear"
        # - everything between begin/end-linear-solve is "linear"
        # - when filtering "nonlinear", exclude any nested "linear" region
        # - "integrator" is everything within begin/end-step-attempt excluding
        #   nonlinear/linear regions.
        #
        # Treat begin/end markers as part of their own region.
        step_active = step_depth > 0 or (region_name == "step-attempt" and event is not None)
        nonlinear_active = nonlinear_depth > 0 or (
            region_name == "nonlinear-solve" and event == "begin"
        )
        linear_active = linear_depth > 0 or (region_name == "linear-solve" and event == "begin")

        cats: set[str] = set()

        if linear_active:
            cats.add("linear")
        elif nonlinear_active:
            cats.add("nonlinear")

        if step_active and ("linear" not in cats) and ("nonlinear" not in cats):
            cats.add("integrator")

        keep = bool(cats & selected)
        if invert:
            keep = not keep

        # Only treat following non-log lines as a continuation block when the
        # log line is an array/vector dump (payload values parsed as lists).
        is_array_dump = any(isinstance(v, list) for v in entry.get("payload", {}).values())
        continuation_end = _find_continuation_end(lines, i) if is_array_dump else i + 1
        if keep:
            for out_line in lines[i:continuation_end]:
                yield out_line

        # Update region nesting after processing/printing so the begin line
        # itself is considered part of the region.
        if region_name == "step-attempt":
            if event == "begin":
                step_depth += 1
            elif event == "end":
                step_depth = max(0, step_depth - 1)

        if region_name == "nonlinear-solve":
            if event == "begin":
                nonlinear_depth += 1
            elif event == "end":
                nonlinear_depth = max(0, nonlinear_depth - 1)

        if region_name == "linear-solve":
            if event == "begin":
                linear_depth += 1
            elif event == "end":
                linear_depth = max(0, linear_depth - 1)

        i = continuation_end


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
            values.append(entry[key])
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


class StepData:
    """Helper class for parsing a step attempt from a log file into a hierarchical
    dictionary where entries may be lists of dictionaries.
    """

    def __init__(self):
        self.stack = [{}]

    def __repr__(self):
        return json.dumps(self.stack[0], indent=2)

    def update(self, data):
        """Update the active dictionary"""
        conflicts = self.stack[-1].keys() & data.keys()
        if conflicts:
            raise KeyError(f"Cannot update: keys already exist: {conflicts}")
        self.stack[-1].update(data)

    def open_dict(self, key):
        """Activate a nested dictionary"""
        if key not in self.stack[-1]:
            self.stack[-1][key] = {}
        self.stack.append(self.stack[-1][key])

    def close_dict(self):
        """Deactivate the active dictionary"""
        self.stack.pop()

    def open_list(self, key):
        """Activate a list of dictionaries"""
        if key not in self.stack[-1]:
            self.stack[-1][key] = []
        new_dict = {}
        self.stack[-1][key].append(new_dict)
        self.stack.append(new_dict)

    def close_list(self):
        """Deactivate the active list"""
        self.stack.pop()

    def get_step(self):
        """Get the step dictionary and reset the container"""
        result = self.stack[0]
        self.stack = [{}]
        return result


def log_file_to_list(filename):
    """Parse a log file and return as a Python list of dictionaries.

    This is the core parsing function that converts SUNDIALS log files into Python
    data structures. Use this function when you need to work with the data directly
    in Python (analysis, filtering, custom processing).

    :param str filename: The name of the log file to parse.
    :returns: A list of dictionaries, one per step attempt.
    :rtype: list[dict]

    The list returned for a time integrator log file will contain a dictionary for each
    step attempt:

    .. code-block:: python

       [
         {
           "step": 1,
           "tn": 0.0,
           "h": 0.01,
           "status": "success",
           "dsm": 2.6e-13,
           "level": 0,
           "stages": [
             {"stage": 0, "implicit": 0, "tcur": 0.0, "status": "success"},
             {"stage": 1, "implicit": 0, "tcur": 0.001, "status": "success"},
             ...
           ],
           "compute-solution": {
             "mass type": 0,
             "status": "success"
           }
         },
         {
           "step": 2,
           "tn": 0.01,
           "h": 0.02,
           ...
         },
         ...
       ]

    **Example usage:**

    .. code-block:: python

       from suntools import logs
       import matplotlib.pyplot as plt

       # Parse the log file
       data = logs.log_file_to_list("sun.log")

       # Extract lists of passed and failed step data
       passed = [s for s in data if "success" in s["status"]]
       failed = [s for s in data if "failed" in s["status"]]

       print("Steps stats: ")
       print(f"  Attempted:  {len(data)}")
       print(f"  Successful: {len(passed)}")
       print(f"  Failed:     {len(failed)}")
       print(f"  Min Step:   {min(s["h"] for s in passed)}")
       print(f"  Max Step:   {max(s["h"] for s in passed)}")
       print(f"  Avg Step:   {sum(s["h"] for s in passed) / len(passed)}")

       # Plot step size history
       s_steps, s_times, s_steps = logs.get_history(passed, "h")
       f_steps, f_times, f_steps = logs.get_history(failed, "h")

       fig, ax = plt.subplots()
       ax.plot(s_times, s_steps, color="b", marker=".", label="passed")
       ax.scatter(f_times, f_steps, color="r", marker="x", label="failed")
       ax.legend(loc="best")
       ax.set(xlabel="time", ylabel="step size", title="Step Size History")
       plt.show()
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
            event, region_name, is_list = _label_region_info(label)

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

            if event == "begin":
                if is_list:
                    s.open_list(region_name)
                else:
                    s.open_dict(region_name)
                s.update(line_dict["payload"])
                continue
            elif event == "end":
                s.update(line_dict["payload"])
                if is_list:
                    s.close_list()
                else:
                    s.close_dict()
                continue

            s.update(line_dict["payload"])

    return step_attempts


def print_log(log, indent=2):
    """Print a log file list as formatted JSON.

    This function takes parsed log data and prints it in a human-readable JSON
    format. Useful for debugging and quick inspection.

    :param list log: The log file list from :py:func:`log_file_to_list()`.
    :param int indent: The number of spaces to indent the JSON output (default: 2).

    **Example usage:**

    .. code-block:: python

       from suntools import logs

       # Parse the log file
       data = logs.log_file_to_list("sun.log")

       # Print just the first few steps
       logs.print_log(data[:5])
    """
    print(json.dumps(log, indent=indent))


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
