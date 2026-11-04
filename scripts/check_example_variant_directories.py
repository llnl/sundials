#!/usr/bin/env python3

"""Require one sundials_add_example executable per example leaf directory."""

import argparse
import re
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("packages", nargs="+")
    args = parser.parse_args()

    violations = []
    for package in args.packages:
        root = Path("examples") / package
        for cmake_file in sorted(root.rglob("CMakeLists.txt")):
            contents = cmake_file.read_text()
            calls = len(re.findall(r"\bsundials_add_example\s*\(", contents, re.I))
            if calls > 1:
                violations.append((cmake_file.parent, calls))

    for directory, calls in violations:
        print(f"{directory}: defines {calls} example variants")
    return bool(violations)


if __name__ == "__main__":
    raise SystemExit(main())
