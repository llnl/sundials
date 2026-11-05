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
            if "SUNDIALS_EXAMPLE_VARIANT" in contents or "_variants.cmake" in contents:
                violations.append((cmake_file.parent, "uses shared variant registration"))
            calls = len(re.findall(r"\bsundials_add_example\s*\(", contents, re.I))
            if calls > 1:
                violations.append((cmake_file.parent, f"defines {calls} example variants"))
            if re.search(r"sundials_add_example\s*\(\s*\$", contents):
                violations.append((cmake_file.parent, "uses a variable executable name"))

        for helper in sorted(root.rglob("*_variants.cmake")):
            violations.append((helper, "is a shared variant registration file"))

    for path, message in violations:
        print(f"{path}: {message}")
    return bool(violations)


if __name__ == "__main__":
    raise SystemExit(main())
