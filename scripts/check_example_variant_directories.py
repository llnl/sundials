#!/usr/bin/env python3

"""Check example leaf directories and their documentation references."""

import argparse
import re
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("packages", nargs="+")
    parser.add_argument("--check-docs", action="store_true")
    args = parser.parse_args()

    violations = []
    for package in args.packages:
        root = Path("examples") / package
        for problem in sorted(root.iterdir()):
            if not problem.is_dir():
                continue
            problem_name = problem.name.split("_", 1)[-1]
            for variant in sorted(problem.iterdir()):
                if variant.is_dir() and re.search(
                    rf"(?:^|_){re.escape(problem_name)}(?:_|$)", variant.name
                ):
                    violations.append(
                        (variant, "repeats the problem name in the variant directory")
                    )
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

        documents = (Path("doc") / package).rglob("*.tex") if args.check_docs else []
        for document in sorted(documents):
            for reference in re.findall(
                rf"(?:\.\./)+examples/{re.escape(package)}/[^\s{{}}]+", document.read_text()
            ):
                if not (document.parent / reference).is_file():
                    violations.append((document, f"references missing file {reference}"))

    for path, message in violations:
        print(f"{path}: {message}")
    return bool(violations)


if __name__ == "__main__":
    raise SystemExit(main())
