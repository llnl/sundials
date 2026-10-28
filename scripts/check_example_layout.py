#!/usr/bin/env python3
"""Check that solver examples follow the CVODE problem/variant layout."""

import argparse
import re
from pathlib import Path

PACKAGES = ("cvode", "arkode", "cvodes", "ida", "idas", "kinsol")
SOURCES = {".c", ".cpp", ".cu", ".f90", ".hpp", ".h", ".py"}


def check(root, packages):
    errors = []
    for package in packages:
        directory = root / "examples" / package
        top = directory / "CMakeLists.txt"
        if not top.is_file():
            errors.append(f"{top}: missing package CMakeLists.txt")
            continue
        text = top.read_text()
        for helper in ("sundials_install_example_cmake", "sundials_add_test_install"):
            if helper not in text:
                errors.append(f"{top}: missing {helper}")
        problems = re.findall(r"add_subdirectory\(([^)]+)\)", text)
        for problem in problems:
            problem_dir = directory / problem
            cmake = problem_dir / "CMakeLists.txt"
            if not cmake.is_file():
                errors.append(f"{cmake}: missing problem CMakeLists.txt")
                continue
            variants = re.findall(r"add_subdirectory\(([^)]+)\)", cmake.read_text())
            for variant in variants:
                variant_dir = problem_dir / variant
                variant_cmake = variant_dir / "CMakeLists.txt"
                if not re.fullmatch(
                    r"(?:c|cpp|fortran|cuda|hip|sycl)(?:_[a-z0-9]+)*", variant
                ):
                    errors.append(f"{variant_dir}: invalid language/TPL variant name")
                if not variant_cmake.is_file():
                    errors.append(f"{variant_cmake}: missing variant CMakeLists.txt")
                elif "sundials_add_example" not in variant_cmake.read_text():
                    errors.append(f"{variant_cmake}: missing sundials_add_example")
        for file in directory.rglob("*"):
            if not file.is_file():
                continue
            relative = file.relative_to(directory)
            if file.suffix in SOURCES:
                if len(relative.parts) != 3:
                    errors.append(f"{file}: expected <problem>/<variant>/<source>")
                if not re.fullmatch(r"[a-z][a-z0-9]*(?:_[a-z0-9]+)*", file.stem):
                    errors.append(
                        f"{file}: expected lowercase underscore-separated name"
                    )
                if len(relative.parts) == 3:
                    problem, variant, _ = relative.parts
                    if problem not in problems:
                        errors.append(
                            f"{file}: problem not registered in package CMake"
                        )
                    problem_cmake = directory / problem / "CMakeLists.txt"
                    if (
                        problem_cmake.is_file()
                        and f"add_subdirectory({variant})"
                        not in problem_cmake.read_text()
                    ):
                        errors.append(
                            f"{file}: variant not registered in problem CMake"
                        )
            if file.name.startswith("Makefile"):
                errors.append(f"{file}: obsolete example Makefile")
        for document in (root / "doc" / package).rglob("*.rst"):
            if not document.is_file():
                continue
            for include in re.findall(
                r"literalinclude::\s+(\S+)", document.read_text()
            ):
                if "examples/" in include and not (document.parent / include).is_file():
                    errors.append(f"{document}: missing example include {include}")
    return errors


def check_generated(build):
    errors = []
    for package in PACKAGES:
        for cmake in (build / "examples" / package).glob("*/*/CMakeLists.txt"):
            if "add_executable(" not in cmake.read_text():
                errors.append(f"{cmake}: empty installed example project")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--packages", nargs="+", choices=PACKAGES, default=PACKAGES)
    parser.add_argument(
        "--build", type=Path, help="Also check generated standalone projects"
    )
    args = parser.parse_args()
    errors = check(args.root, args.packages)
    if args.build:
        if not (args.build / "CMakeCache.txt").is_file():
            parser.error("--build must name a configured CMake build")
        errors.extend(check_generated(args.build))
    for error in errors:
        print(error)
    if errors:
        raise SystemExit(1)
    print("Example layout checks passed: " + ", ".join(args.packages))


if __name__ == "__main__":
    main()
