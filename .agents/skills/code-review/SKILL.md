---
name: code-review
description: Review SUNDIALS pull requests. Use this during GitHub Copilot code review, when Copilot is requested as a PR reviewer, or when asked to review changed source code or documentation against SUNDIALS developer-guide requirements.
---

# SUNDIALS Code Review

Review only the pull request changes unless the changed code depends on,
extends, or worsens a noncompliant legacy pattern. Prioritize actionable defects
over broad summaries.

Before leaving review comments, read the relevant guide files:

- `doc/superbuild/source/developers/source_code/Rules.rst`
- `doc/superbuild/source/developers/source_code/Style.rst`
- `doc/superbuild/source/developers/source_code/Naming.rst`
- `doc/superbuild/source/developers/documentation/Style.rst`
- `doc/superbuild/source/developers/documentation/Setup.rst` when reviewing documentation layout,
  build instructions, includes, figures, shared documentation, or citations

## Source Review Checklist

Check source changes for SUNDIALS-specific compliance issues:

- Public and private identifiers follow the namespacing, casing, suffix, enum,
  and C++ class/member rules in `doc/developers/source_code/Naming.rst`.
- New C/C++ code is compatible with C99, C++14, and MSVC v1900+ unless a
  third-party API explicitly requires newer features.
- `if` statements and loops use braces, code uses spaces instead of tabs,
  comments use proper spelling and grammar, and TODO comments use the
  `TODO(<identifier>): ...` form.
- SUNDIALS data structures retain a `SUNContext` unless they are an allowed
  exception.
- SUNDIALS functions that can access a context call `SUNFunctionBegin()` as
  early as possible, using the first parameter that owns a `SUNContext`.
- Code outside the `SUNContext` implementation refers to the active context via
  `SUNCTX_`.
- Calls returning `SUNErrCode` are checked with the appropriate `SUNCheckCall`
  family macro.
- Calls to SUNDIALS functions that do not return `SUNErrCode` are followed by
  `SUNCheckLastErr` when required.
- Programmer-error checks use `SUNAssert`.
- New return statements avoid unnecessary parentheses.
- Dimension variables use `sunindextype`, counters use `suncountertype`, memory
  sizes use `size_t`, and other unsigned integer types are avoided, especially
  in public APIs exposed to Fortran.
- Allocations and similar size calculations use `sizeof(variable)` rather than
  `sizeof(type)`.
- Public headers avoid anonymous enum typedefs and use the SWIG-safe enum
  pattern in `doc/developers/source_code/Rules.rst`.
- Output formatting, statistics printing, and logging use the macros and message
  conventions in `doc/developers/source_code/Style.rst`.

## Documentation Review Checklist

Check documentation changes for SUNDIALS-specific compliance issues:

- reStructuredText headings use the levels in
  `doc/developers/documentation/Style.rst`; `#` headings are reserved for the
  documentation superbuild and must not appear in package documentation
  directories.
- Footnotes are not used; use `note` or `warning` directives when appropriate.
- External links usually use anonymous-link syntax with two trailing
  underscores.
- User-callable C and C++ APIs are documented with `c:function` or
  `cpp:function`.
- Function parameters use `:param <name>:` and returns use either
  `:retval <value>:` entries or one `:returns:` entry as appropriate.
- API additions, behavior changes, and deprecations use `versionadded`,
  `versionchanged`, or `deprecated` directives with placeholder version `x.y.z`.
- Shared documentation, figures, package-specific assets, and citations follow
  the organization in `doc/developers/documentation/Setup.rst`.

## Review Comment Format

For each finding, include:

- the affected file and line number
- the violated SUNDIALS developer-guide rule and guide file
- the impact on maintainability, portability, bindings, documentation build, or
  user-visible behavior
- a concise suggested fix

If no SUNDIALS developer-guide compliance issues are found, say so and mention
any review gaps, such as checks that require running formatters, documentation
builds, or tests.
