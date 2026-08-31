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

Check source changes for SUNDIALS-specific compliance issues.
Check documentation changes for SUNDIALS-specific compliance issues.

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
