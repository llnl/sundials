# Historical SUNDIALS reviewer patterns

This reference calibrates the `code-review` skill. It is evidence about the
repository’s historical review behavior, not a set of mandatory labels or a
substitute for reading the current developer guides.

## Audit scope

On 2026-09-01, `gh api graphql` was used against `LLNL/sundials` with all
`OPEN`, `CLOSED`, and `MERGED` pull requests, including pagination for review
and review-thread connections:

- 736 PRs: 601 merged, 118 closed, 17 open
- 707 general conversation comments
- 5,228 submitted reviews
- 12,236 inline review comments, including thread replies
- 13 bot comments were excluded from the human-reviewer behavior analysis
- no nested review-thread connection was truncated

Topic counts below are overlapping lexical indicators, not hand-labeled
categories. They show where to look first, not how to decide whether a finding
is valid.

| Reviewer lens | Comments | PRs |
| --- | ---: | ---: |
| Error handling and edge cases | 1,970 | 157 |
| Style and maintainability | 1,254 | 172 |
| Documentation and examples | 1,149 | 188 |
| Clarification and requested changes | 803 | 179 |
| Testing and reproducibility | 714 | 188 |
| Memory and resource ownership | 642 | 132 |
| Portability and build integration | 610 | 163 |
| API and compatibility | 570 | 104 |
| Performance and parallelism | 361 | 86 |
| Design and scope | 203 | 69 |

## What reviewers repeatedly do

### Correctness and contracts

Reviewers trace return values and ask whether a setter’s success or failure is
being ignored, whether an error is translated consistently, and whether a
cached value becomes stale after its inputs change. They pay close attention to
null-safe cleanup, partially initialized objects, allocation ownership, and
behavior for insufficient or invalid user input.

### API and maintainability

They question whether a symbol belongs in a public header, whether a new
routine follows the existing prefix/default/`Set` convention, and whether a
shortcut creates a surprising `NULL` or default object. They often ask for a
design explanation when a change introduces a new abstraction, duplicates
module-specific logic, or hides a context/API dependency.

### Tests, numerical behavior, and integration

They request a focused regression test, reproduction details for platform-only
failures, and an explanation of answer-file or roundoff changes. They inspect
precision, MPI/GPU, compiler, Fortran/SWIG, CMake, installation, and package
integration when a change can affect those variants.

### Documentation and examples

They expect documentation to state defaults, parameter meaning, ownership,
errors, constraints, backend limitations, and notation consistently with the
code. User-visible changes commonly trigger requests for both changelog and
`RecentChanges.rst` updates. Reviewers distinguish examples intended to teach
users from examples that primarily provide regression coverage.

## High-signal detailed review patterns

A focused subanalysis of recurring maintainer comments found these useful
checks, which are folded into the main skill rather than exposed as a separate
reviewer profile:

- Follow return flags through wrappers and setters. Check null and invalid
  input, partial cleanup, and whether cached RHS, Jacobian, or linearization
  data is invalidated after its inputs change.
- For integrators, inspect stage loops, forcing calls, post-processing,
  interpolation points, final-state handling, step/order transitions, nested
  steppers, and the exact values passed to user callbacks.
- Compare names, defaults, `Set` routines, registration tables, accessors,
  public/private exposure, and user-data semantics across sibling modules.
  Avoid copy-paste mistakes, unnecessary APIs, and changes that reduce future
  implementation flexibility.
- Keep source, headers, and parallel documentation pages synchronized. When a
  default changes, update all corresponding tables and older/newer module
  variants with `versionchanged` notes.
- Distinguish instructional examples from regression tests. Keep local example
  output separate from CI-generated answer files and check deprecated accessors
  in new examples.
- Consider TPLs off, C90/header checks, CMake installation, compilers, MPI,
  GPU/device pointers, UVM, precision modes, Fortran, and Python/SWIG when
  relevant; a local default build is not enough evidence.

In the focused sample, 417 PRs contained the relevant reviewer activity:
188 conversation comments, 1,512 submitted reviews, and 6,656 inline comments.
Among 7,009 substantive comments, question indicators appeared in 1,311,
request language in 1,312, softeners in 1,249, concern language in 369,
positive language in 257, and code-suggestion markers in 3,316. These
indicators overlap heavily; the important signal is the combination of a
precise observation, a reason, and a proposed correction.

Characteristic concrete checks included:

- reject a null RHS and check the initialization return flag before continuing
  in a binding wrapper;
- make destruction safe for a null coefficient object and partially allocated
  nested arrays;
- replace deprecated `NV_Ith_S` access in examples with
  `N_VGetArrayPointer`;
- add `versionchanged` notes and update every corresponding documentation
  table when a default changes;
- warn users when a permitted UVM or device-pointer configuration may not be
  supported by their system;
- run a TPL-off/C90 header compatibility test and keep CI-generated outputs
  separate from local example output;
- avoid resolving review threads before the replacement has been pushed.

These are examples of reasoning style, not instructions to flag every
occurrence of a matching keyword.

## Tone and comment shape

The historical review style is collaborative and question-led. Across human
substantive comments, questions and softeners are common, but concrete concern
and blocker language is also frequent. A useful pattern is:

1. point to the exact behavior or line;
2. ask whether the intended contract is clear;
3. explain the failure mode or future maintenance cost;
4. propose a small fix, code suggestion, test, or documentation update.

Keep uncertainty explicit, but classify a demonstrated correctness or safety
problem as blocking. Treat the topic counts as prioritization evidence only;
do not flag a line merely because it contains words such as `error`, `NULL`, or
`test`.
