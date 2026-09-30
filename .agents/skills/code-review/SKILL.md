---
name: code-review
description: Review SUNDIALS pull requests using the repository’s correctness, numerical, API, documentation, testing, portability, and maintainability review patterns. Use for GitHub PR reviews or when asked to review changed SUNDIALS source, bindings, build files, examples, or documentation.
---

# SUNDIALS PR Review

Review the pull request as a SUNDIALS maintainer. The historical reviewer
patterns are summarized in [references/reviewer-patterns.md](references/reviewer-patterns.md);
read it when choosing review lenses and calibrating tone.

Review the diff first. Inspect surrounding code, public headers, callers,
tests, documentation, and build metadata only when they establish the impact
of a changed line. Do not report pre-existing problems unless the change
introduces them, relies on them, or makes them worse.

Before commenting, read the relevant developer guidance:

- `doc/developers/source_code/Rules.rst`
- `doc/developers/source_code/Style.rst`
- `doc/developers/source_code/Naming.rst`
- `doc/developers/documentation/Style.rst`
- `doc/developers/documentation/Setup.rst` for documentation layout, builds,
  includes, figures, shared documentation, or citations

## Review order

Use the lenses that match the changed files. Prioritize correctness and
user-visible behavior, then numerical and lifecycle contracts, API and
maintainability, integration and portability, tests and reproducibility,
documentation, and finally style. Do not spend review attention on cosmetic
nits while a functional concern is unresolved.

### Correctness, errors, and edge cases

- Trace success and failure paths, especially every `SUNErrCode` return. Check
  calls with the appropriate `SUNCheckCall` family macro and propagate or
  translate failures consistently. Check `SUNCheckLastErr` where required for
  SUNDIALS calls without a `SUNErrCode` return.
- Check null, zero, negative, empty, overflow, insufficient-history, and
  allocation-failure cases when the new path can receive them. Cleanup must be
  safe for partially initialized and null inputs.
- Follow return flags from wrappers, setters, and helper functions. Confirm
  that errors are checked before continuing and that new warnings or
  debug/release behavior are documented for users.
- Look for stale cached state after setters or changes to `(t, y)`, vectors,
  tolerances, linearization points, or user data. Confirm that defaults and
  reset/reuse behavior remain coherent.
- Use `SUNAssert` for programmer-error checks and preserve the documented
  release/debug semantics of assertions and error checking.
- Ensure files are placed in the paths consistent with existing practices.

### Numerical and integrator behavior

- For integrators, trace stage loops, forcing calls, post-processing,
  interpolation points, final-state handling, step/order transitions, and
  nested-stepper interactions. Confirm that final values are neither
  accidentally processed nor skipped before a later RHS call.
- Check norms, tolerances, precision conversions, integer/real types,
  roundoff sensitivity, answer-file changes, and platform assumptions. A
  passing result on one machine is not enough for a numerically sensitive path.
- Check that RHS, Jacobian, or linearization caches are invalidated whenever
  their inputs change, and verify the values passed to user callbacks.

### API, context, compatibility, and maintainability

- Check public versus private exposure, installed headers, bindings, and ABI
  consequences. New public behavior needs a clear default, ownership contract,
  error contract, and compatibility story.
- Check names, prefixes, defaults, `Set` routine conventions, enum patterns,
  accessor usage, deprecations, and registration tables against the existing
  module and `doc/developers/source_code/Naming.rst`. Avoid introducing a
  second convention for the same abstraction.
- Question unnecessary APIs, copy-paste mistakes, duplicated module logic,
  hidden context dependencies, and changes that unnecessarily limit future
  implementation flexibility.
- Preserve `SUNContext` in SUNDIALS data structures unless the documented
  exception applies. Functions that can access a context should call
  `SUNFunctionBegin()` as early as possible using the first owning parameter;
  code outside the context implementation should use `SUNCTX_` as specified.
- Check that public headers remain C99, C++14, MSVC v1900+, SWIG, and Fortran
  compatible where applicable. Dimension values use `sunindextype`, counters
  use `suncountertype`, sizes use `size_t`, and allocations prefer
  `sizeof(variable)`.

### Ownership, build, and portability integration

- Follow allocations through all exits. Check ownership, lifetime, user-data
  pointers, nested objects, null-safe destruction, and cleanup after partial
  construction.
- Check CMake targets, installation rules, generated interfaces, examples,
  MPI/GPU/precision variants, compiler portability, and package/Spack layout
  when those paths are touched. A feature is incomplete if it works only in
  the default build.
- When relevant, consider TPLs enabled and disabled, C90/header compatibility
  checks, device pointers, UVM, Fortran, Python/SWIG bindings, and precision
  modes. A local default build is not sufficient evidence for these paths.
- For examples, distinguish instructional examples from regression tests.
  Keep instructional examples understandable and keep regression coverage
  reproducible across supported precisions and platforms.

### Tests and documentation

- Expect a focused test or a clear reason one is not practical for changed
  behavior. Check answer files and output tolerances when numerical output can
  change. When a failure is platform-specific, request a reproducible case or
  evidence from the affected CI configuration.
- Keep local example output separate from CI-generated answer files, and do
  not treat a local output update as proof that the answer files are correct.
- Check API docs for accurate parameters, returns, defaults, ownership,
  constraints, and failure behavior. Keep notation and callable signatures
  consistent across source, headers, and all parallel documentation pages.
- When defaults change, update every corresponding table and older/newer
  module variant. User-visible API additions, behavior changes, and
  deprecations need the appropriate `versionadded`, `versionchanged`, or
  `deprecated` directive with placeholder version `x.y.z`, plus the project’s
  changelog and `doc/shared/RecentChanges.rst` updates when applicable.
- In documentation, use the heading hierarchy and link, footnote, directive,
  and organization rules from the documentation guides. Document prohibited
  build configurations and backend/precision limitations in installation
  guidance when relevant.

## Findings and reviewer voice

Report only actionable findings. A finding should identify a concrete changed
line, explain the failure mode or maintenance cost, and give the smallest
reasonable fix. Rank findings by impact:

- Blocking: correctness, memory/resource safety, broken API contract, lost
  error propagation, or a build/portability regression.
- Important: missing tests for changed behavior, incomplete public
  documentation, compatibility risk, stale state, or a maintainability issue
  likely to cause future defects.
- Nit: an explicit developer-guide violation with no meaningful functional
  impact. Do not manufacture nits.

Prefer high-confidence inline comments over a broad checklist. Use the
historically common question-led style: identify the exact behavior, ask
whether the intended contract is clear, explain the consequence, and suggest
a specific correction. Softening language is useful for uncertainty, but never
hide a blocker behind “maybe.” Use a code suggestion when the replacement is
small and unambiguous. Keep one issue per comment, cite related code or a
prior thread when useful, and do not ask the author to resolve a thread before
the replacement has been pushed.

Separate blocking fixes from follow-up improvements. It is acceptable to say
that a concern does not need to be addressed in the current PR while recording
the future risk. Do not imitate automated bot comments, cite historical
reviewer names as authority, or infer a defect from keyword matches alone.

Each comment must include:

- the affected file and line
- the relevant SUNDIALS rule or contract, with its guide file when applicable
- the impact on behavior, portability, bindings, documentation, or maintenance
- a concise suggested fix or a focused question that resolves the uncertainty

If no actionable findings remain, say so and list meaningful review gaps such
as unrun tests, formatters, platform builds, or documentation builds.
