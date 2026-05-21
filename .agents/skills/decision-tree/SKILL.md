---
name: decision-tree
description: Guide SUNDIALS users to the right package, time-stepping module, nonlinear/linear solver strategy, and initial tuning choices. Use when a request asks which SUNDIALS package to use (CVODE, CVODES, IDA, IDAS, ARKODE, KINSOL), whether to choose explicit vs implicit vs IMEX or multirate methods, how to pick dense/band/sparse/Krylov solvers and preconditioners, or what tolerances, step limits, initial-condition handling, and other starting settings to try.
---

# Choose a SUNDIALS Strategy

## Overview

Map a user problem to the right SUNDIALS package and give a defensible first-pass configuration. Keep the recommendation concrete: name the package, method family, linear/nonlinear solver approach, starting settings, and the next doc/example file to inspect.

Start by extracting only the facts that actually control package and solver choice. If key facts are missing, ask for the minimum needed to classify the problem: equation type, stiffness, sensitivities, natural splitting, system size, sparsity, and any special structure. Do not treat the `CVODE` versus `ARKODE` choice as purely a question of model structure; numerical stability and damping requirements can decide it even for an unsplit ODE.

## Intake Checklist

Collect these before recommending a package or method:

- problem form: nonlinear algebraic system, ODE IVP, or DAE IVP
- need for forward or adjoint sensitivities
- expected stiffness: nonstiff, stiff, or mixed stiff/nonstiff
- natural model split: explicit/implicit, slow/fast, Hamiltonian, operator splitting
- whether one-step RK stability properties matter: e.g., need for A-stability or L-stability at order above 2, strong damping of fast modes, or concern about higher-order BDF stability-angle limits
- size and matrix structure: small dense, banded, sparse, or very large matrix-free
- Jacobian and preconditioner availability
- important constraints: positivity, algebraic variables, event/rootfinding, fixed output cadence, structure preservation
- target hardware and vector/backend constraints if they matter to the recommendation

If the user does not know whether the model is stiff, infer it from context but label the assumption explicitly.

## Workflow

### 1. Choose the package family

- Use `KINSOL` for nonlinear algebraic systems `F(u) = 0` with no time integration.
- Use `IDA` or `IDAS` for DAEs of the form `F(t, y, y') = 0`.
- Use `CVODE` or `CVODES` for general ODE IVPs when a multistep solver is a good fit and there is no important reason to prefer one-step RK stability or splitting structure.
- Use `ARKODE` when the problem benefits from one-step Runge-Kutta structure or stability properties: explicit-only, implicit-only, IMEX/additive splitting, multirate evolution, low-storage explicit stepping, Hamiltonian structure, operator splitting, or a stiff problem where DIRK stability and damping are preferable to higher-order BDF behavior.
- If `CVODE` or `IDA` is the right choice, prefer the sensitivity-enabled supersets `CVODES` and `IDAS` when the user needs forward or adjoint sensitivities.

Open [package-selection.md](references/package-selection.md) when the package choice is the main question.

### 2. Choose the method or stepper

- For `CVODE` or `CVODES`, recommend `CV_ADAMS` for nonstiff problems and `CV_BDF` for stiff problems when a multistep BDF method is still a good numerical fit.
- For `IDA` or `IDAS`, the main integration method is variable-order BDF; focus the decision on initial-condition handling, Jacobian strategy, and linear solver choice.
- For `ARKODE`, select the narrowest stepper that matches the model:
  - `ERKStep` for fully explicit ODEs.
  - `ARKStep` for fully implicit DIRK or IMEX/additive explicit-implicit splits, especially when one-step A-stable or L-stable behavior is more important than the variable-order BDF workflow.
  - `MRIStep` for genuine slow/fast multirate problems.
  - `SPRKStep` for separable Hamiltonian systems where structure preservation matters.
  - `LSRKStep`, `ForcingStep`, or `SplittingStep` only when the user explicitly benefits from those formulations.
- For `KINSOL`, choose `KIN_LINESEARCH` when robustness matters, `KIN_NONE` for plain Newton with a good initial guess, `KIN_FP` for natural fixed-point maps, and `KIN_PICARD` when a Picard iteration is natural and a linear solver is available.
- When comparing `CV_BDF` and implicit `ARKStep`, remember that BDF methods above order 2 are not A-stable. If the user needs higher-order stiff decay, stronger damping of parasitic fast modes, or a parabolic diffusion problem is behaving poorly under high-order BDF, either cap BDF order at 2 or move the recommendation toward an appropriate DIRK table in `ARKStep`.

Open [method-and-linear-solvers.md](references/method-and-linear-solvers.md) when the package is already known but the method or solver stack is not.

### 3. Choose the linear and nonlinear solver strategy

- Use dense direct solvers for small dense systems.
- Use band direct solvers when the Jacobian bandwidth is small and known.
- Use sparse direct solvers such as KLU when the Jacobian is sparse and factorization cost is still acceptable.
- Use Krylov solvers plus preconditioning for large stiff systems or matrix-free settings.
- Recommend GMRES first when the user needs a generic Krylov choice.
- Recommend FGMRES when the preconditioner changes between iterations.
- Recommend PCG only for symmetric positive definite linear systems.
- Mention BiCGStab or TFQMR when storage is tighter than GMRES or when GMRES restart behavior is a concern.
- If the recommendation depends on a good preconditioner and the user has none, say that clearly instead of overselling Krylov methods.

### 4. Choose initial settings

- Recommend scalar tolerances only when component scales are comparable; otherwise prefer vector absolute tolerances.
- Tie `atol` to the smallest meaningful magnitude per component, not to machine precision by default.
- Tie `rtol` to the requested relative accuracy; if the user has no target, recommend a moderate starting value and tell them to run a convergence study.
- Suggest `IDASetId` and `IDACalcIC` for semi-explicit index-one DAEs when consistent initial conditions are not already available.
- Suggest user Jacobians or Jacobian-vector routines when the problem is sparse, banded, expensive, or noisy under finite differences.
- Suggest `SetMaxStep` only when there is a known physical or event scale that should cap the step size.
- If the solver hits the default internal-step limit, point out that the package default is often `500` steps before the next output time and explain whether the fix is larger limits, different tolerances, or a different method.

Open [settings-checklist.md](references/settings-checklist.md) when the user mainly needs tolerances, IC handling, step limits, Jacobian/preconditioner guidance, or failure triage.

## Output Style

When giving a recommendation:

- name the chosen package and why alternatives were rejected
- name the method or stepper and whether the problem is being treated as stiff, nonstiff, mixed, multirate, or structure-preserving
- name the linear/nonlinear solver strategy
- give a short starting settings block
- point to the closest repo docs and examples to inspect next

Prefer wording like:

- "Use `IDAS` because the model is a DAE and you also need sensitivities."
- "Start with `CV_BDF + Newton + GMRES` because the system appears stiff and too large for dense factorization."
- "Use `ARKStep` rather than `CVODE` because the RHS already has a meaningful stiff/nonstiff split."
- "Use `ARKStep` rather than `CVODE` because the problem is stiff and the one-step DIRK stability properties matter, even though there is no natural IMEX split."

Ground recommendations in these repo docs when needed:

- `doc/cvode/guide/source/Introduction.rst`
- `doc/cvode/guide/source/Usage/index.rst`
- `doc/cvodes/guide/source/Introduction.rst`
- `doc/ida/guide/source/Introduction.rst`
- `doc/ida/guide/source/Mathematics.rst`
- `doc/idas/guide/source/Introduction.rst`
- `doc/idas/guide/source/Usage/SIM.rst`
- `doc/arkode/guide/source/Introduction.rst`
- `doc/arkode/guide/source/Usage/index.rst`
- `doc/kinsol/guide/source/Introduction.rst`
- `doc/kinsol/guide/source/Usage/index.rst`

When a user wants code, inspect the closest example under `examples/<package>/` and adapt that pattern instead of describing an abstract setup.

## Pitfalls

Just because a SUNDIALS example uses a method/solver/setting doesn't mean its the right choice. Decisions should be grounded in doc recommendations and published literature on time integrator and solver methods. In particular, do not reduce `CVODE` versus `ARKODE` to "no split" versus "has split": stability region, stiff decay, stage order, and order restrictions also matter.
