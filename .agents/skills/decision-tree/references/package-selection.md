# Package Selection

## Contents

- Quick map
- ODE packages
- DAE packages
- Nonlinear system package
- ARKODE stepper map
- Notes that matter in practice

## Quick map

Choose the package from the mathematical problem first.

| Problem shape | Choose | Why |
| --- | --- | --- |
| `F(u) = 0` nonlinear algebraic system | `KINSOL` | No time integration; focus is nonlinear solve strategy |
| `y' = f(t, y)` general ODE IVP | `CVODE` | General-purpose variable-step multistep ODE solver |
| `y' = f(t, y)` with sensitivities | `CVODES` | `CVODE` plus forward and adjoint sensitivities |
| `F(t, y, y') = 0` DAE IVP | `IDA` | Variable-step BDF DAE solver |
| `F(t, y, y') = 0` with sensitivities | `IDAS` | `IDA` plus forward and adjoint sensitivities |
| ODE with explicit/implicit split, multirate structure, Hamiltonian structure, or strong one-step RK preference | `ARKODE` | Exposes stepper families matched to problem structure |

## ODE packages

Use `CVODE` or `CVODES` when the model is a standard ODE IVP and there is no strong reason to exploit a special split or structure.

Prefer `CVODES` over `CVODE` when the request includes:

- forward sensitivities with respect to parameters
- adjoint sensitivities or gradient calculations
- quadratures coupled to sensitivity workflows

Use `ARKODE` instead of `CVODE(S)` when one of these is true:

- the RHS is naturally split into nonstiff and stiff parts and the split is worth exploiting
- the problem is genuinely multirate
- the user wants a one-step Runge-Kutta method instead of a multistep method
- the problem is Hamiltonian or otherwise structure-preserving integration matters
- low-storage explicit RK or super-time-stepping is a main requirement

## DAE packages

Use `IDA` or `IDAS` for DAEs written as `F(t, y, y') = 0`.

Prefer `IDAS` when sensitivities are needed. `IDAS` is a superset of `IDA`.

Practical cues that point to `IDA(S)`:

- algebraic variables are present
- the residual depends on both `y` and `y'`
- consistent initial conditions are a real issue
- the user mentions an index-one semi-explicit DAE

## Nonlinear system package

Use `KINSOL` when the task is to solve a steady-state or nonlinear algebraic system, not to integrate in time.

Typical triggers:

- "solve a nonlinear system"
- "find the steady state"
- "Newton-Krylov"
- "Picard iteration"
- "fixed-point solve"

## ARKODE stepper map

Choose the narrowest ARKODE stepper that matches the model.

| Stepper | Use when |
| --- | --- |
| `ERKStep` | The ODE is fully explicit and nonstiff enough for explicit RK |
| `ARKStep` | The ODE is fully implicit DIRK or has an explicit/implicit additive split for IMEX |
| `MRIStep` | The problem has true slow/fast multirate structure |
| `SPRKStep` | The system is separable Hamiltonian and structure preservation matters |
| `LSRKStep` | Low-storage explicit RK, SSP, or super-time-stepping is the main goal |
| `SplittingStep` | The problem is posed around operator splitting and that formulation is intentional |
| `ForcingStep` | The user specifically needs that forcing-based stepping formulation |

If the user just says "I have an ODE, which SUNDIALS package should I use?", do not jump to niche ARKODE steppers. Start with `CVODE(S)` or the main ARKODE steppers only when the model structure clearly warrants it.

## Notes That Matter In Practice

- Do not recommend both `CVODE` and `CVODES` together for one application; `CVODES` already covers the `CVODE` functionality.
- Do not recommend both `IDA` and `IDAS` together for one application; `IDAS` is the sensitivity-enabled superset.
- For DAEs, package selection is usually easier than consistent-IC setup. If the problem is index-one and the initial conditions are not consistent, plan to discuss `IDASetId` and `IDACalcIC`.
- For large stiff systems, the package choice and the linear solver choice are tightly coupled. If the user cannot supply a useful preconditioner, say that solver performance may be limited.
