# Method and Linear Solver Selection

## Contents

- Time integration method choice
- KINSOL nonlinear strategy choice
- Linear solver choice
- When to push for Jacobians or preconditioners

## Time Integration Method Choice

### CVODE and CVODES

`CVODE` and `CVODES` expose two multistep choices:

- `CV_ADAMS` for nonstiff problems
- `CV_BDF` for stiff problems

Use `CV_BDF` when you see any of these:

- severe step-size restrictions with explicit methods
- disparate time scales with rapidly damped modes
- a Jacobian or preconditioner is already part of the model workflow
- the user is already expecting Newton or Krylov solves

Use `CV_ADAMS` when the problem is plainly nonstiff and the user wants a simpler nonstiff ODE workflow.

For stiff `CV_BDF` runs, the default Newton iteration is the natural recommendation. For nonstiff `CV_ADAMS`, fixed-point iteration can be reasonable if the problem is mild enough, but do not force it when the user is already prepared to use Newton.

### IDA and IDAS

`IDA` and `IDAS` use variable-order BDF methods. The main decision is not the integration formula but the solve stack around it:

- how to obtain consistent initial conditions
- whether the Jacobian should be dense, banded, sparse, or matrix-free
- whether direct factorization is affordable
- whether a preconditioned Krylov method is needed

### ARKODE

Choose the stepper from the model structure:

- `ERKStep` for explicit RK
- `ARKStep` for implicit DIRK or IMEX additive RK
- `MRIStep` for multirate
- `SPRKStep` for separable Hamiltonian systems

Within `ARKStep`, decide which part belongs in the explicit vs implicit operator. Put the genuinely stiff portion in the implicit operator and keep the explicit side for the terms that would otherwise restrict stability.

For many users, the built-in defaults and adaptivity are acceptable starting points. Only reach for custom Butcher tables or custom adaptivity when the user has a concrete reason.

## KINSOL Nonlinear Strategy Choice

Choose `KINSOL` strategy from robustness needs and model form.

| Strategy | Use when | Notes |
| --- | --- | --- |
| `KIN_LINESEARCH` | Robust default Newton choice | Best first recommendation when convergence risk is nontrivial |
| `KIN_NONE` | Plain Newton is acceptable | Requires a good initial guess and reasonable local behavior |
| `KIN_FP` | A natural fixed-point map exists | No linear solver required; Anderson acceleration is available |
| `KIN_PICARD` | A natural Picard map exists and a linear solver is available | Also supports Anderson acceleration |

If the user is unsure, recommend `KIN_LINESEARCH` first for robustness, then simplify only if there is a reason.

## Linear Solver Choice

Pick the linear solver from matrix structure and problem size.

### Direct Solvers

Use a direct solver when factorization is affordable.

- dense direct: small dense Jacobians
- band direct: narrow-band Jacobians with known bandwidth
- sparse direct such as KLU: sparse Jacobians where fill-in is still manageable

Direct methods are usually the cleanest recommendation for small to moderate problems because they remove preconditioner design from the initial setup.

### Krylov Solvers

Use Krylov methods for large stiff problems, large sparse problems, or matrix-free workflows.

Recommended first choice:

- GMRES / SPGMR as the best generic default

Alternatives:

- FGMRES / SPFGMR when the preconditioner changes between iterations
- BiCGStab / SPBCG when storage matters more than GMRES flexibility
- TFQMR / SPTFQMR when GMRES storage or restart behavior is a concern
- PCG when the linear system is symmetric positive definite

If you recommend a Krylov solver for a hard stiff problem, also discuss the preconditioner. A Krylov recommendation without a preconditioner plan is often incomplete.

## When To Push For Jacobians Or Preconditioners

Push for a user Jacobian, Jacobian-vector routine, or preconditioner when:

- finite differences are too expensive
- the Jacobian is sparse or banded and the structure is known
- the problem is large and stiff
- convergence failures suggest the internal approximation is not good enough
- residual evaluations are noisy enough that finite-difference Jacobians will be poor

For direct solvers:

- user Jacobians matter most when sparsity or bandwidth is important

For Krylov solvers:

- preconditioning is usually the difference between a plausible plan and an impractical one

Mention package-specific preconditioners when they match the user setup:

- `CVBBDPRE` / `CVBANDPRE` patterns for CVODE-family problems
- `IDABBDPRE` for distributed DAE cases
- `ARKBBDPRE` / `ARKBANDPRE` for matching ARKODE configurations
- `KINBBDPRE` for parallel KINSOL problems when the built-in block-banded approximation fits
