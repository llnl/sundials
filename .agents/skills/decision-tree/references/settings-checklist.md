# Settings Checklist

## Contents

- Tolerances
- Step-size controls
- DAE-specific setup
- KINSOL-specific setup
- Diagnostics and triage
- Docs and examples to inspect next

## Tolerances

Base tolerances on the physics and the requested accuracy, not on arbitrary tiny numbers.

For time integrators:

- use scalar absolute tolerance only when all state components live on similar scales
- use vector absolute tolerances when component magnitudes differ materially
- set each `atol` near the smallest meaningful magnitude for that component
- set `rtol` from the required relative accuracy

The SUNDIALS guides explicitly recommend experimenting with tolerances and checking how the solution changes as tolerances are reduced. If the user knows the target global error, use the package-specific rules of thumb from the guides as a starting point:

- `CVODE`, `CVODES`, `IDA`, and `IDAS`: set local tolerances roughly two orders of magnitude tighter than the desired global error
- `ARKODE`: set local tolerances roughly one order of magnitude tighter than the desired global error

Treat those as starting guesses, not guarantees, and tell the user to confirm them with a convergence check.

If the user has no accuracy target:

- recommend a moderate first pass, not a heroic one
- tell them to verify the answer by tightening tolerances and comparing outputs

## Step-Size Controls

Do not set every knob by default.

Use these only when the problem justifies them:

- `SetMaxStep`: cap the step size when there is a known physical, event, or output scale that should not be skipped
- `SetMinStep`: use sparingly; it can turn recoverable adaptivity problems into hard failures
- `SetMaxNumSteps`: raise this when the solver legitimately needs more than the package default before the next output time
- fixed-step options in ARKODE: use only when the user truly wants fixed-step integration or is matching an external scheme

Useful default facts from the package guides:

- `CVODE` and `IDA(S)` default to a maximum of `500` internal steps before the next output time

If that limit is hit, explain the possible causes:

- output times are too far apart
- tolerances are tighter than necessary
- the problem is stiffer than assumed
- the current method or package is a poor fit

## DAE-Specific Setup

For `IDA` and `IDAS`, ask early about consistent initial conditions.

If the model is a semi-explicit index-one DAE and the initial conditions are only approximate:

- identify differential vs algebraic components with `IDASetId`
- use `IDACalcIC` to compute consistent initial conditions when appropriate

Also consider:

- `IDASetConstraints` when inequality constraints such as positivity matter
- user Jacobians or Jacobian-vector routines when the residual is expensive or structurally sparse

## KINSOL-Specific Setup

For `KINSOL`, tune only the settings that clearly matter.

Commonly relevant knobs:

- `KINSetFuncNormTol` for the function-norm stopping target
- `KINSetScaledStepTol` for step-length stopping
- `KINSetMaxNewtonStep` when uncontrolled Newton steps are a real risk
- `KINSetNumMaxIters` when the default iteration budget is clearly wrong
- Anderson acceleration settings such as `KINSetMAA`, `KINSetDampingAA`, and `KINSetDelayAA` when using `KIN_FP` or `KIN_PICARD`

Practical recommendation:

- start with `KIN_LINESEARCH` before changing many tolerances
- only push on Anderson acceleration after confirming that a fixed-point or Picard formulation is actually appropriate

## Diagnostics and Triage

Use solver statistics and failure modes to refine the recommendation.

If you see many error-test failures:

- check whether tolerances are unrealistic
- check whether `atol` matches state scales
- check whether the model is less smooth than assumed

If you see nonlinear convergence failures:

- move toward a better Jacobian or preconditioner
- use a more robust nonlinear strategy
- revisit whether the problem should be treated as stiff

If you see many linear iterations:

- improve the preconditioner
- change the Krylov method only after checking preconditioning

If direct solves are too expensive:

- move from dense or sparse direct factorization to Krylov plus preconditioning

If the user is fighting consistent-IC failures in a DAE:

- verify the DAE classification
- verify the differential/algebraic partition
- use `IDACalcIC` only when its assumptions fit the problem

## Docs and Examples To Inspect Next

Point the user to the nearest package guide and example after making a recommendation.

Good starting points in this repo:

- `doc/cvode/guide/source/Usage/index.rst`
- `doc/cvodes/guide/source/Usage/SIM.rst`
- `doc/ida/guide/source/Mathematics.rst`
- `doc/idas/guide/source/Usage/SIM.rst`
- `doc/arkode/guide/source/Usage/index.rst`
- `doc/kinsol/guide/source/Usage/index.rst`
- `examples/cvode/`
- `examples/cvodes/`
- `examples/ida/`
- `examples/idas/`
- `examples/arkode/`
- `examples/kinsol/`
