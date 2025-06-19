# ARKODE
### Version 6.4.0 (Jun 2025)

**Daniel R. Reynolds,
  Department of Mathematics, SMU**

**David J. Gardner, Carol S. Woodward, and Cody J. Balos,
  Center for Applied Scientific Computing, LLNL**

ARKODE is a package for the solution of stiff, nonstiff, and multirate ordinary
differential equation (ODE) systems (initial value problems) given in linearly
implicit the form
```
M y' = f1(t,y) + f2(t,y), y(t0) = y0.
```
The integration methods implemented in ARKODE include explicit and implicit
Runge-Kutta methods, implicit-explicit (IMEX) additive Runge-Kutta methods, and
multirate infinitesimal (MRI) methods.

ARKODE is part of a the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKODE, CVODE, CVODES, IDA, IDAS, and KINSOL.
It is written in ANSI standard C and can be used in a variety of computing
environments including serial, shared memory, distributed memory, and
accelerator-based (e.g., GPU) systems. This flexibility is obtained from a
modular design that leverages the shared vector, matrix, linear solver, and
nonlinear solver APIs used across SUNDIALS packages.

## Documentation

See the ARKODE documentation at [Read the Docs](https://sundials.readthedocs.io/en/latest/arkode)
for more information about ARKODE usage.

## Installation

For installation instructions see the
[SUNDIALS Installation Guide](https://sundials.readthedocs.io/en/latest/Install_link.html).

## Release History

Information on recent changes to ARKODE can be found in the "Introduction"
chapter of the ARKODE User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the ARKODE User Guide.

## References

* D. R. Reynolds, D. J. Gardner, C. S. Woodward, and C. J. Balos,
  "User Documentation for ARKODE v6.4.0," LLNL technical report
  LLNL-SM-668082, Jun 2025.

* D. R. Reynolds, "Example Programs for ARKODE v6.4.0," Technical Report,
  Southern Methodist University Center for Scientific Computation, Jun 2025.
