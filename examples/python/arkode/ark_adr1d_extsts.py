#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Python version of the C++ example
#   examples/arkode/CXX_serial/ark_adr1d_extsts.cpp
#
# This example solves a one-dimensional advection-diffusion-reaction
# system with three chemical species,
#
#   u_t = -c u_x + d u_xx + A - (w + 1) * u + v * u^2
#   v_t = -c v_x + d v_xx + w * u - v * u^2
#   w_t = -c w_x + d w_xx + (B - w) / eps - w * u
#
# on x in [0, 1], with stationary boundary conditions. The time
# integration uses an ExtSTS method:
#
#   diffusion: explicit inner STS method
#   advection: explicit slow RHS
#   reaction : implicit slow RHS
#
# The spatial derivatives use second-order centered differences on a
# uniform mesh.  The example stores the state as a serial N_Vector, but
# the callbacks view that vector as an (nx, 3) NumPy array so that each
# row is one mesh point and the columns are the u, v, and w species.
# -----------------------------------------------------------------

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np
import sundials4py.arkode as ark
import sundials4py.core as sun


class ADR1DExtSTSProblem:
    """Advection-diffusion-reaction problem data and ARKODE callbacks."""

    def __init__(
        self, c=1.0e-2, d=1.0e-1, A=0.6, B=2.0, eps=1.0e-2, tf=3.0, xl=0.0, xu=1.0, nx=512
    ):
        self.c = c
        self.d = d
        self.A = A
        self.B = B
        self.eps = eps
        self.tf = tf
        self.xl = xl
        self.xu = xu
        self.nx = nx

        self.dx = (xu - xl) / (nx - 1)
        self.neq = 3 * nx

    def view2d(self, yvec):
        """Return a zero-copy (mesh point, species) view of an N_Vector."""
        return sun.N_VGetArrayPointer(yvec).reshape((self.nx, 3))

    def set_initial_condition(self, yvec):
        y = self.view2d(yvec)
        x = np.linspace(self.xl, self.xu, self.nx)
        perturbation = 0.1 * np.sin(np.pi * x)

        y[:, 0] = self.A + perturbation
        y[:, 1] = self.B / self.A + perturbation
        y[:, 2] = self.B + perturbation
        return 0

    def f_advection(self, t, yvec, fvec, user_data):
        y = self.view2d(yvec)
        f = self.view2d(fvec)

        # Centered first derivative for the advection term.  The boundary
        # rows are left at zero to impose stationary boundary conditions.
        f[0, :] = 0.0
        f[-1, :] = 0.0
        f[1:-1, :] = -self.c * (y[2:, :] - y[:-2, :]) / (2.0 * self.dx)
        return 0

    def f_diffusion(self, t, yvec, fvec, user_data):
        y = self.view2d(yvec)
        f = self.view2d(fvec)

        # Centered second derivative for diffusion.  This is the fast RHS
        # evolved by the inner STS method.
        f[0, :] = 0.0
        f[-1, :] = 0.0
        f[1:-1, :] = self.d * (y[:-2, :] - 2.0 * y[1:-1, :] + y[2:, :]) / (self.dx * self.dx)
        return 0

    def f_reaction(self, t, yvec, fvec, user_data):
        y = self.view2d(yvec)
        f = self.view2d(fvec)

        # Reaction terms are local to each mesh point and are treated
        # implicitly by the outer ExtSTS method.
        f[0, :] = 0.0
        f[-1, :] = 0.0
        u = y[1:-1, 0]
        v = y[1:-1, 1]
        w = y[1:-1, 2]

        f[1:-1, 0] = self.A - (w + 1.0) * u + v * u * u
        f[1:-1, 1] = w * u - v * u * u
        f[1:-1, 2] = (self.B - w) / self.eps - w * u
        return 0

    def diffusion_domeig(self, t, yvec, fnvec, user_data, temp1, temp2, temp3):
        # Dominant eigenvalue estimate for the finite-difference diffusion
        # operator. The STS method uses this to choose stable inner steps.
        lambdaR = -4.0 * self.d / (self.dx * self.dx)
        lambdaI = 0.0
        return 0, lambdaR, lambdaI


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="1D advection-diffusion-reaction problem using ARKODE ExtSTS"
    )

    parser.add_argument("--c", type=float, default=1.0e-2, help="advection speed")
    parser.add_argument("--d", type=float, default=1.0e-1, help="diffusion constant")
    parser.add_argument("--A", type=float, default=0.6, help="reaction parameter A")
    parser.add_argument("--B", type=float, default=2.0, help="reaction parameter B")
    parser.add_argument("--eps", type=float, default=1.0e-2, help="stiffness parameter")
    parser.add_argument("--tf", type=float, default=3.0, help="final time")
    parser.add_argument("--xl", type=float, default=0.0, help="domain lower boundary")
    parser.add_argument("--xu", type=float, default=1.0, help="domain upper boundary")
    parser.add_argument("--nx", type=int, default=512, help="number of mesh points")

    parser.add_argument(
        "--sts_method", choices=("RKC", "RKL"), default="RKC", help="inner STS method"
    )
    parser.add_argument(
        "--mri_method", default="ARKODE_IMEX_MRI_GARK_GKC21", help="MRI coupling table name"
    )
    parser.add_argument("--rtol", type=float, default=1.0e-4, help="relative tolerance")
    parser.add_argument("--atol", type=float, default=1.0e-9, help="absolute tolerance")
    parser.add_argument("--fixed_h", type=float, default=0.0, help="fixed outer step size")
    parser.add_argument(
        "--maxsteps", type=int, default=10000, help="maximum internal steps between output times"
    )
    parser.add_argument("--nout", type=int, default=10, help="number of output times")

    args = parser.parse_args(argv)

    if args.nx < 3:
        parser.error("--nx must be at least 3")
    if args.xu <= args.xl:
        parser.error("--xu must be greater than --xl")
    if args.tf <= 0.0:
        parser.error("--tf must be positive")
    if args.eps <= 0.0:
        parser.error("--eps must be positive")
    if args.rtol <= 0.0 or args.atol <= 0.0:
        parser.error("--rtol and --atol must be positive")
    if args.fixed_h < 0.0:
        parser.error("--fixed_h must be nonnegative")
    if args.maxsteps <= 0:
        parser.error("--maxsteps must be positive")
    if args.nout <= 0:
        parser.error("--nout must be positive")

    return args


def print_setup(problem, args):
    print("\nProblem parameters and options:")
    print(" --------------------------------- ")
    print(f"  c                = {problem.c}")
    print(f"  d                = {problem.d}")
    print(f"  A                = {problem.A}")
    print(f"  B                = {problem.B}")
    print(f"  eps              = {problem.eps}")
    print(" --------------------------------- ")
    print(f"  tf               = {problem.tf}")
    print(f"  xl               = {problem.xl}")
    print(f"  xu               = {problem.xu}")
    print(f"  nx               = {problem.nx}")
    print(f"  dx               = {problem.dx}")
    print(" --------------------------------- ")
    print("  integrator       = ExtSTS")
    print("  advection        = Explicit")
    print("  diffusion        = Explicit STS")
    print("  reaction         = Implicit")
    print(f"  rtol             = {args.rtol}")
    print(f"  atol             = {args.atol}")
    print(f"  fixed h          = {args.fixed_h}")
    print(" --------------------------------- ")
    print(f"  MRI method       = {args.mri_method}")
    print(f"  STS method       = {args.sts_method}")
    print(" --------------------------------- \n")


def rms_norm(yvec, neq):
    return np.sqrt(sun.N_VDotProd(yvec, yvec) / neq)


def print_solution_header():
    print("          t                    ||y||_rms")
    print(" -----------------------------------------------")


def print_solution(t, yvec, neq):
    print(f" {t:22.15e} {rms_norm(yvec, neq):22.15e}")


def print_solution_footer():
    print(" -----------------------------------------------\n")


def stats_text(arkode_mem):
    with tempfile.NamedTemporaryFile(
        prefix="ark_adr1d_extsts_", suffix=".txt", delete=False
    ) as tmp:
        stats_path = Path(tmp.name)

    try:
        status, file_ptr = sun.SUNFileOpen(str(stats_path), "w")
        assert status == sun.SUN_SUCCESS
        status = ark.ARKodePrintAllStats(arkode_mem.get(), file_ptr, sun.SUN_OUTPUTFORMAT_TABLE)
        assert status == ark.ARK_SUCCESS

        file_ptr = None
        return stats_path.read_text()
    finally:
        stats_path.unlink(missing_ok=True)


def print_stats(arkode, sts):
    print("\nExtSTS Integrator:")
    print(stats_text(arkode), end="")

    print("\nInner STS Method:")
    print(stats_text(sts), end="")


def main(argv=None):
    args = parse_args(argv)

    if not hasattr(ark, "MRIStepCreateExtSTS") or not hasattr(ark, "MRIStepGetSTS"):
        raise RuntimeError("This example requires sundials4py with ARKODE ExtSTS bindings enabled")

    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    problem = ADR1DExtSTSProblem(
        c=args.c,
        d=args.d,
        A=args.A,
        B=args.B,
        eps=args.eps,
        tf=args.tf,
        xl=args.xl,
        xu=args.xu,
        nx=args.nx,
    )

    print_setup(problem, args)

    y = sun.N_VNew_Serial(problem.neq, sunctx)
    assert y is not None
    problem.set_initial_condition(y)

    # Create the ExtSTS stepper.  Diffusion is the inner STS RHS, while
    # advection and reaction are explicit and implicit slow RHS functions.
    arkode = ark.MRIStepCreateExtSTS(
        problem.f_diffusion, problem.f_advection, problem.f_reaction, 0.0, y, sunctx
    )
    assert arkode is not None

    # Retrieve and configure the inner STS stepper created by ExtSTS.
    status, sts = ark.MRIStepGetSTS(arkode.get())
    assert status == ark.ARK_SUCCESS
    assert sts is not None

    sts_method = ark.ARKODE_LSRK_RKC_2 if args.sts_method == "RKC" else ark.ARKODE_LSRK_RKL_2
    status = ark.LSRKStepSetSTSMethod(sts.get(), sts_method)
    assert status == ark.ARK_SUCCESS

    status = ark.LSRKStepSetDomEigFn(sts.get(), problem.diffusion_domeig)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSStolerances(arkode.get(), args.rtol, args.atol)
    assert status == ark.ARK_SUCCESS

    if args.fixed_h > 0.0:
        status = ark.ARKodeSetFixedStep(arkode.get(), args.fixed_h)
        assert status == ark.ARK_SUCCESS

    # The reaction RHS is implicit.  Providing a band matrix and band linear
    # solver is enough here; ARKODE will approximate the Jacobian internally.
    J = sun.SUNBandMatrix(problem.neq, 2, 2, sunctx)
    assert J is not None

    LS = sun.SUNLinSol_Band(y, J, sunctx)
    assert LS is not None

    status = ark.ARKodeSetLinearSolver(arkode.get(), LS, J)
    assert status == ark.ARK_SUCCESS

    coupling = ark.MRIStepCoupling_LoadTableByName(args.mri_method)
    assert coupling is not None

    status = ark.MRIStepSetCoupling(arkode.get(), coupling)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetMaxNumSteps(arkode.get(), args.maxsteps)
    assert status == ark.ARK_SUCCESS

    t = 0.0
    dtout = problem.tf / args.nout
    tout = dtout

    print_solution_header()
    print_solution(t, y, problem.neq)

    for _ in range(args.nout):
        status, t = ark.ARKodeEvolve(arkode.get(), tout, y, ark.ARK_NORMAL)
        if status < 0:
            raise RuntimeError("ARKodeEvolve failed with status {status}")

        print_solution(t, y, problem.neq)
        tout = min(tout + dtout, problem.tf)

    print_solution_footer()

    print("Final integrator statistics:")
    print_stats(arkode, sts)


if __name__ == "__main__":
    main(sys.argv[1:])
