#!/usr/bin/env python3
# --------------------------------------------------------------------------------------
# Programmer(s): David J. Gardner
# --------------------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------
# Example: Lotka-Volterra Predator-Prey ODE System
#
# This example demonstrates how to use CVODES to solve a simple ODE system. The
# Lotka-Volterra equations model predator-prey dynamics in ecology:
#
#    du/dt =  p_0*u - p_1*u*v    (prey growth minus predation)
#    dv/dt = -p_2*v + p_3*u*v    (predator death plus growth from eating prey)
#
# where:
#   u = prey population
#   v = predator population
#   p_0 = prey birth rate
#   p_1 = predation rate
#   p_2 = predator death rate
#   p_3 = predator growth rate from predation
#
# We use the parameters p = [1.5, 1.0, 3.0, 1.0], initial condition y(0) = [1.0, 1.0],
# and integration interval t in [0, 10].
# --------------------------------------------------------------------------------------
# !!! ATTENTION SUNDIALS DEVELOPERS: This example is imported into the Python binding
# documentation to illustrate sundials4py usage. When updating this file, make sure to
# update documentation as needed e.g., line highlight numbers. !!!
# --------------------------------------------------------------------------------------

# --- start example ---
import numpy as np
import sys
import matplotlib.pyplot as plt
from sundials4py.core import *  # Always import the core submodule
from sundials4py.cvodes import *  # Import the desired SUNDIALS package


class LotkaVolterraODE:
    """
    Encapsulates the Lotka-Volterra ODE problem.

    This class defines the ODE system and provides the functions that CVODE needs to
    evolve it in time: the right-hand side (RHS) function and the Jacobian.
    """

    def __init__(self, p):
        """
        Initialize with model parameters.

        Args:
            p: list or array of 4 parameters [p_0, p_1, p_2, p_3]
        """
        self.p = np.array(p, dtype=sunrealtype)
        self.NEQ = 2  # Number of equations in the system

    def set_init_cond(self, yvec):
        """
        Set the initial conditions in the solution vector.

        Args:
            yvec: N_Vector to store initial values

        Returns:
            0 on success
        """
        y = N_VGetArrayPointer(yvec)  # Returns a numpy ndarray view of the data
        y[0] = 1.0  # Initial prey population
        y[1] = 1.0  # Initial predator population
        return 0

    def rhs(self, t, yvec, ydotvec, _):
        """
        Right-hand side function: computes dy/dt = f(t, y).

        Args:
            t: current time (not used in this autonomous system)
            yvec: N_Vector with the current solution values y = [u, v]
            ydotvec: N_Vector output vector for time derivatives f = [du/dt, dv/dt]
            _: user data pointer MUST NOT be used in Python

        Returns:
            0 on success
            positive value for a recoverable error (integrator will retry the step)
            negative value for an unrecoverable error (integration will halt)
        """
        p = self.p
        y = N_VGetArrayPointer(yvec)  # Returns a numpy ndarray view of the data
        ydot = N_VGetArrayPointer(ydotvec)

        # Compute the derivatives
        ydot[0] = p[0] * y[0] - p[1] * y[0] * y[1]  # du/dt: prey dynamics
        ydot[1] = -p[2] * y[1] + p[3] * y[0] * y[1]  # dv/dt: predator dynamics
        return 0

    def jac(self, t, yvec, fyvec, J, _, tmp1, tmp2, tmp3):
        """
        Jacobian function: computes J = df/dy.

        Args:
            t: current time
            yvec: N_Vector with the current solution values
            fyvec: N_Vector with the current RHS values (not used here)
            J: SUNmatrix to store the output Jacobian matrix
            _: user data pointer MUST NOT be used in Python
            tmp1, tmp2, tmp3: temporary workspace N_Vectors (not used here)

        Returns:
            0 on success
            positive value for a recoverable error (integrator will retry the step)
            negative value for an unrecoverable error (integration will halt)
        """
        p = self.p
        y = N_VGetArrayPointer(yvec)  # Returns a numpy ndarray view of the data
        Jdata = SUNDenseMatrix_Data(J)  # Returns a numpy ndarray view of the data

        # Compute partial derivatives
        # J[0,0] = d(du/dt)/du, J[0,1] = d(du/dt)/dv
        Jdata[0, 0] = p[0] - p[1] * y[1]
        Jdata[0, 1] = -p[1] * y[0]

        # J[1,0] = d(dv/dt)/du, J[1,1] = d(dv/dt)/dv
        Jdata[1, 0] = p[3] * y[1]
        Jdata[1, 1] = -p[2] + p[3] * y[0]
        return 0


def main():
    """
    Main function demonstrating the SUNDIALS/CVODE workflow.

    A typical workflow for solving an ODE with CVODE is:
    1. Create the SUNDIALS context
    2. Create the solution vector and set initial conditions
    3. Create and initialize CVODE
    4. Configure CVODE (set tolerances, linear solver, Jacobian, etc.)
    5. Advance the solution in time
    6. Retrieve CVODE statistics
    7. Destroy objects and free memory (handled automatically in Python)
    """

    # ----------------------------------------------------------------------------------
    # Step 1: Create the SUNDIALS context
    # ----------------------------------------------------------------------------------

    # The context provides support for features such as error handling, logging, and
    # profiling and is required to construct all SUNDIALS objects. SUN_COMM_NULL is used
    # for serial (non-parallel) problems.

    # In C, sunctx is an output parameter (return-by-pointer):
    # SUNContext_Create(SUN_COMM_NULL, &sunctx). In Python, it is returned in a tuple
    # where the first entry is the function return value and the second is the context.
    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    # ----------------------------------------------------------------------------------
    # Step 2: Set up the ODE problem
    # ----------------------------------------------------------------------------------

    # Create the ODE problem with parameters p = [1.5, 1.0, 3.0, 1.0]
    params = [1.5, 1.0, 3.0, 1.0]
    ode = LotkaVolterraODE(params)

    # Create a serial N_Vector to hold the solution
    y = N_VNew_Serial(ode.NEQ, sunctx)
    assert y is not None

    # Set initial conditions: y(0) = [1.0, 1.0]
    ode.set_init_cond(y)

    # ----------------------------------------------------------------------------------
    # Step 3: Create and initialize CVODE
    # ----------------------------------------------------------------------------------

    # Create a CVODE instance using Backward Differentiation Formulas (BDF)
    cvode = CVodeCreate(CV_BDF, sunctx)
    assert cvode is not None

    # Initialize CVODE with the RHS function, initial time (t0=0.0), and initial state
    status = CVodeInit(cvode.get(), ode.rhs, 0.0, y)  # access CVODE memory with .get()
    assert status == CV_SUCCESS

    # ----------------------------------------------------------------------------------
    # Step 4: Set CVODE options (tolerances, linear solver, etc.)
    # ----------------------------------------------------------------------------------

    # CVODE will adapt the step size and method order so the estimated local error meets
    # the user defined tolerances.
    reltol = 1e-6  # Relative tolerance
    abstol = 1e-10  # Absolute tolerance
    status = CVodeSStolerances(cvode.get(), reltol, abstol)
    assert status == CV_SUCCESS

    # CVODE defaults to using a Newton iteration to solve nonlinear systems at each time
    # step which requires solving a linear system each iteration. For small problems, a
    # dense direct solver is efficient and easy to use.

    # Create a dense matrix (NEQ x NEQ) to store the Jacobian
    A = SUNDenseMatrix(ode.NEQ, ode.NEQ, sunctx)
    assert A is not None

    # Create a dense linear solver that works with the matrix A and vector y
    LS = SUNLinSol_Dense(y, A, sunctx)
    assert LS is not None

    # Attach the linear solver to CVODE
    status = CVodeSetLinearSolver(cvode.get(), LS, A)
    assert status == CV_SUCCESS

    # Providing an analytical Jacobian can significantly improve performance. If not
    # provided, CVODE will approximate it using finite differences.
    status = CVodeSetJacFn(cvode.get(), ode.jac)
    assert status == CV_SUCCESS

    # ----------------------------------------------------------------------------------
    # Step 5: Advance the ODE in time
    # ----------------------------------------------------------------------------------

    # Set the current output time and get access to the underlying array data for output
    tret = 0.0
    yarr = N_VGetArrayPointer(y)  # Returns a numpy ndarray view of the data

    # Print header with problem information
    print("\nLotka-Volterra ODE (CVODE):")
    print(f"    initial condition, y0 = [1.0, 1.0]")
    print(f"    parameters = {params}")
    print(f"    reltol = {reltol}, abstol = {abstol}\n")
    print("        t           u           v")
    print("   ---------------------------------")
    print(f"  {tret:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")

    # Lists to store solution for plotting
    t_vals = [tret]
    u_vals = [yarr[0]]
    v_vals = [yarr[1]]

    # Integrate from t=0 to t=10, outputting every 0.05 time units using CV_NORMAL mode.
    # In normal mode, CVODE will take internal steps until it has reached or just passed
    # the output time and then return a time interpolated solution at the output time.
    dtout = 0.05
    iout = 0
    while tret < 10.0:
        # Advance the system and return the solution at requested output time
        status, tret = CVode(cvode.get(), tret + dtout, y, CV_NORMAL)
        assert status == CV_SUCCESS

        # Store solution for plotting (every output)
        t_vals.append(tret)
        u_vals.append(yarr[0])
        v_vals.append(yarr[1])

        # Print every 10th output (every 1.0 time unit) to keep the output readable
        iout += 1
        if iout % 10 == 0 or tret >= 10.0:
            print(f"  {tret:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")
    print("   ---------------------------------")

    # ----------------------------------------------------------------------------------
    # Step 6: Get CVODE statistics
    # ----------------------------------------------------------------------------------

    # CVODE provides various statistics about the integration that can help assess
    # performance and diagnose issues. These include values such as the step counts,
    # function evaluations, and information about the linear solver.

    status, nst = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS

    status, netf = CVodeGetNumErrTestFails(cvode.get())
    assert status == CV_SUCCESS

    status, nfe = CVodeGetNumRhsEvals(cvode.get())
    assert status == CV_SUCCESS

    status, nsetups = CVodeGetNumLinSolvSetups(cvode.get())
    assert status == CV_SUCCESS

    status, nje = CVodeGetNumJacEvals(cvode.get())
    assert status == CV_SUCCESS

    # For C functions with multiple output parameters (return-by-pointer),
    # CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn), the first element is always the
    # function return value (error code) and the subsequent elements follow the same
    # ordering as in the C function signature.
    status, nni, ncfn = CVodeGetNonlinSolvStats(cvode.get())
    assert status == CV_SUCCESS

    print("\nIntegrator Statistics:")
    print(f"  Number of steps taken                    = {nst}")
    print(f"  Number of error test fails               = {netf}")
    print(f"  Number of RHS evaluations                = {nfe}")
    print(f"  Number of linear solver setups           = {nsetups}")
    print(f"  Number of Jacobian evaluations           = {nje}")
    print(f"  Number of nonlinear solver iterations    = {nni}")
    print(f"  Number of nonlinear convergence failures = {ncfn}")

    # ----------------------------------------------------------------------------------
    # Plot the solution
    # ----------------------------------------------------------------------------------

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left plot: Time series of prey and predator populations
    ax1.plot(t_vals, u_vals, "b-", label="Prey (u)")
    ax1.plot(t_vals, v_vals, "r--", label="Predator (v)")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Population")
    ax1.set_title("Lotka-Volterra: Population vs Time")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Right plot: Phase portrait (predator vs prey)
    ax2.plot(u_vals, v_vals, "g-")
    ax2.plot(u_vals[0], v_vals[0], "ko", label="Start")
    ax2.plot(u_vals[-1], v_vals[-1], "ks", label="End")
    ax2.set_xlabel("Prey (u)")
    ax2.set_ylabel("Predator (v)")
    ax2.set_title("Phase Portrait")
    ax2.legend()
    ax2.grid(alpha=0.3)

    # Display the plot if not running in test mode
    if "pytest" not in sys.modules:
        plt.tight_layout()
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
# --- end example ---


def test_cvs_lotkavolterra():
    """Function for pytest to discover the example as a test"""
    main()
