# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
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

import numpy as np
from sundials4py.core import *
from .problem import ODEProblem


class AnalyticODE(ODEProblem):
    """
    * The following is a simple example problem with analytical
    * solution,
    *    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
    * for t in the interval [0.0, 10.0], with initial condition: y=0.
    *
    * The stiffness of the problem is directly proportional to the
    * value of "lambda".  The value of lambda should be negative to
    * result in a well-posed ODE; for values with magnitude larger
    * than 100 the problem becomes quite stiff.
    """

    def __init__(self, lamb=-10.0):
        self.lamb = lamb
        self.inner_stepper = None

    def f(self, t, yvec, ydotvec, _):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - self.lamb * np.arctan(t)
        return 0

    def fS(self, Ns, t, y, ydot, yS, ySdot, _, tmp1, tmp2):
        # Sensitivity RHS: df/dy * yS + df/dp (here, p = y0, so df/dp = 0)
        yarr = N_VGetArrayPointer(y)
        ySarr = N_VGetArrayPointer(yS[0])
        ySdotarr = N_VGetArrayPointer(ySdot[0])
        # df/dy = lambda
        lamb = self.lamb
        ySdotarr[0] = lamb * ySarr[0]
        return 0

    def fB(self, t, y, yB, yBdot, _):
        # Adjoint RHS: -df/dy^T * lambdaB
        yarr = N_VGetArrayPointer(y)
        yBarr = N_VGetArrayPointer(yB)
        yBdotarr = N_VGetArrayPointer(yBdot)
        lamb = self.lamb
        yBdotarr[0] = -lamb * yBarr[0]
        return 0

    def fBS(self, t, y, yS, yB, yBdot, _):
        # Adjoint RHS: -df/dy^T * lambdaB
        yarr = N_VGetArrayPointer(y)
        yBarr = N_VGetArrayPointer(yB)
        yBdotarr = N_VGetArrayPointer(yBdot)
        lamb = self.lamb
        yBdotarr[0] = -lamb * yBarr[0]
        return 0

    def dom_eig(self, t, yvec, fnvec, tempv1, tempv2, tempv3):
        lamdbaR = self.lamb
        lamdbaI = 0.0
        return 0, lamdbaR, lamdbaI

    def jac_fn(self, t, yvec, fyvec, J, _, tmp1, tmp2, tmp3):
        # For this scalar problem, dF/dy = lambda
        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = self.lamb
        return 0

    def solution(self, y0vec, yvec, t):
        y = N_VGetArrayPointer(yvec)
        y[0] = np.atan(t)
        return 0

    def set_init_cond(self, y0vec):
        y0 = N_VGetArrayPointer(y0vec)
        y0[0] = 0.0
        return 0


class AnalyticMultiscaleODE(ODEProblem):
    """
    * We consider the initial value problem
    *    y' + lambda*y = y^2, y(0) = 1
    * proposed in
    *
    * Estep, D., et al. "An a posteriori–a priori analysis of multiscale operator
    * splitting." SIAM Journal on Numerical Analysis 46.3 (2008): 1116-1146.
    *
    * The parameter lambda is positive, t is in [0, 1], and the exact solution is
    *
    *    y(t) = lambda*y / (y(0) - (y(0) - lambda)*exp(lambda*t))
    """

    T0 = 0.0
    TF = 1.0

    def __init__(self, lamb=2.0):
        self.lamb = lamb

    def f_linear(self, t, yvec, ydotvec, _):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[:] = -self.lamb * y
        return 0

    def f_nonlinear(self, t, yvec, ydotvec, _):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[:] = y * y
        return 0

    def solution(self, y0vec, yvec, tf):
        y0 = N_VGetArrayPointer(y0vec)[0]
        y = N_VGetArrayPointer(yvec)
        y[0] = self.lamb * y0 / (y0 - (y0 - self.lamb) * np.exp(self.lamb * tf))
        return 0

    def set_init_cond(self, y0vec):
        y0 = N_VGetArrayPointer(y0vec)
        y0[0] = 1.0
        return 0


class AnalyticDAE:
    """
    * The following is a simple example problem with analytical
    * solution adapted from example 10.2 of Ascher & Petzold, "Computer
    * Methods for Ordinary Differential Equations and
    * Differential-Algebraic Equations," SIAM, 1998, page 267:
    *    x1'(t) = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t)
    *         0 = (t+2)*x1 - (t+2)*exp(t)
    * for t in the interval [0.0, 1.0], with initial condition:
    *    x1(0) = 1   and   x2(0) = -1/2.
    * The problem has true solution
    *    x1(t) = exp(t)  and  x2(t) = exp(t)/(t-2)
    """

    T0 = 0.0
    TF = 1.0

    def __init__(self, alpha=10.0):
        self.alpha = alpha

    def res(self, t, yyvec, ypvec, resvec, _):
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        res = N_VGetArrayPointer(resvec)
        alpha = self.alpha

        # System residual function:
        #   0 = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t) - x1'(t)
        #   0 = (t+2)*x1 - (t+2)*exp(t)
        res[0] = (
            (1.0 - alpha) / (t - 2.0) * yy[0]
            - yy[0]
            + (alpha - 1.0) * yy[1]
            + 2.0 * np.exp(t)
            - yp[0]
        )
        res[1] = (t + 2.0) * yy[0] - (t + 2.0) * np.exp(t)
        return 0

    def solution(self, yyvec, ypvec, t):
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        yy[0] = np.exp(t)
        yy[1] = np.exp(t) / (t - 2.0)
        yp[0] = np.exp(t)
        yp[1] = np.exp(t) / (t - 2.0) - np.exp(t) / (t - 2.0) / (t - 2.0)
        return 0

    def set_init_cond(self, yyvec, ypvec, t0):
        return self.solution(yyvec, ypvec, t0)

    def psolve(self, t, yyvec, ypvec, rrvec, rvec, zvec, cj, delta):
        """
        Exact solution as preconditioner
        P = df/dy + cj*df/dyp
                =>
                    P = [ - cj - (alpha - 1)/(t - 2) - 1, alpha - 1]
                        [                          t + 2,         0]

        z = P^{-1} r
        */
        """
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        r = N_VGetArrayPointer(rvec)
        z = N_VGetArrayPointer(zvec)
        alpha = self.alpha
        a11 = -cj - (alpha - 1.0) / (t - 2.0) - 1.0
        a12 = alpha - 1.0
        a21 = t + 2.0
        z[0] = r[1] / a21
        z[1] = -(a11 * r[1] - a21 * r[0]) / (a12 * a21)
        return 0


class AnalyticNonlinearSys:
    """
    * This implements the nonlinear system
    *
    * 3x - cos((y-1)z) - 1/2 = 0
    * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
    * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
    *
    * The nonlinear fixed point function is
    *
    * g1(x,y,z) = 1/3 cos((y-1)z) + 1/6
    * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
    * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
    *
    * Corrector form g(x,y,z):
    *
    * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6 - x0
    * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9 - y0
    * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60 - z0
    *
    * This system has the analytic solution x = 1/2, y = 1, z = -pi/6.
    """

    NEQ = 3

    def __init__(self, u0vec=None):
        self.u0vec = u0vec

    # CJB: __enter__ and __exit__ are defined so that this class can be
    # use with python "with" statements. This is a workaround for the following scenario:
    #   u0 = N_VNew_Serial(...)
    #   problem = AnalyticNonlinearSys(u0)
    #   def g_fn(self, u, g, _):
    #      return problem.fixed_point_fn(u, g)
    # Without using a "with" block, this code will cause nanobind to complain about reference
    # leaks because `problem`, which holds a reference to u0, seems to not be garbage collected
    # until the nanobind shutdown callback.
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.u0vec = None

    def fixed_point_fn(self, uvec, gvec, _):
        u = N_VGetArrayPointer(uvec)
        g = N_VGetArrayPointer(gvec)
        x, y, z = u[0], u[1], u[2]

        g[0] = (1.0 / 3.0) * np.cos((y - 1.0) * z) + (1.0 / 6.0)
        g[1] = (1.0 / 9.0) * np.sqrt(x * x + np.sin(z) + 1.06) + 0.9
        g[2] = -(1.0 / 20.0) * np.exp(-x * (y - 1.0)) - (10.0 * np.pi - 3.0) / 60.0

        return 0

    def corrector_fp_fn(self, uvec, gvec):
        self.fixed_point_fn(uvec, gvec, None)
        N_VLinearSum(1.0, gvec, -1.0, self.u0vec, gvec)
        return 0

    def conv_test(self, nls, yvec, delvec, tol, ewtvec):
        delnrm = N_VMaxNorm(delvec)
        if delnrm <= tol:
            return SUN_SUCCESS
        else:
            return SUN_NLS_CONTINUE

    def solution(self, uvec):
        u = N_VGetArrayPointer(uvec)
        u[0] = 0.5
        u[1] = 1.0
        u[2] = -np.pi / 6.0
        return 0
