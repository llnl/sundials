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


class HarmonicOscillatorODE(ODEProblem):
    def __init__(self, A=10.0, phi=0.0, omega=1.0):
        self.A = A
        self.phi = phi
        self.omega = omega

    def xdot(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = y[1]
        return 0

    def vdot(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[1] = -self.omega * self.omega * y[0]
        return 0

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        y[0] = self.A * np.cos(self.phi)
        y[1] = -self.A * self.omega * np.sin(self.phi)

    def solution(self, y0vec, yvec, t):
        y = N_VGetArrayPointer(yvec)
        y[0] = self.A * np.cos(self.omega * t + self.phi)
        y[1] = -self.A * self.omega * np.sin(self.omega * t + self.phi)
        return 0
