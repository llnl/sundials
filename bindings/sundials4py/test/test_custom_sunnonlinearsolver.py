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

import pytest
from fixtures import *
from numpy.testing import assert_allclose
from sundials4py.core import *


class CopyNonlinearSolver(CustomSUNNonlinearSolver):
    # The solve operation copies y0 into y, giving a deterministic check that
    # C-level nonlinear solver calls reach the Python subclass.
    def __init__(self, sunctx, solver_type=SUNNONLINEARSOLVER_ROOTFIND):
        self.calls = {"initialize": 0, "setup": 0, "solve": 0, "set_max_iters": 0}
        self.max_iters = None
        super().__init__(sunctx, solver_type)

    def initialize(self):
        self.calls["initialize"] += 1
        return SUN_SUCCESS

    def setup(self, y):
        self.calls["setup"] += 1
        self.last_setup_y = y
        return SUN_SUCCESS

    def solve(self, y0, y, w, tol, call_lsetup):
        self.calls["solve"] += 1
        self.last_tol = tol
        self.last_call_lsetup = call_lsetup
        N_VGetArrayPointer(y)[:] = N_VGetArrayPointer(y0)
        return SUN_SUCCESS

    def set_max_iters(self, maxiters):
        self.calls["set_max_iters"] += 1
        self.max_iters = maxiters
        return SUN_SUCCESS

    def get_num_iters(self):
        return SUN_SUCCESS, 4

    def get_cur_iter(self):
        return SUN_SUCCESS, 2

    def get_num_conv_fails(self):
        return SUN_SUCCESS, 1


class CallbackNonlinearSolver(CustomSUNNonlinearSolver):
    # Records the Python adapters SUNDIALS hands to each callback setter. Rather
    # than running a real nonlinear iteration, solve() calls the adapters so a
    # test can observe a callback installed from Python travelling out through
    # the C layer and back into Python.
    def __init__(self, sunctx):
        self.sys_fn = None
        self.root_sys_fn = None
        self.fixed_point_sys_fn = None
        self.lsetup_fn = None
        self.lsolve_fn = None
        self.conv_test_fn = None
        self.norm_fn = None
        self.get_update_norm_fn = None
        self.get_conv_rate_fn = None
        self.results = {}
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    def set_sys_fn(self, sys_fn):
        self.sys_fn = sys_fn
        return SUN_SUCCESS

    def set_sys_fns(self, root_fn, fixed_point_fn):
        self.root_sys_fn = root_fn
        self.fixed_point_sys_fn = fixed_point_fn
        return SUN_SUCCESS

    def set_lsetup_fn(self, lsetup_fn):
        self.lsetup_fn = lsetup_fn
        return SUN_SUCCESS

    def set_lsolve_fn(self, lsolve_fn):
        self.lsolve_fn = lsolve_fn
        return SUN_SUCCESS

    def set_conv_test_fn(self, conv_test_fn):
        self.conv_test_fn = conv_test_fn
        return SUN_SUCCESS

    def set_norm_fn(self, norm_fn):
        self.norm_fn = norm_fn
        return SUN_SUCCESS

    def set_get_update_norm_fn(self, get_update_norm_fn):
        self.get_update_norm_fn = get_update_norm_fn
        return SUN_SUCCESS

    def set_get_conv_rate_fn(self, get_conv_rate_fn):
        self.get_conv_rate_fn = get_conv_rate_fn
        return SUN_SUCCESS

    def setup(self, y):
        # Called from inside solve() below to check that the active-memory scope
        # nests: the inner scope must not destroy the outer one.
        self.results["setup_sys"] = self.sys_fn(y, y)
        return SUN_SUCCESS

    def solve(self, y0, y, w, tol, call_lsetup):
        # Each adapter below may only be invoked while this solve() is on the
        # stack; that is exactly what the active-memory scope enforces.
        if self.sys_fn is not None:
            self.results["sys"] = self.sys_fn(y0, y)
        if self.lsetup_fn is not None:
            self.results["lsetup"] = self.lsetup_fn(True)
        if self.lsolve_fn is not None:
            self.results["lsolve"] = self.lsolve_fn(y)
        return SUN_SUCCESS


class IncompleteNonlinearSolver(CustomSUNNonlinearSolver):
    # Required-method validation should reject subclasses that inherit solve().
    def __init__(self, sunctx):
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)


def test_custom_sunnonlinearsolver_type_and_lazy_materialization(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver type and lazy materialization.
    NLS = CopyNonlinearSolver(sunctx, SUNNONLINEARSOLVER_HYBRID)

    assert NLS._materialization_count() == 0
    assert SUNNonlinSolGetType(NLS) == SUNNONLINEARSOLVER_HYBRID
    assert NLS._materialization_count() == 1


def test_custom_sunnonlinearsolver_required_solve_is_validated(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver required solve is validated.
    NLS = IncompleteNonlinearSolver(sunctx)

    with pytest.raises(TypeError, match="SUNNonlinSolGetType"):
        SUNNonlinSolGetType(NLS)

    assert NLS._materialization_count() == 0


def test_custom_sunnonlinearsolver_setup_and_solve(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver setup and solve.
    NLS = CopyNonlinearSolver(sunctx)
    y0 = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    w = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(y0)[:] = [3.0, 7.0]

    assert SUNNonlinSolInitialize(NLS) == SUN_SUCCESS
    assert SUNNonlinSolSetup(NLS, y0) == SUN_SUCCESS
    assert SUNNonlinSolSolve(NLS, y0, y, w, 1.0e-8, True) == SUN_SUCCESS

    assert NLS.calls["initialize"] == 1
    assert NLS.calls["setup"] == 1
    assert NLS.calls["solve"] == 1
    assert NLS.last_tol == float(sunrealtype(1.0e-8))
    assert NLS.last_call_lsetup
    assert_allclose(N_VGetArrayPointer(y), [3.0, 7.0])


def test_custom_sunnonlinearsolver_optional_methods(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver optional methods.
    NLS = CopyNonlinearSolver(sunctx)

    assert SUNNonlinSolSetMaxIters(NLS, 12) == SUN_SUCCESS
    assert SUNNonlinSolGetNumIters(NLS) == (SUN_SUCCESS, 4)
    assert SUNNonlinSolGetCurIter(NLS) == (SUN_SUCCESS, 2)
    assert SUNNonlinSolGetNumConvFails(NLS) == (SUN_SUCCESS, 1)
    assert NLS.calls["set_max_iters"] == 1
    assert NLS.max_iters == 12


def test_custom_sunnonlinearsolver_scoped_adapters_round_trip(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver scoped adapters round trip.
    # sys, lsetup, and lsolve receive their memory pointer from the enclosing
    # setup()/solve() call, so their adapters are usable only from within one.
    NLS = CallbackNonlinearSolver(sunctx)
    y0 = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    w = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(y0)[:] = [1.0, 2.0]
    calls = []

    def sys_fn(u, g, _):
        calls.append("sys")
        N_VGetArrayPointer(g)[:] = 2.0 * N_VGetArrayPointer(u)
        return SUN_SUCCESS

    def lsetup_fn(jbad, _):
        calls.append(("lsetup", bool(jbad)))
        return SUN_SUCCESS, True

    def lsolve_fn(b, _):
        calls.append("lsolve")
        return SUN_SUCCESS

    assert SUNNonlinSolSetSysFn(NLS, sys_fn) == SUN_SUCCESS
    assert SUNNonlinSolSetLSetupFn(NLS, lsetup_fn) == SUN_SUCCESS
    assert SUNNonlinSolSetLSolveFn(NLS, lsolve_fn) == SUN_SUCCESS

    # The subclass is handed ordinary Python callables, never SUNDIALS' opaque
    # function and data pointers.
    assert callable(NLS.sys_fn)
    assert callable(NLS.lsetup_fn)
    assert callable(NLS.lsolve_fn)

    assert SUNNonlinSolSolve(NLS, y0, y, w, 1.0e-9, True) == SUN_SUCCESS

    assert calls == ["sys", ("lsetup", True), "lsolve"]
    assert NLS.results["sys"] == SUN_SUCCESS
    assert NLS.results["lsetup"] == (SUN_SUCCESS, True)
    assert NLS.results["lsolve"] == SUN_SUCCESS
    assert_allclose(N_VGetArrayPointer(y), [2.0, 4.0])


def test_custom_sunnonlinearsolver_scoped_adapter_rejects_use_outside_solve(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver scoped adapter rejects use outside solve.
    NLS = CallbackNonlinearSolver(sunctx)
    y = N_VNew_Serial(2, sunctx)

    def sys_fn(u, g, _):
        return SUN_SUCCESS

    assert SUNNonlinSolSetSysFn(NLS, sys_fn) == SUN_SUCCESS

    # No setup() or solve() is in progress, so there is no memory pointer to
    # pass through and the adapter must refuse rather than guess.
    with pytest.raises(RuntimeError, match="setup\\(\\) or solve\\(\\)"):
        NLS.sys_fn(y, y)


def test_custom_sunnonlinearsolver_scoped_adapters_nest(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver scoped adapters nest.
    # A package may call setup() from inside solve(); the inner scope must
    # restore the outer one rather than clear it.
    NLS = CallbackNonlinearSolver(sunctx)
    y0 = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    w = N_VNew_Serial(2, sunctx)

    def sys_fn(u, g, _):
        return SUN_SUCCESS

    assert SUNNonlinSolSetSysFn(NLS, sys_fn) == SUN_SUCCESS

    # Re-entering through the C API opens a second scope inside the first. The
    # adapter must still work after that inner scope has closed, which it can
    # only do if the inner scope restored the outer one instead of clearing it.
    def solve(y0, y, w, tol, call_lsetup):
        assert SUNNonlinSolSetup(NLS, y0) == SUN_SUCCESS
        NLS.results["sys"] = NLS.sys_fn(y0, y)
        return SUN_SUCCESS

    NLS.solve = solve
    assert SUNNonlinSolSolve(NLS, y0, y, w, 1.0e-9, False) == SUN_SUCCESS
    assert NLS.results["setup_sys"] == SUN_SUCCESS
    assert NLS.results["sys"] == SUN_SUCCESS


def test_custom_sunnonlinearsolver_replacing_a_callback_revokes_the_old_adapter(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver replacing a callback revokes the old adapter.
    NLS = CallbackNonlinearSolver(sunctx)
    y = N_VNew_Serial(2, sunctx)

    def first(u, g, _):
        return SUN_SUCCESS

    def second(u, g, _):
        return SUN_SUCCESS

    assert SUNNonlinSolSetSysFn(NLS, first) == SUN_SUCCESS
    stale = NLS.sys_fn

    assert SUNNonlinSolSetSysFn(NLS, second) == SUN_SUCCESS
    assert NLS.sys_fn is not stale

    # Calling the displaced adapter would dispatch through a function pointer
    # SUNDIALS has already replaced, so it is revoked instead.
    with pytest.raises(RuntimeError, match="no longer valid"):
        stale(y, y)


def test_custom_sunnonlinearsolver_data_carrying_adapters(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver data carrying adapters.
    # These four callbacks arrive with their own data pointer, so their adapters
    # do not consult the active-memory scope and work outside setup()/solve().
    NLS = CallbackNonlinearSolver(sunctx)
    y = N_VNew_Serial(2, sunctx)
    delta = N_VNew_Serial(2, sunctx)
    ewt = N_VNew_Serial(2, sunctx)
    w = N_VNew_Serial(2, sunctx)
    calls = []

    def conv_test(nls, u, d, tol, weight, _):
        calls.append("conv_test")
        return SUN_SUCCESS

    def norm(d, weight, _):
        calls.append("norm")
        return SUN_SUCCESS, 2.5

    def get_update_norm(_):
        calls.append("get_update_norm")
        return SUN_SUCCESS, 0.25

    def get_conv_rate(_):
        calls.append("get_conv_rate")
        return SUN_SUCCESS, 0.5

    assert SUNNonlinSolSetConvTestFn(NLS, conv_test) == SUN_SUCCESS
    assert SUNNonlinSolSetNormFn(NLS, norm) == SUN_SUCCESS
    assert SUNNonlinSolSetGetUpdateNormFn(NLS, get_update_norm) == SUN_SUCCESS
    assert SUNNonlinSolSetGetConvRateFn(NLS, get_conv_rate) == SUN_SUCCESS

    assert NLS.conv_test_fn(y, delta, 1.0e-3, ewt) == SUN_SUCCESS
    assert NLS.norm_fn(delta, w) == (SUN_SUCCESS, 2.5)
    assert NLS.get_update_norm_fn() == (SUN_SUCCESS, 0.25)
    assert NLS.get_conv_rate_fn() == (SUN_SUCCESS, 0.5)
    assert calls == ["conv_test", "norm", "get_update_norm", "get_conv_rate"]


def test_custom_sunnonlinearsolver_set_sys_fns_receives_both_adapters(sunctx):
    # Purpose:
    # Custom sunnonlinearsolver set sys fns receives both adapters.
    NLS = CallbackNonlinearSolver(sunctx)
    y0 = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    w = N_VNew_Serial(2, sunctx)
    calls = []

    def root_fn(u, g, _):
        calls.append("root")
        return SUN_SUCCESS

    def fixed_point_fn(u, g, _):
        calls.append("fixed_point")
        return SUN_SUCCESS

    assert SUNNonlinSolSetSysFns(NLS, root_fn, fixed_point_fn) == SUN_SUCCESS
    assert callable(NLS.root_sys_fn)
    assert callable(NLS.fixed_point_sys_fn)

    # Both slots are independent adapters over the same memory scope.
    def solve(y0, y, w, tol, call_lsetup):
        NLS.root_sys_fn(y0, y)
        NLS.fixed_point_sys_fn(y0, y)
        return SUN_SUCCESS

    NLS.solve = solve
    assert SUNNonlinSolSolve(NLS, y0, y, w, 1.0e-9, False) == SUN_SUCCESS
    assert calls == ["root", "fixed_point"]


def test_native_sunnonlinearsolver_conversion_still_works(sunctx):
    # Purpose:
    # Native sunnonlinearsolver conversion still works.
    y = N_VNew_Serial(2, sunctx)
    NLS = SUNNonlinSol_Newton(y, sunctx)

    assert SUNNonlinSolGetType(NLS) == SUNNONLINEARSOLVER_ROOTFIND
