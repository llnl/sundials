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
from sundials4py.core import *


class HController(CustomSUNHController):
    # H controllers return one output value, so estimate_step returns
    # (status, hnew) to mirror the Python binding convention.
    def __init__(self, sunctx):
        self.calls = {"estimate_step": 0, "reset": 0, "set_defaults": 0, "update_h": 0}
        self.bias = None
        super().__init__(sunctx)

    def estimate_step(self, h, p, dsm):
        self.calls["estimate_step"] += 1
        return SUN_SUCCESS, h / (p + dsm)

    def reset(self):
        self.calls["reset"] += 1
        return SUN_SUCCESS

    def set_defaults(self):
        self.calls["set_defaults"] += 1
        return SUN_SUCCESS

    def set_error_bias(self, bias):
        self.bias = bias
        return SUN_SUCCESS

    def update_h(self, h, dsm):
        self.calls["update_h"] += 1
        self.last_update_h = (h, dsm)
        return SUN_SUCCESS


class MRIController(CustomSUNMRIController):
    # MRI controllers return two output values, giving separate coverage for the
    # estimate_step_tol trampoline and optional update hook.
    def __init__(self, sunctx):
        self.calls = {"estimate_step_tol": 0, "update_mri_h_tol": 0}
        super().__init__(sunctx)

    def estimate_step_tol(self, H, tolfac, P, DSM, dsm):
        self.calls["estimate_step_tol"] += 1
        return SUN_SUCCESS, H / P, tolfac * (DSM + dsm)

    def update_mri_h_tol(self, H, tolfac, DSM, dsm):
        self.calls["update_mri_h_tol"] += 1
        self.last_update = (H, tolfac, DSM, dsm)
        return SUN_SUCCESS


class IncompleteHController(CustomSUNHController):
    # Missing estimate_step() should prevent native handle materialization.
    pass


def test_custom_hcontroller_type_and_estimate(sunctx):
    # Purpose:
    # Custom hcontroller type and estimate.
    C = HController(sunctx)

    assert C._materialization_count() == 0
    assert SUNAdaptController_GetType(C) == SUN_ADAPTCONTROLLER_H
    status, hnew = SUNAdaptController_EstimateStep(C, 4.0, 3, 1.0)

    assert status == SUN_SUCCESS
    assert hnew == 1.0
    assert C.calls["estimate_step"] == 1
    assert C._materialization_count() == 1


def test_custom_hcontroller_required_method_is_validated(sunctx):
    # Purpose:
    # Custom hcontroller required method is validated.
    C = IncompleteHController(sunctx)

    with pytest.raises(TypeError, match="SUNAdaptController_GetType"):
        SUNAdaptController_GetType(C)

    assert C._materialization_count() == 0


def test_custom_hcontroller_optional_methods(sunctx):
    # Purpose:
    # Custom hcontroller optional methods.
    C = HController(sunctx)

    assert SUNAdaptController_Reset(C) == SUN_SUCCESS
    assert SUNAdaptController_SetDefaults(C) == SUN_SUCCESS
    assert SUNAdaptController_SetErrorBias(C, 1.25) == SUN_SUCCESS
    assert SUNAdaptController_UpdateH(C, 0.5, 0.75) == SUN_SUCCESS

    assert C.calls["reset"] == 1
    assert C.calls["set_defaults"] == 1
    assert C.bias == 1.25
    assert C.last_update_h == (0.5, 0.75)


def test_custom_mricontroller_type_estimate_and_optional_update(sunctx):
    # Purpose:
    # Custom mricontroller type estimate and optional update.
    C = MRIController(sunctx)

    assert SUNAdaptController_GetType(C) == SUN_ADAPTCONTROLLER_MRI_H_TOL
    status, Hnew, tolfacnew = SUNAdaptController_EstimateStepTol(C, 8.0, 2.0, 4, 0.5, 0.25)
    assert status == SUN_SUCCESS
    assert Hnew == 2.0
    assert tolfacnew == 1.5
    assert SUNAdaptController_UpdateMRIHTol(C, 8.0, 2.0, 0.5, 0.25) == SUN_SUCCESS
    assert C.calls["estimate_step_tol"] == 1
    assert C.calls["update_mri_h_tol"] == 1


def test_native_sunadaptcontroller_conversion_still_works(sunctx):
    # Purpose:
    # Native sunadaptcontroller conversion still works.
    C = SUNAdaptController_Soderlind(sunctx)

    assert SUNAdaptController_GetType(C) == SUN_ADAPTCONTROLLER_H
