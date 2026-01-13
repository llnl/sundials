#!/usr/bin/env python3

import pytest
import numpy as np
from sundials4py.core import *
from sundials4py.test import SUNContext_TestErrHandler


def test_with_null_comm():
    # Create a new context with a null comm
    err, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert err == SUN_SUCCESS
    assert sunctx is not None

    # Try calling a SUNContext_ function
    last_err = SUNContext_GetLastError(sunctx)

    assert last_err == SUN_SUCCESS


def test_push_pop_err_handlers():
    # Create a new context with a null comm
    err, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert err == SUN_SUCCESS
    assert sunctx is not None

    called = {"err_fn1": False, "err_fn2": False}

    def err_fn1(line, func_name, file_name, msg, err_code, _, sunctx):
        # err_fn2 should already be called since it was pushed second
        assert called["err_fn2"]
        called["err_fn1"] = True

    def err_fn2(line, func_name, file_name, msg, err_code, _, sunctx):
        called["err_fn2"] = True

    status = SUNContext_PushErrHandler(sunctx, err_fn1)
    assert status == SUN_SUCCESS

    status = SUNContext_PushErrHandler(sunctx, err_fn2)
    assert status == SUN_SUCCESS

    SUNContext_TestErrHandler(sunctx)
    assert called["err_fn1"]

    called = {"err_fn1": False, "err_fn2": False}

    status = SUNContext_PopErrHandler(sunctx)
    assert status == SUN_SUCCESS

    status = SUNContext_PopErrHandler(sunctx)
    assert status == SUN_SUCCESS

    SUNContext_TestErrHandler(sunctx)
    assert not called["err_fn1"]
    assert not called["err_fn2"]

    # Popping again should do nothing
    status = SUNContext_PopErrHandler(sunctx)
    assert status == SUN_SUCCESS
