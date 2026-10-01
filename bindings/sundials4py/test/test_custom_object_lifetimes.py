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
# Lifetime and shutdown behavior of the custom (Python-subclassable) SUNDIALS
# object families.
#
# These tests are about *ownership*, not numerics: who keeps whom alive, when a
# native handle is released, and what happens at interpreter shutdown. They are
# deliberately separate from the per-family operation tests so that a change in
# ownership policy shows up as a failure here rather than as a crash elsewhere.
# -----------------------------------------------------------------

import gc
import subprocess
import sys
import threading
import weakref

import numpy as np
import pytest
from fixtures import *
from sundials4py.core import *


class DiagonalMatrix(CustomSUNMatrix):
    # The smallest complete SUNMatrix implementation, duplicated here rather than
    # imported so that a change to the operations tests cannot silently alter
    # what the lifetime tests are measuring.
    def __init__(self, diagonal, sunctx):
        self.diagonal = np.asarray(diagonal, dtype=sunrealtype).copy()
        super().__init__(sunctx)

    def clone(self):
        return type(self)(self.diagonal.copy(), self.sunctx)

    def zero(self):
        self.diagonal[:] = 0
        return SUN_SUCCESS

    def copy(self, dst):
        dst.diagonal[:] = self.diagonal
        return SUN_SUCCESS

    def scaleadd(self, c, other):
        self.diagonal[:] = c * self.diagonal + other.diagonal
        return SUN_SUCCESS

    def scaleaddi(self, c):
        self.diagonal[:] = c * self.diagonal + 1
        return SUN_SUCCESS

    def matvec(self, x, y):
        N_VGetArrayPointer(y)[:] = self.diagonal * N_VGetArrayPointer(x)
        return SUN_SUCCESS


class TrackedMatrix(DiagonalMatrix):
    # Records a weak reference to every object clone() produces, which is the
    # only handle a test can get on an implementation that native SUNDIALS owns.
    def __init__(self, diagonal, sunctx, clones=None):
        self.clones = clones if clones is not None else []
        super().__init__(diagonal, sunctx)

    def clone(self):
        cloned = TrackedMatrix(self.diagonal.copy(), self.sunctx, self.clones)
        self.clones.append(weakref.ref(cloned))
        return cloned


class TrivialLinearSolver(CustomSUNLinearSolver):
    def __init__(self, sunctx):
        super().__init__(sunctx, SUNLINEARSOLVER_DIRECT)

    def solve(self, A, x, b, tol):
        N_VScale(1.0, b, x)
        return SUN_SUCCESS


class FailOnceDescriptor:
    def __init__(self, method, *, reenter=False):
        self.method = method
        self.fail = True
        self.reenter = reenter
        self.instance = None

    def __get__(self, instance, owner):
        if self.fail:
            self.fail = False
            if self.reenter:
                SUNMatGetID(self.instance)
            raise RuntimeError("deliberate materialization failure")
        return self.method.__get__(instance, owner)


def collect():
    # CPython frees a plain object as soon as its refcount drops, but a cycle
    # through the instance dictionary would need the collector, so ask for it.
    gc.collect()


def test_failed_materialization_is_transactional_and_retryable(sunctx):
    # Purpose:
    # Failed materialization is transactional and retryable.
    class Matrix(DiagonalMatrix):
        matvecsetup = FailOnceDescriptor(lambda self: SUN_SUCCESS)

    A = Matrix([1.0, 2.0], sunctx)
    with pytest.raises(TypeError, match="SUNMatGetID"):
        SUNMatGetID(A)

    assert A._materialization_count() == 0
    assert SUNMatGetID(A) == SUNMATRIX_CUSTOM
    assert A._materialization_count() == 1


def test_reentrant_materialization_resets_state_and_allows_retry(sunctx):
    # Purpose:
    # Reentrant materialization resets state and allows retry.
    descriptor = FailOnceDescriptor(lambda self: SUN_SUCCESS, reenter=True)

    class Matrix(DiagonalMatrix):
        matvecsetup = descriptor

    A = Matrix([1.0, 2.0], sunctx)
    descriptor.instance = A

    with pytest.raises(TypeError, match="SUNMatGetID"):
        SUNMatGetID(A)

    assert A._materialization_count() == 0
    assert SUNMatGetID(A) == SUNMATRIX_CUSTOM
    assert A._materialization_count() == 1


def test_dropping_an_unused_custom_object_releases_its_handle(sunctx):
    # Purpose:
    # Dropping an unused custom object releases its handle.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    assert SUNMatGetID(A) == SUNMATRIX_CUSTOM
    assert A._materialization_count() == 1

    # The native handle holds only a weak reference back to the Python object, so
    # materializing must not turn the object into a permanent resident.
    ref = weakref.ref(A)
    del A
    collect()
    assert ref() is None


def test_retaining_the_python_object_for_the_full_use_period_succeeds(sunctx):
    # Purpose:
    # Retaining the python object for the full use period succeeds.
    # The supported pattern: the Python object outlives every native use of it.
    A = DiagonalMatrix([2.0, 3.0], sunctx)
    x = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(x)[:] = [1.0, 1.0]

    for _ in range(5):
        assert SUNMatMatvec(A, x, y) == SUN_SUCCESS

    # Repeated conversion reuses the one cached handle rather than rebuilding it.
    assert A._materialization_count() == 1
    assert list(N_VGetArrayPointer(y)) == [2.0, 3.0]


def test_unretained_temporary_is_invalid_application_usage(sunctx):
    # Purpose:
    # Unretained temporary is invalid application usage.
    # DOCUMENTS UNSUPPORTED USAGE. Passing a custom object that nothing retains
    # is safe only for the duration of the call it is passed to: the object is
    # collected as soon as the expression ends, and the native handle dies with
    # it. A binding that stored the handle for later use -- as an integrator's
    # linear solver does -- would then be holding a dangling pointer.
    #
    # The test asserts the collection rather than the crash, because triggering
    # the crash is undefined behavior and not something a test should perform.
    tracker = []

    class Recording(DiagonalMatrix):
        def __init__(self, diagonal, sunctx):
            super().__init__(diagonal, sunctx)
            tracker.append(weakref.ref(self))

    assert SUNMatGetID(Recording([1.0, 2.0], sunctx)) == SUNMATRIX_CUSTOM
    collect()

    assert len(tracker) == 1
    assert tracker[0]() is None


def test_c_owned_clones_retain_and_release_their_implementation_once(sunctx):
    # Purpose:
    # C owned clones retain and release their implementation once.
    A = TrackedMatrix([1.0, 2.0], sunctx)
    B = SUNMatClone(A)

    # Python holds no reference to the cloned implementation; the clone's native
    # content holds a strong one, because SUNDIALS owns the clone outright.
    assert len(A.clones) == 1
    collect()
    assert A.clones[0]() is not None
    assert SUNMatScaleAddI(3.0, B) == SUN_SUCCESS

    # Dropping the last handle to the clone runs SUNMatDestroy, which releases
    # that strong reference -- exactly once, or the implementation would leak.
    del B
    collect()
    assert A.clones[0]() is None


def test_clone_of_a_clone_releases_both_implementations(sunctx):
    # Purpose:
    # Clone of a clone releases both implementations.
    # A chain of C-owned clones must unwind completely; a missed release
    # anywhere in the chain would keep every earlier implementation alive.
    A = TrackedMatrix([1.0, 2.0], sunctx)
    B = SUNMatClone(A)
    C = SUNMatClone(B)

    assert len(A.clones) == 2
    del B, C
    collect()
    assert [ref() for ref in A.clones] == [None, None]


def test_custom_object_keeps_its_context_alive():
    # Purpose:
    # Custom object keeps its context alive.
    # The context is reached through a shared_ptr stored both on the Python
    # object and in the native content, so dropping the caller's own reference
    # must not invalidate a matrix that is still in use.
    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    A = DiagonalMatrix([1.0, 2.0], sunctx)
    assert SUNMatGetID(A) == SUNMATRIX_CUSTOM

    del sunctx
    collect()

    assert SUNMatZero(A) == SUN_SUCCESS
    assert A.sunctx is not None


def test_destruction_from_a_worker_thread(sunctx):
    # Purpose:
    # Destruction from a worker thread.
    # The destructor drops Python references, so it must be run with the GIL
    # held. Creating and destroying the object entirely on a worker thread
    # exercises that path on a thread state other than the main one.
    errors = []

    def body():
        try:
            A = DiagonalMatrix([1.0, 2.0], sunctx)
            assert SUNMatGetID(A) == SUNMATRIX_CUSTOM
            B = SUNMatClone(A)
            assert SUNMatScaleAddI(1.0, B) == SUN_SUCCESS
        except BaseException as exc:  # pragma: no cover - reported below
            errors.append(exc)

    worker = threading.Thread(target=body)
    worker.start()
    worker.join()
    collect()

    assert errors == []


def test_custom_linear_solver_handle_is_released_with_its_object(sunctx):
    # Purpose:
    # Custom linear solver handle is released with its object.
    LS = TrivialLinearSolver(sunctx)
    x = N_VNew_Serial(2, sunctx)
    b = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(b)[:] = [5.0, 6.0]

    assert SUNLinSolSolve(LS, None, x, b, 0.0) == SUN_SUCCESS
    assert list(N_VGetArrayPointer(x)) == [5.0, 6.0]

    ref = weakref.ref(LS)
    del LS
    collect()
    assert ref() is None


# Run in a subprocess so that the interpreter really does finalize with live
# custom objects still reachable from module globals. nanobind's contract is
# that leaked objects at finalization are reported but do not crash, so the exit
# status is what this test checks.
_FINALIZATION_SCRIPT = """
import numpy as np
from sundials4py.core import *


class DiagonalMatrix(CustomSUNMatrix):
    def __init__(self, diagonal, sunctx):
        self.diagonal = np.asarray(diagonal, dtype=sunrealtype).copy()
        super().__init__(sunctx)

    def clone(self):
        return type(self)(self.diagonal.copy(), self.sunctx)

    def zero(self):
        self.diagonal[:] = 0
        return SUN_SUCCESS

    def copy(self, dst):
        dst.diagonal[:] = self.diagonal
        return SUN_SUCCESS

    def scaleadd(self, c, other):
        self.diagonal[:] = c * self.diagonal + other.diagonal
        return SUN_SUCCESS

    def scaleaddi(self, c):
        self.diagonal[:] = c * self.diagonal + 1
        return SUN_SUCCESS

    def matvec(self, x, y):
        N_VGetArrayPointer(y)[:] = self.diagonal * N_VGetArrayPointer(x)
        return SUN_SUCCESS


status, sunctx = SUNContext_Create(SUN_COMM_NULL)

# Deliberately kept alive to the end of the process: a materialized custom
# object, and a clone that native SUNDIALS owns.
KEPT_ALIVE = DiagonalMatrix([1.0, 2.0], sunctx)
assert SUNMatGetID(KEPT_ALIVE) == SUNMATRIX_CUSTOM
KEPT_CLONE = SUNMatClone(KEPT_ALIVE)
assert SUNMatScaleAddI(1.0, KEPT_CLONE) == SUN_SUCCESS
print("ok")
"""


def test_interpreter_finalization_with_live_custom_objects():
    # Purpose:
    # Interpreter finalization with live custom objects.
    result = subprocess.run(
        [sys.executable, "-c", _FINALIZATION_SCRIPT], capture_output=True, text=True
    )

    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout
