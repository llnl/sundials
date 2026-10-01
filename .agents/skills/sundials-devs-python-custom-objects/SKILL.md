---
name: sundials-devs-python-custom-objects
description: Add or extend a Python-subclassable implementation of a SUNDIALS object in sundials4py. Use when exposing a C object with an operations table—such as SUNMatrix, SUNLinearSolver, SUNNonlinearSolver, or SUNAdaptController—so Python users can implement its operations and pass instances anywhere the native pointer type is accepted. Do not use for merely wrapping a concrete native constructor or a standalone callback.
---

# Build a Python-custom SUNDIALS object interface

Implement the new family by following the existing custom-object architecture,
not by creating a parallel ownership or conversion mechanism. Start from the
closest existing family in `bindings/sundials4py/include/*_custom.hpp` and reuse
the common machinery in:

- `bindings/sundials4py/include/sundials4py_custom_object.hpp`
- `bindings/sundials4py/include/sundials4py_custom_casters.hpp`
- `bindings/sundials4py/include/sundials4py_types.hpp`

Use `CustomSUNMatrix` as the reference for cloning and mixed ownership,
`CustomSUNLinearSolver` for borrowed return values and independently nullable
callbacks, `CustomSUNNonlinearSolver` for callbacks that depend on an active
opaque memory pointer, and `CustomSUNAdaptController` for related C object types
sharing one Python base implementation.

## 1. Inventory the native object contract

Read the object's public header and implementation before designing its Python
API. Record:

- the opaque handle and concrete struct layout;
- its empty constructor and empty/free/destroy functions;
- every operation-table slot and exact C signature;
- which operations are mandatory, optional, deprecated, or conditional on an
  object type;
- the semantics of absent optional slots;
- ownership of every pointer argument and return value;
- whether SUNDIALS creates additional instances through clone-like operations;
- whether setter operations accept function/data pointer pairs;
- whether callbacks receive their data directly or through a later `mem`
  argument;
- where the owning `SUNContext` is stored and which failure sentinel each
  operation requires.

Inspect callers in `src/` as well as declarations. A slot that appears optional
in a struct may still be assumed non-null by a package, and a returned pointer
may be borrowed, transferred, or dereferenced immediately.

Define the Python method signatures in terms of useful Python objects. Follow
the existing binding convention for output parameters: return a tuple containing
the status followed by output values. Keep opaque C data pointers out of the
Python-facing protocol.

## 2. Add a subclassable C++ base

Create `bindings/sundials4py/include/<object>_custom.hpp` with a C++ base class
registered in the appropriate handwritten binding source.

The base should hold:

- a strong `shared_ptr` owner for its `SUNContext`;
- immutable type information needed by enum-returning operations;
- a cached `shared_ptr` to the lazily created native handle;
- explicit `unmaterialized`, `materializing`, and `materialized` state;
- a materialization counter when useful for focused tests.

Register default methods for the complete Python interface. Required defaults
should raise `NotImplementedError`; optional defaults exist so override
detection can compare descriptors reliably.

Validate required operations when materializing, not in the C++ base
constructor. This lets a Python subclass finish construction before any
override can run. Use `custom_method_overridden<Base>()`, which checks the
concrete type's MRO and supports implementations inherited from user-defined
intermediate bases.

Snapshot optional-operation support during materialization. Later monkey
patching must not add or remove native vtable slots. An already-installed
trampoline may continue using normal dynamic Python lookup, so replacing the
implementation of an installed method can affect later calls.

## 3. Make materialization transactional

The custom type caster calls `_get_sundials_handle(self)` on first conversion.
Implement this as a state machine:

1. Reject a missing base constructor/context.
2. Reject reentry while state is `materializing`.
3. Return the cached handle when already `materialized`.
4. Set state to `materializing`.
5. Build the handle into temporary RAII-managed state.
6. Publish the handle, increment the count, and mark it `materialized` only
   after construction succeeds.
7. On failure, clear the cached handle, restore `unmaterialized`, and rethrow.

Perform required-method validation, optional-method discovery, weak-reference
creation, and other potentially throwing Python work before allocating the C
shell. Hold content in `std::unique_ptr` until it is attached. After attachment,
leave only nonthrowing assignments, or add an explicit cleanup guard for every
later failure path.

Use the established empty-object constructor and matching deleter with
`sundials::experimental::our_make_shared`. Do not invent a second native
lifetime mechanism.

## 4. Model ownership explicitly

Use `SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS()` in the private content struct and
give each family a unique tag. Recover content only through
`custom_content_cast`; native implementations of the same C type can flow
through the same binding and must not be interpreted as custom content.

For a user-created object, store a weak reference from native content back to
the Python implementation. The Python object owns its cached native handle, so
a strong back-reference would form an uncollectable cycle. Document that users
must retain the Python object while SUNDIALS retains its raw pointer.

If SUNDIALS creates a new handle itself, such as through `clone`, give that
handle a strong reference to its Python implementation. The native handle then
owns the implementation and releases it during native destruction.

For borrowed pointers returned by an operation, determine which object owns the
pointee after the trampoline returns. If the Python method can return a
temporary wrapper, retain that wrapper in the custom content and release or
replace it according to the C API's ownership contract.

Destroy content with `custom_content_destroy` or `shutdown_safe_delete`. Native
destruction may occur without the GIL or during interpreter finalization; do not
unconditionally destroy `nb::object` members.

## 5. Implement operation trampolines safely

Each C operation trampoline should:

1. obtain the owning context from the native handle;
2. acquire the GIL;
3. validate and recover tagged custom content;
4. promote the weak implementation reference for the duration of the call;
5. convert arguments with explicit reference policies where appropriate;
6. invoke the Python method;
7. validate and unpack its return value;
8. initialize output parameters to deterministic values on failure;
9. catch every exception with `SUNDIALS4PY_CATCH_AND_REPORT` and return the
   operation's correct failure sentinel.

The shared exception helper reports `SUN_ERR_EXT_FAIL` through the owning
`SUNContext`. For Python exceptions, nanobind's `python_error::what()` preserves
the exception type, message, and traceback. Never permit either the operation
exception or a secondary error-handler exception to cross the C ABI.

Enum/type getters should normally return values cached in native content and
avoid Python calls. Destroy/free operations should be null-safe and should not
raise.

## 6. Add transparent argument conversion

Extend `sundials4py_custom_casters.hpp` for the native pointee type. The caster
must:

- try nanobind's ordinary native-wrapper conversion first;
- recognize the new custom base only if native conversion fails;
- lazily materialize and borrow the cached raw handle for the call;
- preserve existing `None` behavior for nullable pointer arguments;
- leave C-to-Python conversion using the native wrapper representation.

Include the new custom header before the caster definitions in
`sundials4py_types.hpp`. The caster must be visible consistently in every
translation unit that binds the pointer type; piecemeal inclusion can create
different conversion behavior across modules.

Materialization errors occur inside a `noexcept` caster. Follow the established
caster's error-state behavior and verify the actual Python exception users see;
do not assume a C++ exception message will survive nanobind overload handling.

## 7. Wrap callback setters when the object receives native callbacks

If an operation accepts a C callback, also use
`../sundials-devs-python-callbacks/SKILL.md` for the basic function-table and
wrapper pattern.

For custom object interfaces, add these safeguards:

- Wrap each C function/data pair as a normal Python callable.
- Store its `NativeCallbackState` in the custom content's registry.
- Invalidate the previous state whenever a slot is replaced or unset.
- Treat independently nullable callbacks independently; test every supported
  combination rather than collapsing “setup only” into “neither.”
- Invalidate all outstanding adapters before destroying custom content.
- Make stale adapters raise instead of invoking a dangling function or data
  pointer.

Some callbacks receive an opaque package-memory pointer only when the object's
`setup` or `solve` operation runs. For those, use `ActiveMemScope` and record
whether entry came from an integrator or a direct Python binding. The adapter
must reject calls outside the active scope and calls made under the wrong mode.
Do not infer provenance solely by comparing opaque pointer values.

If callback tables are also attached to native implementations through a
`python` field, ensure every native free path destroys that table. When
`SUN*Free()` delegates to `ops->free` and bypasses generic cleanup, interpose on
the original free operation, clear the `python` pointer before delegation,
destroy the table safely under the GIL, and then call the saved operation.

## 8. Integrate generated and handwritten bindings

Register the custom Python base and its methods in the handwritten
`bindings/sundials4py/sundials/*.cpp` file. Keep callback setters, nullable
arguments that litgen cannot express, and ownership-sensitive returns
handwritten.

Update the relevant `generate.yaml` exclusions whenever a generated wrapper is
replaced. Regenerate and commit generated headers using the repository's
documented sundials4py generation workflow; do not leave configuration and
generated output inconsistent.

Audit pointer-returning wrappers for failure. Before wrapping a raw result in a
`shared_ptr`, check for `nullptr` and return Python `None` if that is the API's
failure result. Also inspect the underlying C dispatcher: it must not
dereference a null operation result before returning it.

Avoid public C-ABI changes solely for the Python implementation. A `python`
field is appropriate only when the C object must retain binding callback state;
initialize and destroy it consistently when one is genuinely needed.

## 9. Documentation and examples

Document:

- constructor arguments and valid object-type values;
- required and optional methods;
- exact Python method signatures and tuple return conventions;
- callback callable signatures and nullability;
- recoverable status returns versus unrecoverable Python exceptions;
- object, clone, callback-adapter, context, and borrowed-return lifetimes;
- the requirement that applications retain objects whose pointers SUNDIALS
  stores.

Add an annotated example that exercises the object through a real consuming
package, not only direct `SUN*` operation calls. Update `CHANGELOG.md` and
`doc/shared/RecentChanges.rst` for the user-visible feature. Comments in
published source and tests must be self-contained and must not refer to design
documents that are not shipped.

## 10. Required tests

Keep tests grouped by purpose: per-family operation tests, shared lifecycle and
error tests, and a small integration file. Cover at least:

- native wrappers, custom instances, and permitted `None` arguments through
  both generated and handwritten bindings;
- lazy one-time materialization, failed-materialization rollback, retry, and
  reentry protection;
- required-method validation and inherited optional overrides;
- optional-vtable immutability after materialization;
- argument conversion, tuple returns, output parameters, and failure sentinels;
- Python exceptions reported to `SUNContext` with operation, type, message, and
  traceback, followed by a successful call proving the Python error state was
  consumed;
- a failing Python error handler remaining contained at the C boundary;
- weak user ownership, strong ownership of C-created clones, context lifetime,
  borrowed-return retention/replacement, worker-thread destruction, and
  interpreter finalization;
- callback replacement, unsetting, independent nullability, stale-adapter
  rejection, active-memory misuse, and nested active-memory restoration;
- at least one real integrator or nonlinear-solver path that retains and invokes
  the custom object repeatedly.

Give every test function a short comment block stating the behavior or contract
it protects. Do not test undefined behavior such as deliberately invoking a raw
pointer after its Python owner has been destroyed.

## 11. Validate

Format all changed C/C++ and Python files using the repository-prescribed
versions. Build the sundials4py extension, run focused new tests first, then run
the complete bindings test suite. Finally run `git diff --check` and inspect the
entire diff for generated-file consistency, accidental artifacts, and ownership
or error paths not exercised by tests.

Typical commands are:

```bash
cmake --build <python-build-dir> --target sundials4py -j4
PYTHONPATH=<built-module>:bindings/sundials4py:bindings/sundials4py/test \
  python -m pytest -q bindings/sundials4py/test
git diff --check
```

If build configuration or test discovery is uncertain, also use the
`sundials-build` skill.
