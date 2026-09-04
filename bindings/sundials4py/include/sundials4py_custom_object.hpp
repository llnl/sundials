/*------------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *------------------------------------------------------------------------------
 * Machinery shared by every Python-implemented ("custom") SUNDIALS object.
 *
 * Custom SUNDIALS objects follow a repeated pattern: a Python subclass stays an
 * ordinary nanobind object until some binding needs a raw SUNDIALS handle, at
 * which point a native "shell" object is built whose operation table dispatches
 * back into Python. The pieces of that pattern which are genuinely common live
 * here so that the policy is stated exactly once:
 *
 *   ShutdownSafeGIL           interpreter-finalization-safe GIL acquisition
 *   NativeCallbackState       revocable capture of C function/data pairs
 *   NativeCallbackRegistry    per-vtable-slot ownership of the above
 *   ActiveMemMode/Scope       explicit tracking of the integrator "mem" pointer
 *   custom_method_overridden  optional-operation detection via the MRO
 *   custom_content_cast       tag-checked recovery of custom content
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_CUSTOM_OBJECT_HPP
#define _SUNDIALS4PY_CUSTOM_OBJECT_HPP

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_types.h>

#include "sundials4py_core_types.hpp"

namespace sundials4py {

/*------------------------------------------------------------------------------
 * Interpreter shutdown
 *----------------------------------------------------------------------------*/

/*
 * GIL acquisition that is safe to use from a destructor.
 *
 * Custom object content is destroyed from SUNDIALS' C free/destroy operations,
 * which may run on a thread that does not hold the GIL, and may also run while
 * the interpreter is being torn down (for example when a SUNContext held by a
 * module-level Python object is collected during finalization). nanobind does
 * not permit acquiring the GIL or touching Python objects once its internals
 * have been dismantled, so this helper reports whether Python is still usable
 * instead of unconditionally acquiring.
 *
 * Callers that find Python unavailable must skip destruction of any nb::object
 * state. Leaking a handful of references as the process exits is the correct
 * trade: the alternative is undefined behavior during finalization.
 */
class ShutdownSafeGIL
{
public:
  ShutdownSafeGIL() : available_(Py_IsInitialized() != 0 && nb::is_alive())
  {
    if (available_) { gil_.emplace(); }
  }

  ShutdownSafeGIL(const ShutdownSafeGIL&)            = delete;
  ShutdownSafeGIL& operator=(const ShutdownSafeGIL&) = delete;

  /* True when it is safe to run Python code and destroy nb::object state. */
  bool python_available() const { return available_; }

private:
  bool available_;
  std::optional<nb::gil_scoped_acquire> gil_;
};

/*
 * Delete state containing Python references only while Python is usable.
 *
 * SUNDIALS free operations can run on threads that did not enter through
 * Python, and can also run during interpreter finalization. In the latter case
 * the state is intentionally leaked because decrementing Python references is
 * no longer safe and the process is already exiting.
 */
template<typename T>
void shutdown_safe_delete(T* state)
{
  if (!state) { return; }

  ShutdownSafeGIL gil;
  if (gil.python_available()) { delete state; }
}

/*------------------------------------------------------------------------------
 * Exception translation at C callback boundaries
 *----------------------------------------------------------------------------*/

/*
 * Report a C++ exception through the owning SUNContext without allowing a
 * second exception (for example, from a Python error handler) to escape the C
 * callback boundary. nanobind::python_error::what() includes the Python
 * exception type, message, and traceback, so the std::exception overload also
 * preserves the useful Python diagnostic.
 */
inline void report_custom_exception(SUNContext sunctx, const char* operation,
                                    const std::exception& error) noexcept
{
  try
  {
    ShutdownSafeGIL gil;
    if (!gil.python_available()) { return; }

    std::string message = std::string("exception in ") + operation + ":\n" +
                          error.what();
    if (sunctx)
    {
      SUNHandleErrWithMsg(__LINE__, operation, __FILE__, message.c_str(),
                          SUN_ERR_EXT_FAIL, sunctx);
    }
    else
    {
      SUNGlobalFallbackErrHandler(__LINE__, operation, __FILE__, "%s",
                                  SUN_ERR_EXT_FAIL, message.c_str());
    }
  }
  catch (...)
  {
    /* Error reporting is itself on a C ABI boundary. In particular, a Python
       SUNContext error handler is allowed to fail, but that failure must not
       replace the original callback failure or unwind into SUNDIALS. */
    if (Py_IsInitialized() != 0 && nb::is_alive())
    {
      nb::gil_scoped_acquire gil;
      PyErr_Clear();
    }
  }
}

inline void report_custom_unknown_exception(SUNContext sunctx,
                                            const char* operation) noexcept
{
  const std::runtime_error error("unknown non-standard C++ exception");
  report_custom_exception(sunctx, operation, error);
}

/* Every trampoline supplies its own failure sentinel, since SUNDIALS callback
   APIs return status codes, pointers, integers, and real scalars. */
#define SUNDIALS4PY_CATCH_AND_REPORT(SUNCTX, OPERATION, FAILURE)       \
  catch (const std::exception& error)                                  \
  {                                                                    \
    ::sundials4py::report_custom_exception(SUNCTX, OPERATION, error);  \
    return FAILURE;                                                    \
  }                                                                    \
  catch (...)                                                          \
  {                                                                    \
    ::sundials4py::report_custom_unknown_exception(SUNCTX, OPERATION); \
    return FAILURE;                                                    \
  }

/*------------------------------------------------------------------------------
 * Explicit active-memory modes
 *----------------------------------------------------------------------------*/

/*
 * Several SUNNonlinearSolver callbacks take an opaque "mem" pointer that is only
 * supplied to the solver when SUNDIALS enters setup() or solve(). There are two
 * distinct providers of that pointer, and confusing them would hand a package's
 * private memory to a callback expecting the binding's function table (or the
 * reverse), so the provider is tracked explicitly rather than inferred.
 */
enum class ActiveMemMode
{
  /* No setup()/solve() call is in progress; mem-dependent callbacks are
     invalid. */
  none,

  /* A SUNDIALS package (CVODES, IDAS, ARKODE, KINSOL) drove the custom vtable
     and supplied its own integrator memory. */
  integrator,

  /* A hand-written sundials4py wrapper was the caller, so the relevant pointer
     is the object's `python` function table, matching the pre-existing native
     callback wrappers. */
  direct_binding
};

inline const char* active_mem_mode_name(ActiveMemMode mode)
{
  switch (mode)
  {
  case ActiveMemMode::integrator: return "integrator";
  case ActiveMemMode::direct_binding: return "direct-binding";
  default: return "none";
  }
}

/*
 * Marks the dynamic extent during which a hand-written sundials4py wrapper --
 * rather than a SUNDIALS package -- is driving a custom object.
 *
 * Two things need this. A setter (SUNNonlinSolSetSysFn and friends) uses it so
 * the custom trampoline can record that the callback being installed expects the
 * binding's function table. An entry point (SUNNonlinSolSetup/Solve) uses it so
 * the custom trampoline establishes direct-binding rather than integrator mode.
 *
 * Determining this by comparing the incoming `mem` against `NLS->python` would
 * also work, but it infers provenance from pointer identity; stating it at the
 * one place that knows the answer is both cheaper and harder to get wrong.
 */
class DirectBindingScope
{
public:
  DirectBindingScope() : saved_(flag()) { flag() = true; }

  ~DirectBindingScope() { flag() = saved_; }

  DirectBindingScope(const DirectBindingScope&)            = delete;
  DirectBindingScope& operator=(const DirectBindingScope&) = delete;

  static bool active() { return flag(); }

private:
  /* The flag lives inside an inline function rather than being a
     `static thread_local` data member: Apple clang emits the thread-local
     wrapper routine for an inline static data member with strong external
     linkage, so every translation unit that includes this header contributes a
     definition and the link fails with a duplicate symbol. A function-local
     thread_local in an inline function has the vague linkage we actually want. */
  static bool& flag()
  {
    static thread_local bool active = false;
    return active;
  }

  bool saved_;
};

/* The mode a callback installed right now should later be invoked under. */
inline ActiveMemMode current_install_mode()
{
  return DirectBindingScope::active() ? ActiveMemMode::direct_binding
                                      : ActiveMemMode::integrator;
}

/*------------------------------------------------------------------------------
 * Revocable native callback state
 *----------------------------------------------------------------------------*/

/*
 * Type-erased base so that a single registry can invalidate every outstanding
 * adapter state without knowing the concrete C callback signature.
 */
struct NativeCallbackStateBase
{
  virtual ~NativeCallbackStateBase() = default;

  /* Cleared when the owning vtable slot is replaced, unset, or destroyed. */
  bool valid{false};
};

/*
 * A C function pointer plus its opaque data pointer, held indirectly so that
 * the binding can revoke it.
 *
 * Adapters exposed to Python capture a shared_ptr to one of these rather than
 * copying the raw pair. If SUNDIALS later replaces the callback, or the custom
 * object is destroyed, the state is marked invalid and any Python callable the
 * user still holds raises a clear exception instead of calling through a stale
 * function pointer with a dangling data pointer.
 */
template<typename Fn>
struct NativeCallbackState : NativeCallbackStateBase
{
  Fn fn{nullptr};
  void* data{nullptr};

  // Which provider's memory this callback expects, for the subset of callbacks
  // whose data pointer only arrives at call time. Callbacks that
  // carry an explicit data pointer ignore this.
  ActiveMemMode required_mode{ActiveMemMode::integrator};
};

/*
 * Owns one live NativeCallbackState per named vtable slot.
 *
 * Installing into a slot invalidates whatever was there before, which is
 * exactly the semantics SUNDIALS setters have: setting a new system function
 * means the previous one must never be called again.
 */
class NativeCallbackRegistry
{
public:
  template<typename Fn>
  std::shared_ptr<NativeCallbackState<Fn>> install(const char* slot, Fn fn,
                                                   void* data)
  {
    invalidate(slot);
    if (!fn) { return nullptr; }

    auto state           = std::make_shared<NativeCallbackState<Fn>>();
    state->fn            = fn;
    state->data          = data;
    state->required_mode = current_install_mode();
    state->valid         = true;
    slots_[slot]         = state;
    return state;
  }

  /* Revoke the state held for one slot, if any. */
  void invalidate(const char* slot)
  {
    auto it = slots_.find(slot);
    if (it == slots_.end()) { return; }
    if (it->second) { it->second->valid = false; }
    slots_.erase(it);
  }

  /* Revoke every slot; called when the custom object is destroyed. */
  void invalidate_all()
  {
    for (auto& entry : slots_)
    {
      if (entry.second) { entry.second->valid = false; }
    }
    slots_.clear();
  }

private:
  std::map<std::string, std::shared_ptr<NativeCallbackStateBase>> slots_;
};

/*
 * Guard used at the top of every adapter body. Raising here turns into a normal
 * Python exception because adapters are always invoked from Python.
 */
inline void require_valid_callback(const NativeCallbackStateBase* state,
                                   const char* what)
{
  if (!state || !state->valid)
  {
    throw std::runtime_error(
      std::string("[sundials4py] the native ") + what +
      " callback adapter is no longer valid; SUNDIALS has replaced or released "
      "it, so calling it now would use a stale function and data pointer");
  }
}

/*------------------------------------------------------------------------------
 * The active-memory scope itself
 *----------------------------------------------------------------------------*/

/*
 * RAII scope that records the active memory pointer and mode for the duration
 * of one custom setup()/solve() call.
 *
 * Nesting is supported (a package may call setup() from inside solve()) because
 * the previous value is saved and restored, and restoration happens on every
 * exit path including exceptions since it lives in the destructor.
 *
 * Templated on the content type so that each object family can keep its own
 * private content struct; the only requirement is that it expose the members
 * `active_mem_mode` and `active_mem`.
 */
template<typename ContentT>
class ActiveMemScope
{
public:
  ActiveMemScope(ContentT* content, ActiveMemMode mode, void* mem) noexcept
    : content_(content)
  {
    if (!content_) { return; }
    saved_mode_               = content_->active_mem_mode;
    saved_mem_                = content_->active_mem;
    content_->active_mem_mode = mode;
    content_->active_mem      = mem;
  }

  ~ActiveMemScope() noexcept
  {
    if (!content_) { return; }
    content_->active_mem_mode = saved_mode_;
    content_->active_mem      = saved_mem_;
  }

  ActiveMemScope(const ActiveMemScope&)            = delete;
  ActiveMemScope& operator=(const ActiveMemScope&) = delete;

private:
  ContentT* content_;
  ActiveMemMode saved_mode_{ActiveMemMode::none};
  void* saved_mem_{nullptr};
};

/*
 * Fetch the active memory pointer, rejecting the call outright when no scope is
 * open or the open scope came from the wrong provider. Returning a wrong-but-
 * plausible pointer here would corrupt an integrator, so this is a hard error.
 */
template<typename ContentT>
void* require_active_mem(ContentT* content, ActiveMemMode required,
                         const char* what)
{
  if (!content)
  {
    throw std::runtime_error(
      std::string("[sundials4py] the ") + what +
      " callback was invoked on a custom object with no native content");
  }

  if (content->active_mem_mode == ActiveMemMode::none)
  {
    throw std::runtime_error(
      std::string("[sundials4py] the ") + what +
      " callback may only be called while SUNDIALS is inside this solver's "
      "setup() or solve(); no such call is currently active");
  }

  if (content->active_mem_mode != required)
  {
    throw std::runtime_error(
      std::string("[sundials4py] the ") + what + " callback requires " +
      active_mem_mode_name(required) + " memory, but this solver was entered in " +
      active_mem_mode_name(content->active_mem_mode) +
      " mode; the callback and the caller disagree about which memory pointer "
      "is in play");
  }

  return content->active_mem;
}

/*------------------------------------------------------------------------------
 * Optional override detection
 *----------------------------------------------------------------------------*/

/*
 * Report whether `name` is overridden somewhere on the concrete type's MRO.
 *
 * The descriptor found by normal attribute lookup on the concrete type is
 * compared against the descriptor registered on the sundials4py base class. Any
 * difference means a user class -- the concrete class or an intermediate base
 * shared by several concrete classes -- supplied an implementation.
 *
 * Deliberately not hasattr(): the base classes define every abstract method so
 * that calling one raises NotImplementedError with a useful name, which means
 * hasattr() is always true. Deliberately not the instance __dict__ either,
 * since methods live on the type.
 */
template<typename BaseClass>
bool custom_method_overridden(nb::handle impl, const char* name)
{
  nb::object concrete = nb::steal<nb::object>(
    PyObject_GetAttrString((PyObject*)Py_TYPE(impl.ptr()), name));
  if (!concrete.is_valid()) { nb::raise_python_error(); }

  nb::object base = nb::type<BaseClass>().attr(name);
  int same        = PyObject_RichCompareBool(concrete.ptr(), base.ptr(), Py_EQ);
  if (same < 0) { nb::raise_python_error(); }
  return same == 0;
}

/*------------------------------------------------------------------------------
 * Tagged custom content
 *----------------------------------------------------------------------------*/

/*
 * Members every custom content struct is expected to provide, so that the
 * helpers below can be shared without any compile-time member detection (the
 * bindings are built as C++17, where detection idioms are noticeably clumsier
 * than the small cost of a uniform layout).
 *
 * Families that never take ownership of a handle simply leave `strong_impl`
 * empty; families with no native callbacks leave `callbacks` untouched.
 */
#define SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS(TAG)                               \
  static constexpr uintptr_t tag_value = TAG;                                 \
                                                                              \
  /* Distinguishes this private payload from native content routed through  \
     the same wrapper type. Must stay first so the check is cheap. */ \
  uintptr_t tag{tag_value};                                                   \
                                                                              \
  /* Set for handles the user owns: SUNDIALS must not keep the Python       \
     object alive, or attaching a solver would create an immortal cycle. */ \
  nb::object weak_impl;                                                       \
                                                                              \
  /* Set for handles SUNDIALS owns outright (e.g. SUNMatClone results),     \
     where the Python implementation must not be collected first. */ \
  nb::object strong_impl;                                                     \
                                                                              \
  /* Keeps the SUNContext alive for as long as the native shell exists. */    \
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner;            \
                                                                              \
  /* Revocable native callback adapters handed to Python. */                  \
  NativeCallbackRegistry callbacks

/*
 * Recover custom content from a native handle, or nullptr.
 *
 * Every custom content struct starts with a distinctive `tag` constant. Native
 * SUNDIALS objects of the same class route through the same wrapper types, so
 * the tag is what lets a trampoline tell "my content" from "somebody else's
 * content" without ever dereferencing the latter as the wrong type.
 */
template<typename ContentT, typename Handle>
ContentT* custom_content_cast(Handle handle)
{
  if (!handle || !handle->content) { return nullptr; }
  auto* content = static_cast<ContentT*>(handle->content);
  if (content->tag != ContentT::tag_value) { return nullptr; }
  return content;
}

/*
 * Resolve the Python implementation behind a native handle.
 *
 * User-created handles hold a weak reference so that attaching a custom object
 * to an integrator does not create an uncollectable cycle; handles that
 * SUNDIALS itself owns (matrix clones) hold a strong reference so the Python
 * object cannot die first. `label` names the class in error messages.
 */
template<typename ContentT>
nb::object custom_content_impl(ContentT* content, const char* label)
{
  if (!content)
  {
    throw nb::type_error(
      (std::string("native handle does not contain ") + label + " content").c_str());
  }

  if (content->strong_impl.is_valid()) { return content->strong_impl; }

  nb::object impl = content->weak_impl();
  if (impl.is_none())
  {
    throw std::runtime_error(
      (std::string(label) + " Python object has been destroyed while SUNDIALS "
                            "still held its native handle")
        .c_str());
  }
  return impl;
}

/*
 * Destroy custom content, skipping the work entirely if Python is already gone.
 * Centralizes the interpreter-shutdown policy for all four families.
 */
template<typename ContentT>
void custom_content_destroy(ContentT*& content)
{
  if (!content) { return; }

  ShutdownSafeGIL gil;
  if (!gil.python_available())
  {
    /* Interpreter finalization is in progress: releasing the retained Python
       references now is not permitted, so the content is intentionally leaked
       as the process exits. */
    content = nullptr;
    return;
  }

  // Revoke every adapter first: a subclass may still hold Python callables that
  // close over this content's callback state, and they must fail loudly rather
  // than call through freed memory.
  content->callbacks.invalidate_all();
  delete content;
  content = nullptr;
}

} // namespace sundials4py

#endif // _SUNDIALS4PY_CUSTOM_OBJECT_HPP
