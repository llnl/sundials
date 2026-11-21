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
 *   ActiveMemScope            tracking of the SUNDIALS "mem" pointer
 *   custom_method_overridden  optional-operation detection via the MRO
 *----------------------------------------------------------------------------*/

#ifndef SUNDIALS4PY_CUSTOM_OBJECT_HPP
#define SUNDIALS4PY_CUSTOM_OBJECT_HPP

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_types.h>

#include "sundials4py_core_types.hpp"

#if defined(NB_FREE_THREADED)
#error \
  "sundials4py custom objects are not yet free-threading safe; build without FREE_THREADED"
#endif
namespace sundials4py {

/*------------------------------------------------------------------------------
 * Interpreter shutdown
 *----------------------------------------------------------------------------*/

/*
 * GIL acquisition safe for destructors during interpreter shutdown. If Python
 * is already unavailable, callers must leave their nb::object state untouched.
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

inline PyObject*& pending_custom_exception();
inline bool& capture_custom_exceptions();

/* Record the first exception in the active Python-facing call. The caller must
   hold the GIL. */
inline void record_pending_exception(const std::exception& error) noexcept
{
  if (!capture_custom_exceptions() || pending_custom_exception()) { return; }

  if (const auto* python_error = dynamic_cast<const nb::python_error*>(&error))
  {
    pending_custom_exception() = python_error->value().inc_ref().ptr();
  }
  else if (const auto* builtin_error =
             dynamic_cast<const nb::builtin_exception*>(&error))
  {
    PyObject* type = PyExc_RuntimeError;
    switch (builtin_error->type())
    {
    case nb::exception_type::stop_iteration: type = PyExc_StopIteration; break;
    case nb::exception_type::index_error: type = PyExc_IndexError; break;
    case nb::exception_type::key_error: type = PyExc_KeyError; break;
    case nb::exception_type::value_error: type = PyExc_ValueError; break;
    case nb::exception_type::type_error: type = PyExc_TypeError; break;
    case nb::exception_type::buffer_error: type = PyExc_BufferError; break;
    case nb::exception_type::import_error: type = PyExc_ImportError; break;
    case nb::exception_type::attribute_error:
      type = PyExc_AttributeError;
      break;
    default: break;
    }
    pending_custom_exception() = PyObject_CallFunction(type, "s",
                                                       builtin_error->what());
    if (!pending_custom_exception())
    {
      PyErr_Clear();
      pending_custom_exception() =
        PyObject_CallFunction(PyExc_RuntimeError,
                              "s", "sundials4py: failed to construct the pending exception");
      if (!pending_custom_exception()) { PyErr_Clear(); }
    }
  }
  else
  {
    pending_custom_exception() = PyObject_CallFunction(PyExc_RuntimeError, "s",
                                                       error.what());
    if (!pending_custom_exception())
    {
      PyErr_Clear();
      pending_custom_exception() =
        PyObject_CallFunction(PyExc_RuntimeError,
                              "s", "sundials4py: failed to construct the pending exception");
      if (!pending_custom_exception()) { PyErr_Clear(); }
    }
  }
}

/*
 * Report a C++ exception through the owning SUNContext without allowing a
 * second exception (for example, from a Python error handler) to escape the C
 * callback boundary. nanobind::python_error::what() includes the Python
 * exception type, message, and traceback, so the std::exception overload also
 * preserves the useful Python diagnostic.
 */
inline void report_custom_exception(SUNContext sunctx, const char* operation,
                                    const std::exception& error,
                                    const char* file, int line,
                                    bool record_pending = true) noexcept
{
  try
  {
    ShutdownSafeGIL gil;
    if (!gil.python_available()) { return; }

    /* Preserve the first exception raised during the enclosing SUNDIALS call.
       Scalar and void operations have no reliable error return, so the
       enclosing CustomExceptionScope re-raises this exception after native
       control returns. */
    if (record_pending) { record_pending_exception(error); }

    std::string message = std::string("exception in ") + operation + ":\n" +
                          error.what();
    SUNHandleErrWithMsg(line, operation, file, message.c_str(),
                        SUN_ERR_EXT_FAIL, sunctx);
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

/* One owned exception reference per thread. Only access this with the GIL. */
inline PyObject*& pending_custom_exception()
{
  static thread_local PyObject* exception = nullptr;
  return exception;
}

inline bool& capture_custom_exceptions()
{
  static thread_local bool capture = false;
  return capture;
}

inline bool custom_exception_pending() noexcept
{
  return capture_custom_exceptions() && pending_custom_exception() != nullptr;
}

/*
 * Exception state is scoped in the binding callable, after nanobind has
 * converted all arguments.  A scope saves the enclosing state so nested
 * solver/binding calls neither clear the outer exception nor disable its
 * capture mode.  Only the first exception in one scope is retained by
 * report_custom_exception().
 */
class CustomExceptionScope
{
public:
  CustomExceptionScope()
    : saved_pending_(std::exchange(pending_custom_exception(), nullptr)),
      saved_capture_(std::exchange(capture_custom_exceptions(), true))
  {}

  ~CustomExceptionScope()
  {
    /* If the native call unwound with a C++ exception, no rethrow consumed the
       captured reference.  Drop it before restoring the enclosing scope. */
    Py_XDECREF(std::exchange(pending_custom_exception(), saved_pending_));
    capture_custom_exceptions() = saved_capture_;
  }

  void rethrow_if_pending()
  {
    if (PyObject* exception = std::exchange(pending_custom_exception(), nullptr))
    {
      PyErr_SetRaisedException(exception);
      throw nb::python_error();
    }
  }

  CustomExceptionScope(const CustomExceptionScope&)            = delete;
  CustomExceptionScope& operator=(const CustomExceptionScope&) = delete;

private:
  PyObject* saved_pending_;
  bool saved_capture_;
};

template<typename F, typename R, typename... Args>
auto scoped_impl(F f, R (F::*)(Args...) const)
{
  return [f = std::move(f)](Args... args) -> R
  {
    CustomExceptionScope scope;
    if constexpr (std::is_void_v<R>)
    {
      f(std::forward<Args>(args)...);
      scope.rethrow_if_pending();
    }
    else
    {
      R result = f(std::forward<Args>(args)...);
      scope.rethrow_if_pending();
      return result;
    }
  };
}

template<typename F>
auto scoped(F f)
{
  return scoped_impl(std::move(f), &F::operator());
}

template<typename R, typename... Args>
auto scoped(R (*fn)(Args...))
{
  return scoped([fn](Args... args) -> R
                { return fn(std::forward<Args>(args)...); });
}

/*
 * All generated and handwritten Python functions go through this adapter.
 * Keeping the wrapping at the module-definition boundary ensures argument
 * conversion happens before CustomExceptionScope is constructed.
 */
template<typename Name, typename F, typename... Args>
void scoped_def(nb::module_& module, Name&& name, F&& fn, Args&&... args)
{
  module.def(std::forward<Name>(name), scoped(std::forward<F>(fn)),
             std::forward<Args>(args)...);
}

inline void report_custom_unknown_exception(SUNContext sunctx,
                                            const char* operation,
                                            const char* file, int line,
                                            bool record_pending = true) noexcept
{
  const std::runtime_error error("unknown non-standard C++ exception");
  report_custom_exception(sunctx, operation, error, file, line, record_pending);
}

/* Every trampoline supplies its own failure sentinel, since SUNDIALS callback
   APIs return status codes, pointers, integers, and real scalars. */
#define SUNDIALS4PY_CATCH_AND_REPORT(SUNCTX, OPERATION, FAILURE)               \
  catch (const std::exception& error)                                          \
  {                                                                            \
    ::sundials4py::report_custom_exception(SUNCTX, OPERATION, error, __FILE__, \
                                           __LINE__);                          \
    return FAILURE;                                                            \
  }                                                                            \
  catch (...)                                                                  \
  {                                                                            \
    ::sundials4py::report_custom_unknown_exception(SUNCTX, OPERATION,          \
                                                   __FILE__, __LINE__);        \
    return FAILURE;                                                            \
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
};

/* A closed set makes misspelled registry slots a compile-time error. */
enum class NativeCallbackSlot
{
  atimes,
  rhs,
  psetup,
  psolve,
  sysfn,
  rootsysfn,
  fixedpointsysfn,
  lsetupfn,
  lsolvefn,
  ctestfn,
  normfn,
  getupdatenormfn,
  getconvratefn
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
  static std::shared_ptr<NativeCallbackState<Fn>> prepare(Fn fn, void* data)
  {
    if (!fn) { return nullptr; }

    auto state   = std::make_shared<NativeCallbackState<Fn>>();
    state->fn    = fn;
    state->data  = data;
    state->valid = true;
    return state;
  }

  void commit(NativeCallbackSlot slot,
              std::shared_ptr<NativeCallbackStateBase> state)
  {
    auto& entry = slots_[index(slot)];
    if (entry) { entry->valid = false; }
    entry = std::move(state);
  }

  /* Revoke the state held for one slot, if any. */
  void invalidate(NativeCallbackSlot slot)
  {
    auto& entry = slots_[index(slot)];
    if (entry) { entry->valid = false; }
    entry.reset();
  }

  /* Revoke every slot; called when the custom object is destroyed. */
  void invalidate_all()
  {
    for (auto& entry : slots_)
    {
      if (entry) { entry->valid = false; }
      entry.reset();
    }
  }

private:
  static constexpr std::size_t index(NativeCallbackSlot slot) noexcept
  {
    return static_cast<std::size_t>(slot);
  }

  static constexpr std::size_t slot_count =
    static_cast<std::size_t>(NativeCallbackSlot::getconvratefn) + 1;
  std::array<std::shared_ptr<NativeCallbackStateBase>, slot_count> slots_{};
};

struct PreparedCallback
{
  explicit PreparedCallback(std::shared_ptr<NativeCallbackStateBase> value)
    : state(std::move(value))
  {}

  ~PreparedCallback()
  {
    if (!committed && state) { state->valid = false; }
  }

  void commit(NativeCallbackRegistry& registry, NativeCallbackSlot slot)
  {
    registry.commit(slot, state);
    committed = true;
  }

  std::shared_ptr<NativeCallbackStateBase> state;
  bool committed{false};
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
 * All custom-object state is accessed while holding the GIL, which provides
 * the required synchronization. This includes ActiveMemScope,
 * CustomObjectBase materialization state, NativeCallbackRegistry, the
 * pending-exception slot, and CustomObjectBase::raw_handle_. Do not add locks
 * that remain held across calls into Python: those calls may release the GIL.
 *
 * RAII scope that records the active memory pointer for the duration of one
 * custom setup()/solve() call.
 *
 * Nesting is supported (a package may call setup() from inside solve()) because
 * the previous value is saved and restored, and restoration happens on every
 * exit path including exceptions since it lives in the destructor.
 *
 * Templated on the content type so that each object family can keep its own
 * private content struct; the only requirement is that it expose the members
 * `active_mem` and `active_mem_owner`.
 */
template<typename ContentT>
class ActiveMemScope
{
public:
  ActiveMemScope(ContentT* content, void* mem) : content_(content)
  {
    if (!content_) { return; }
    if (content_->active_mem &&
        content_->active_mem_owner != std::this_thread::get_id())
    {
      throw std::runtime_error(
        "[sundials4py] a custom SUNNonlinearSolver may not be used from two "
        "threads at once");
    }
    saved_mem_                 = content_->active_mem;
    saved_owner_               = content_->active_mem_owner;
    content_->active_mem       = mem;
    content_->active_mem_owner = std::this_thread::get_id();
  }

  ~ActiveMemScope() noexcept
  {
    if (!content_) { return; }
    content_->active_mem       = saved_mem_;
    content_->active_mem_owner = saved_owner_;
  }

  ActiveMemScope(const ActiveMemScope&)            = delete;
  ActiveMemScope& operator=(const ActiveMemScope&) = delete;

private:
  ContentT* content_{nullptr};
  void* saved_mem_{nullptr};
  std::thread::id saved_owner_{};
};

/*
 * Fetch the active memory pointer, rejecting the call when no scope is open or
 * when another thread owns it.
 */
template<typename ContentT>
void* require_active_mem(ContentT* content, const char* what)
{
  if (!content)
  {
    throw std::runtime_error(
      std::string("[sundials4py] the ") + what +
      " callback was invoked on a custom object with no native content");
  }

  if (!content->active_mem)
  {
    throw std::runtime_error(
      std::string("[sundials4py] the ") + what +
      " callback may only be called while SUNDIALS is inside this solver's "
      "setup() or solve(); no such call is currently active");
  }

  if (content->active_mem_owner != std::this_thread::get_id())
  {
    throw std::runtime_error(
      std::string("[sundials4py] the ") + what +
      " callback was invoked from a thread other than the active solver call");
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
  return !concrete.is(base);
}

inline void require_same_context(
  SUNContext expected,
  const std::shared_ptr<std::remove_pointer_t<SUNContext>>& actual,
  const char* operation)
{
  if (actual.get() != expected)
  {
    throw nb::value_error((std::string(operation) +
                           "() must construct its result with the same "
                           "SUNContext as the object being cloned")
                            .c_str());
  }
}

/*------------------------------------------------------------------------------
 * Shared custom-object materialization state and native content
 *----------------------------------------------------------------------------*/

template<typename Derived, typename Pointee>
class CustomObjectBase
{
public:
  using Context = std::shared_ptr<std::remove_pointer_t<SUNContext>>;
  using Handle  = std::shared_ptr<Pointee>;

  CustomObjectBase(Context sunctx, const char* label)
    : sunctx_owner_(std::move(sunctx)), label_(label)
  {}

  Handle _get_sundials_handle(nb::handle self)
  {
    if (raw_handle_)
    {
      return Handle(std::shared_ptr<void>{}, static_cast<Pointee*>(raw_handle_));
    }
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        (std::string(label_) + " base constructor was not initialized").c_str());
    }
    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(std::string("reentrant ") + label_ +
                               " native handle materialization");
    }
    if (state_ == HandleState::materialized) { return handle_; }

    state_ = HandleState::materializing;
    try
    {
      handle_ = static_cast<Derived&>(*this).make_handle(self);
      state_  = HandleState::materialized;
    }
    catch (...)
    {
      handle_.reset();
      state_ = HandleState::unmaterialized;
      throw;
    }
    return handle_;
  }

  bool _is_materialized() const { return state_ == HandleState::materialized; }

  Context sunctx() const { return sunctx_owner_; }

protected:
  Context sunctx_owner_;
  void* raw_handle_{nullptr};

private:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  Handle handle_;
  HandleState state_{HandleState::unmaterialized};
  const char* label_;
};

/* Shared native-content ownership and callback state for every family. */
struct CustomContentBase
{
  nb::object weak_impl;
  nb::object strong_impl;
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner;
  NativeCallbackRegistry callbacks;

  // For C-owned clone shells: where the implementation caches this shell,
  // and the shell itself. Cleared under the GIL in custom_content_destroy.
  void** raw_handle_slot{nullptr};
  void* raw_handle{nullptr};
};

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
 * Centralizes the interpreter-shutdown policy for all custom object families.
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

  if (content->raw_handle_slot && *content->raw_handle_slot == content->raw_handle)
  {
    *content->raw_handle_slot = nullptr;
  }

  // Revoke every adapter first: a subclass may still hold Python callables that
  // close over this content's callback state, and they must fail loudly rather
  // than call through freed memory.
  content->callbacks.invalidate_all();
  delete content;
  content = nullptr;
}

} // namespace sundials4py

#endif // SUNDIALS4PY_CUSTOM_OBJECT_HPP
