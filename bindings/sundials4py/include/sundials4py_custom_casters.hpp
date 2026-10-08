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
 * nanobind type casters that let a Python subclass be passed anywhere the
 * generated bindings expect a raw SUNDIALS pointer.
 *
 * Every sundials4py binding signature names the pointee type (_generic_SUNMatrix
 * and friends), so specializing type_caster for those types is enough to
 * cover the entire generated API without touching a single wrapper. This header
 * must be included by every translation unit that binds such a signature; if one
 * translation unit sees the caster and another does not, the two disagree about
 * how the type is passed. It is therefore pulled in through
 * sundials4py_types.hpp rather than being included piecemeal.
 *----------------------------------------------------------------------------*/

#ifndef SUNDIALS4PY_CUSTOM_CASTERS_HPP
#define SUNDIALS4PY_CUSTOM_CASTERS_HPP

#include <cstdint>
#include <exception>
#include <type_traits>
#include <utility>

#include "sundials_adaptcontroller_custom.hpp"
#include "sundials_domeigestimator_custom.hpp"
#include "sundials_linearsolver_custom.hpp"
#include "sundials_matrix_custom.hpp"
#include "sundials_nonlinearsolver_custom.hpp"
#include "sundials_nvector_custom.hpp"

namespace nanobind::detail {

/*
 * Define the caster for one custom object family.
 *
 * GENERIC_TYPE is the SUNDIALS pointee struct (_generic_SUNMatrix, ...),
 * CUSTOM_CLASS is the sundials4py base class Python subclasses inherit from, and
 * WHAT names the family in the fallback error message.
 *
 * The bodies are identical apart from those three tokens, so they
 * are generated rather than copied: a divergence between them would be a subtle
 * behavioral inconsistency between object families.
 *
 * Deriving from type_caster_base_tag (rather than type_caster_base<Type>) and
 * delegating to a locally constructed type_caster_base keeps the native path
 * bit-for-bit identical to stock nanobind, so pre-existing wrapper objects are
 * unaffected; only the "not a native wrapper" case is new.
 */
#define SUNDIALS4PY_DEFINE_CUSTOM_CASTER(GENERIC_TYPE, CUSTOM_CLASS, WHAT)          \
  template<>                                                                        \
  struct type_caster<GENERIC_TYPE> : type_caster_base_tag                           \
  {                                                                                 \
    using Type                 = GENERIC_TYPE;                                      \
    static constexpr auto Name = const_name<Type>();                                \
    template<typename T>                                                            \
    using Cast = precise_cast_t<T>;                                                 \
                                                                                    \
    /* Nanobind performs two conversion passes only when resolving overloaded \
       functions. Tests that count reported errors should use a binding with a \
       single overload. */     \
    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept     \
    {                                                                               \
      /* Prefer the stock nanobind conversion so existing native wrapper         \
         objects keep their current behavior. */  \
      type_caster_base<Type> native_caster;                                         \
      if (native_caster.from_python(src, flags, cleanup))                           \
      {                                                                             \
        value = native_caster.operator Type*();                                     \
        return true;                                                                \
      }                                                                             \
                                                                                    \
      if (nb::isinstance<sundials4py::CUSTOM_CLASS>(src))                           \
      {                                                                             \
        sundials4py::CUSTOM_CLASS* custom = nullptr;                                \
        try                                                                         \
        {                                                                           \
          /* The returned raw pointer is owned by the Python object's \
             shared_ptr; nanobind only borrows it for this C++ call. */             \
          custom = nb::cast<sundials4py::CUSTOM_CLASS*>(src);                       \
          value  = custom->_get_sundials_handle(src).get();                         \
        }                                                                           \
        catch (const std::exception& error)                                         \
        {                                                                           \
          /* A false result means that this overload cannot accept src, and \
             nanobind will report its normal argument-conversion TypeError. \
             Report the real reason first: it distinguishes a missing required \
             operation from an allocation failure. Route it through the \
             object's own context so a handler installed there sees it; fall \
             back to the global handler only if the cast itself failed. The \
             Python error indicator is left clear. */       \
          ::sundials4py::report_custom_exception(custom ? custom->sunctx().get()    \
                                                        : nullptr,                  \
                                                 WHAT, error, __FILE__,             \
                                                 __LINE__);                         \
          return false;                                                             \
        }                                                                           \
        catch (...)                                                                 \
        {                                                                           \
          ::sundials4py::report_custom_unknown_exception(custom                     \
                                                           ? custom->sunctx().get() \
                                                           : nullptr,               \
                                                         WHAT, __FILE__,            \
                                                         __LINE__);                 \
          return false;                                                             \
        }                                                                           \
        return value != nullptr;                                                    \
      }                                                                             \
                                                                                    \
      return false;                                                                 \
    }                                                                               \
                                                                                    \
    template<typename T>                                                            \
    static handle from_cpp(T&& value, rv_policy policy,                             \
                           cleanup_list* cleanup) noexcept                          \
    {                                                                               \
      Type* ptr = nullptr;                                                          \
      if constexpr (std::is_pointer_v<std::remove_reference_t<T>>)                  \
      {                                                                             \
        ptr = value;                                                                \
      }                                                                             \
      else                                                                          \
      {                                                                             \
        ptr = &value;                                                               \
      }                                                                             \
                                                                                    \
      /* Preserve the Python subclass when C++ passes a borrowed custom handle   \
         back through a callback. Ownership transfers must continue through the  \
         native caster so the handle is released exactly once. */  \
      if (ptr && policy != rv_policy::take_ownership &&                             \
          policy != rv_policy::move)                                                \
      {                                                                             \
        nb::object impl = sundials4py::CUSTOM_CLASS::_python_object_for(ptr);       \
        if (impl.is_valid()) { return impl.release(); }                             \
      }                                                                             \
      return type_caster_base<Type>::from_cpp(std::forward<T>(value), policy,       \
                                              cleanup);                             \
    }                                                                               \
                                                                                    \
    template<typename T_>                                                           \
    bool can_cast() const noexcept                                                  \
    { return std::is_pointer_v<T_> || (value != nullptr); }                         \
                                                                                    \
    operator Type*() { return value; }                                              \
                                                                                    \
    operator Type&()                                                                \
    {                                                                               \
      if (!value) { throw next_overload(); }                                        \
      return *value;                                                                \
    }                                                                               \
                                                                                    \
    operator Type&&()                                                               \
    {                                                                               \
      if (!value) { throw next_overload(); }                                        \
      return (Type&&)*value;                                                        \
    }                                                                               \
                                                                                    \
  private:                                                                          \
    Type* value{nullptr};                                                           \
  };

SUNDIALS4PY_DEFINE_CUSTOM_CASTER(_generic_SUNMatrix, CustomSUNMatrix, "SUNMatrix")

SUNDIALS4PY_DEFINE_CUSTOM_CASTER(_generic_SUNLinearSolver,
                                 CustomSUNLinearSolver, "SUNLinearSolver")

SUNDIALS4PY_DEFINE_CUSTOM_CASTER(_generic_SUNNonlinearSolver,
                                 CustomSUNNonlinearSolver, "SUNNonlinearSolver")

SUNDIALS4PY_DEFINE_CUSTOM_CASTER(_generic_SUNAdaptController,
                                 CustomSUNAdaptController, "SUNAdaptController")

SUNDIALS4PY_DEFINE_CUSTOM_CASTER(SUNDomEigEstimator_, CustomSUNDomEigEstimator,
                                 "SUNDomEigEstimator")

SUNDIALS4PY_DEFINE_CUSTOM_CASTER(_generic_N_Vector, CustomNVector, "N_Vector")

#undef SUNDIALS4PY_DEFINE_CUSTOM_CASTER

} // namespace nanobind::detail

#endif // SUNDIALS4PY_CUSTOM_CASTERS_HPP
