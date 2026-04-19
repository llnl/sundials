/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 *
 * This header provides the NVectorBase<Derived> base class that bridges C++
 * vector implementations to the SUNDIALS C NVector interface. The three key
 * aspects of the bridge class design are:
 *
 * 1. The derived class inherits from NVectorBase<Derived>, giving the base
 *    class compile-time knowledge of the concrete type.
 *
 * 2. For each N_Vector operation, the base class defines a static function with
 *    the C signature SUNDIALS expects. Each function extracts derived objects
 *    from the N_Vectors and calls the corresponding derived class method, so
 *    the derived class only uses its own type.
 *
 * 3. Every operation is guarded by a compile-time trait that checks if the
 *    derived class defines a method with the expected name and signature. If
 *    the method is not defined, the corresponding function pointer in the C ops
 *    table stays NULL. This lets the SUNDIALS detect missing operations at
 *    runtime and issue an error, rather than failing at compile.
 *
 * Each of these points are described in more detail below.
 *
 * -----------------------------------------------------------------------------
 * EXAMPLE
 * -----------------------------------------------------------------------------
 *
 * class MyVec : public sundials::NVectorBase<MyVec>
 * {
 * public:
 *   MyVec(int n, SUNContext ctx)
 *     : NVectorBase<MyVec>(ctx), vec_(n, 0.0) {}
 *
 *   MyVec* cloneNew(SUNContext ctx) const
 *   { return new MyVec(vec_.size(), ctx); }
 *
 *   sunindextype length() const { return d_.size(); }
 *
 *   void linearSum(sunrealtype a, const MyVec& x,
 *                  sunrealtype b, const MyVec& y)
 *   {
 *     for (size_t i = 0; i < vec_.size(); ++i)
 *       vec_[i] = a * x.vec_[i] + b * y.vec_[i];
 *   }
 *   // ... remaining operations ...
 * private:
 *   std::vector<sunrealtype> vec_;
 * };
 *
 * -----------------------------------------------------------------------------
 * BASE CLASS DESIGN
 * -----------------------------------------------------------------------------
 *
 * 1. CRTP (Curiously Recurring Template Pattern)
 *
 *    The derived class passes itself as a template parameter to the base class:
 *
 *    class MyVec : public NVectorBase<MyVec> { ... };
 *
 *    This gives the base class compile-time knowledge of the derived type, so
 *    it can cast the NVector void* content pointer to MyVec* and call MyVec
 *    methods without virtual functions or runtime type queries.
 *
 * 2. STATIC TRAMPOLINES (function pointer bridge)
 *
 *    The NVector ops table is a struct of C function pointers, e.g.:
 *
 *    void (*nvscale)(sunrealtype, N_Vector, N_Vector);
 *
 *    A C function pointer cannot point to a C++ non-static member function, so
 *    we define static member functions ("trampolines") with the C signature,
 *    extract the C++ object, and call the corresponding method. For example,
 *    the NVector ops table stores
 *
 *    ops->nvscale = &scale_;
 *
 *    where the static function is defined as
 *
 *    static void scale_(sunrealtype a, N_Vector x, N_Vector z)
 *    {
 *      extract(z)->linearSum(a, *extractConst(x));
 *    }
 *
 *    The convention is to use the output vector as "this" and the input vectors
 *    as const MyVec& arguments.
 *
 *    1D N_Vector array arguments (N_Vector*) are packed as std::vector<MyVec*>
 *    for outputs and std::vector<const MyVec*> for inputs.
 *
 *    2D N_Vector array arguments (N_Vector**) are packed as a vector of
 *    vectors, std::vector<std::vector<...>>, where the outer vector is rows,
 *    the inner vector is columns.
 *
 * 3. SFINAE (Substitution Failure Is Not An Error)
 *
 *    For each operation, a compile-time trait checks if the derived class
 *    defines the corresponding method.
 *
 *    template<class T, class = void>
 *    struct has_foo : std::false_type {};     // default: method missing
 *
 *    template<class T>
 *    struct has_foo<T, std::void_t<decltype(std::declval<T>().foo())>>
 *      : std::true_type {};                   // specialization: method exists
 *
 *    std::declval<T>() pretends to create a T without actually calling any
 *    constructor and decltype(...) asks "what type would this expression
 *    return?".
 *
 *    If T has a foo() method, the expression is valid, void_t produces void,
 *    the specialization matches, and the trait is true.
 *
 *    If T does NOT have foo(), the expression is invalid, but instead of a hard
 *    error, the compiler silently falls back to the primary template (false).
 *
 *    The macros SUNDIALS_DETECT_* macros these traits for the a given trait
 *    name, function signature, and (if applicable) return value. For example,
 *
 *    SUNDIALS_DETECT_METHOD(has_setConst, setConst(sunrealtype{}))
 *
 *    defines the trait has_setConst. The has_* traits are then used in "if
 *    constexpr (condition)" checks that are evaluated at compile time in
 *    installOps() to conditionally install trampolines, e.g.:
 *
 *    if constexpr (impl::has_dotProd<Derived>::value)
 *    {
 *      ops->nvdotprod = &dotprod_;
 *    }
 *
 *    The compiler only compiles the matching branch, so if Derived doesn't
 *    define dotProd, the body is discarded and the trampoline function dotprod_
 *    is never instantiated.
 * ---------------------------------------------------------------------------*/

#ifndef NVECTOR_BASE_HPP
#define NVECTOR_BASE_HPP

#include <cstdio>
#include <type_traits>
#include <utility>
#include <vector>

#include <sundials/sundials_context.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

namespace sundials {

// =============================================================================
// SFINAE detection traits determine at compile time if a type T has a method
// with a specific name, signature, and return type
//
// SUNDIALS_DETECT_METHOD and SUNDIALS_DETECT_CONST_METHOD checks that the call
// expression is valid (used for void-returning methods like linearSum,
// setConst, etc., where the return type doesn't matter). The CONST version is
// used when the expression is called on const T (e.g., length, dotProd, etc.).
//
// SUNDIALS_DETECT_METHOD_R and SUNDIALS_DETECT_CONST_METHOD_R additionally
// check that the return type is convertible to the expected type.
// =============================================================================

namespace impl {

#define SUNDIALS_DETECT_METHOD(trait_name, expr)                      \
  template<class T, class = void>                                     \
  struct trait_name : std::false_type                                 \
  {};                                                                 \
  template<class T>                                                   \
  struct trait_name<T, std::void_t<decltype(std::declval<T>().expr)>> \
    : std::true_type                                                  \
  {};

#define SUNDIALS_DETECT_CONST_METHOD(trait_name, expr)                      \
  template<class T, class = void>                                           \
  struct trait_name : std::false_type                                       \
  {};                                                                       \
  template<class T>                                                         \
  struct trait_name<T, std::void_t<decltype(std::declval<const T>().expr)>> \
    : std::true_type                                                        \
  {};

#define SUNDIALS_DETECT_METHOD_R(trait_name, ret_type, expr)                        \
  template<class T, class = void>                                                   \
  struct trait_name : std::false_type                                               \
  {};                                                                               \
  template<class T>                                                                 \
  struct trait_name<T, std::void_t<decltype(std::declval<T>().expr),                \
                                   std::enable_if_t<std::is_convertible_v<          \
                                     decltype(std::declval<T>().expr), ret_type>>>> \
    : std::true_type                                                                \
  {};

#define SUNDIALS_DETECT_CONST_METHOD_R(trait_name, ret_type, expr)        \
  template<class T, class = void>                                         \
  struct trait_name : std::false_type                                     \
  {};                                                                     \
  template<class T>                                                       \
  struct trait_name<                                                      \
    T, std::void_t<decltype(std::declval<const T>().expr),                \
                   std::enable_if_t<std::is_convertible_v<                \
                     decltype(std::declval<const T>().expr), ret_type>>>> \
    : std::true_type                                                      \
  {};

SUNDIALS_DETECT_CONST_METHOD_R(has_vectorID, N_Vector_ID, vectorID())

SUNDIALS_DETECT_CONST_METHOD_R(has_cloneNew, T*,
                               cloneNew(std::declval<SUNContext>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_length, sunindextype, length())

SUNDIALS_DETECT_CONST_METHOD_R(has_localLength, sunindextype, localLength())

SUNDIALS_DETECT_CONST_METHOD_R(has_communicator, SUNComm, communicator())

// -----------------------------------------------------------------------------
// Standard vector operations
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_METHOD(has_linearSum,
                       linearSum(sunrealtype{}, std::declval<const T&>(),
                                 sunrealtype{}, std::declval<const T&>()))

SUNDIALS_DETECT_METHOD(has_setConst, setConst(sunrealtype{}))

SUNDIALS_DETECT_METHOD(has_prod,
                       prod(std::declval<const T&>(), std::declval<const T&>()))

SUNDIALS_DETECT_METHOD(has_div,
                       div(std::declval<const T&>(), std::declval<const T&>()))

SUNDIALS_DETECT_METHOD(has_scale, scale(sunrealtype{}, std::declval<const T&>()))

SUNDIALS_DETECT_METHOD(has_abs, abs(std::declval<const T&>()))

SUNDIALS_DETECT_METHOD(has_inv, inv(std::declval<const T&>()))

SUNDIALS_DETECT_METHOD(has_addConst,
                       addConst(std::declval<const T&>(), sunrealtype{}))

SUNDIALS_DETECT_CONST_METHOD_R(has_dotProd, sunrealtype,
                               dotProd(std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_maxNorm, sunrealtype, maxNorm())

SUNDIALS_DETECT_CONST_METHOD_R(has_wrmsNorm, sunrealtype,
                               wrmsNorm(std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_wrmsNormMask, sunrealtype,
                               wrmsNormMask(std::declval<const T&>(),
                                            std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_min, sunrealtype, min())

SUNDIALS_DETECT_CONST_METHOD_R(has_wl2Norm, sunrealtype,
                               wl2Norm(std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_l1Norm, sunrealtype, l1Norm())

SUNDIALS_DETECT_METHOD(has_compare,
                       compare(sunrealtype{}, std::declval<const T&>()))

SUNDIALS_DETECT_METHOD_R(has_invTest, sunbooleantype,
                         invTest(std::declval<const T&>()))

SUNDIALS_DETECT_METHOD_R(has_constrMask, sunbooleantype,
                         constrMask(std::declval<const T&>(),
                                    std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_minQuotient, sunrealtype,
                               minQuotient(std::declval<const T&>()))

// -----------------------------------------------------------------------------
// Utility
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_METHOD_R(has_arrayPointer, sunrealtype*, arrayPointer())

SUNDIALS_DETECT_METHOD_R(has_deviceArrayPointer, sunrealtype*,
                         deviceArrayPointer())

SUNDIALS_DETECT_METHOD(has_setArrayPointer,
                       setArrayPointer(std::declval<sunrealtype*>()))

// -----------------------------------------------------------------------------
// Fused vector operations
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_METHOD_R(
  has_linearCombination, SUNErrCode,
  linearCombination(std::declval<const sunrealtype*>(),
                    std::declval<const std::vector<const T*>&>()))

SUNDIALS_DETECT_METHOD_R(has_scaleAddMulti, SUNErrCode,
                         scaleAddMulti(std::declval<const sunrealtype*>(),
                                       std::declval<const std::vector<const T*>&>(),
                                       std::declval<const std::vector<T*>&>()))

SUNDIALS_DETECT_CONST_METHOD_R(
  has_dotProdMulti, SUNErrCode,
  dotProdMulti(std::declval<const std::vector<const T*>&>(),
               std::declval<sunrealtype*>()))

// -----------------------------------------------------------------------------
// Vector array operations
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_METHOD_R(
  has_linearSumVectorArray, SUNErrCode,
  linearSumVectorArray(sunrealtype{}, std::declval<const std::vector<const T*>&>(),
                       sunrealtype{}, std::declval<const std::vector<const T*>&>(),
                       std::declval<const std::vector<T*>&>()))

SUNDIALS_DETECT_METHOD_R(
  has_scaleVectorArray, SUNErrCode,
  scaleVectorArray(std::declval<const sunrealtype*>(),
                   std::declval<const std::vector<const T*>&>(),
                   std::declval<const std::vector<T*>&>()))

SUNDIALS_DETECT_METHOD_R(has_constVectorArray, SUNErrCode,
                         constVectorArray(sunrealtype{},
                                          std::declval<const std::vector<T*>&>()))

SUNDIALS_DETECT_CONST_METHOD_R(
  has_wrmsNormVectorArray, SUNErrCode,
  wrmsNormVectorArray(std::declval<const std::vector<const T*>&>(),
                      std::declval<const std::vector<const T*>&>(),
                      std::declval<sunrealtype*>()))

SUNDIALS_DETECT_CONST_METHOD_R(
  has_wrmsNormMaskVectorArray, SUNErrCode,
  wrmsNormMaskVectorArray(std::declval<const std::vector<const T*>&>(),
                          std::declval<const std::vector<const T*>&>(),
                          std::declval<const T&>(), std::declval<sunrealtype*>()))

SUNDIALS_DETECT_METHOD_R(
  has_scaleAddMultiVectorArray, SUNErrCode,
  scaleAddMultiVectorArray(std::declval<const sunrealtype*>(),
                           std::declval<const std::vector<const T*>&>(),
                           std::declval<const std::vector<std::vector<const T*>>&>(),
                           std::declval<const std::vector<std::vector<T*>>&>()))

SUNDIALS_DETECT_METHOD_R(has_linearCombinationVectorArray, SUNErrCode,
                         linearCombinationVectorArray(
                           std::declval<const sunrealtype*>(),
                           std::declval<const std::vector<std::vector<const T*>>&>(),
                           std::declval<const std::vector<T*>&>()))

// -----------------------------------------------------------------------------
// Local reduction kernels (no parallel communication)
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_CONST_METHOD_R(has_dotProdLocal, sunrealtype,
                               dotProdLocal(std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_maxNormLocal, sunrealtype, maxNormLocal())

SUNDIALS_DETECT_CONST_METHOD_R(has_minLocal, sunrealtype, minLocal())

SUNDIALS_DETECT_CONST_METHOD_R(has_l1NormLocal, sunrealtype, l1NormLocal())

SUNDIALS_DETECT_METHOD_R(has_invTestLocal, sunbooleantype,
                         invTestLocal(std::declval<const T&>()))

SUNDIALS_DETECT_METHOD_R(has_constrMaskLocal, sunbooleantype,
                         constrMaskLocal(std::declval<const T&>(),
                                         std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_minQuotientLocal, sunrealtype,
                               minQuotientLocal(std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_wsqrSumLocal, sunrealtype,
                               wsqrSumLocal(std::declval<const T&>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_wsqrSumMaskLocal, sunrealtype,
                               wsqrSumMaskLocal(std::declval<const T&>(),
                                                std::declval<const T&>()))

// -----------------------------------------------------------------------------
// Single buffer reduction
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_CONST_METHOD_R(
  has_dotProdMultiLocal, SUNErrCode,
  dotProdMultiLocal(std::declval<const std::vector<const T*>&>(),
                    std::declval<sunrealtype*>()))

SUNDIALS_DETECT_CONST_METHOD_R(has_dotProdMultiAllReduce, SUNErrCode,
                               dotProdMultiAllReduce(int{},
                                                     std::declval<sunrealtype*>()))

// -----------------------------------------------------------------------------
// XBraid interface
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_CONST_METHOD(has_bufSize, bufSize(std::declval<sunindextype*>()))

SUNDIALS_DETECT_CONST_METHOD(has_bufPack, bufPack(std::declval<void*>()))

SUNDIALS_DETECT_METHOD(has_bufUnpack, bufUnpack(std::declval<void*>()))

// -----------------------------------------------------------------------------
// Debugging
// -----------------------------------------------------------------------------

SUNDIALS_DETECT_CONST_METHOD(has_print, print())

SUNDIALS_DETECT_CONST_METHOD(has_printFile, printFile(std::declval<FILE*>()))

} // namespace impl

// =============================================================================
// NVectorBase<Derived>
//
// CRTP base class. Manages the lifetime of the N_Vector.
// =============================================================================

template<class Derived>
class NVectorBase
{
public:
  // ---------------------------------------------------------------------------
  // Construction
  //
  // Allocates the N_Vector with N_VNewEmpty (nullifies the content and all ops
  // pointers), sets the content pointer to this object, and calls installOps()
  // to set the trampolines for every operation the derived class provides.
  // ---------------------------------------------------------------------------

  explicit NVectorBase(SUNContext sunctx)
  {
    nv_ = N_VNewEmpty(sunctx);

    // Store a pointer to this C++ object in the C struct's void* content
    // field. The double static_cast (this -> Derived* -> void*) ensures
    // we store the Derived pointer, not a base-class pointer - this matters
    // if multiple inheritance shifts the address. The trampolines will
    // reverse this cast: static_cast<Derived*>(nv_->content).
    nv_->content = static_cast<void*>(static_cast<Derived*>(this));

    // Check the Derived class for methods and install trampolines for every
    // method that exists
    installOps();
  }

  // ---------------------------------------------------------------------------
  // Move construction and assignment
  //
  // Transfers ownership of the C struct. The source is left in a valid but
  // empty state (nv_ == nullptr). The content pointer is updated to point at
  // the new location of the C++ object.
  // ---------------------------------------------------------------------------

  NVectorBase(NVectorBase&& other) noexcept : nv_(other.nv_)
  {
    other.nv_ = nullptr;
    if (nv_) nv_->content = static_cast<void*>(static_cast<Derived*>(this));
  }

  NVectorBase& operator=(NVectorBase&& other) noexcept
  {
    if (this != &other)
    {
      if (nv_) N_VFreeEmpty(nv_);
      nv_       = other.nv_;
      other.nv_ = nullptr;
      if (nv_) nv_->content = static_cast<void*>(static_cast<Derived*>(this));
    }
    return *this;
  }

  // ---------------------------------------------------------------------------
  // Copy construction and assignment
  //
  // Copy construction allocates a fresh C struct (via N_VNewEmpty) for the new
  // object. The ops table is populated identically since it is the same Derived
  // type. The derived class copies its own data through its normal
  // compiler-generated (or user-defined) copy constructor.
  //
  // Copy assignment is a no-op at the base level: the target already owns a
  // valid C struct with the correct ops. The derived class's copy assignment
  // handles data member copying as usual.
  // ---------------------------------------------------------------------------

  NVectorBase(const NVectorBase& other)
  {
    nv_          = N_VNewEmpty(other.nv_->sunctx);
    nv_->content = static_cast<void*>(static_cast<Derived*>(this));
    installOps();
  }

  NVectorBase& operator=(const NVectorBase&) { return *this; }

  // ---------------------------------------------------------------------------
  // Destruction
  //
  // Detaches the content and ops pointers (to prevent N_VFreeEmpty from freeing
  // memory that belongs to the C++ object or was allocated by N_VNewEmpty as
  // part of the ops table) and then frees the C struct.
  //
  // For clones created by the clone trampoline, destruction happens via
  // N_VDestroy -> destroy_ trampoline -> delete, which nulls nv_ before
  // reaching this destructor, so the body is safely skipped.
  // ---------------------------------------------------------------------------

  virtual ~NVectorBase()
  {
    if (nv_)
    {
      nv_->content = nullptr;
      N_VFreeEmpty(nv_);
      nv_ = nullptr;
    }
  }

  // ---------------------------------------------------------------------------
  // C handle access
  //
  // Returns the underlying N_Vector that SUNDIALS solvers operate on. The
  // implicit conversion operators allow passing a NVectorBase-derived object
  // directly to any C function expecting an N_Vector.
  // ---------------------------------------------------------------------------

  N_Vector getNVector() noexcept { return nv_; }

  N_Vector getNVector() const noexcept { return nv_; }

  operator N_Vector() noexcept { return nv_; }

  operator N_Vector() const noexcept { return nv_; }

  // ---------------------------------------------------------------------------
  // Extraction helpers
  //
  // Recover a typed Derived pointer from the void* content field of an
  // N_Vector. These are public so that user code outside the class hierarchy
  // can use them e.g., in callback functions that receive N_Vector arguments
  // and need to access the underlying C++ data (ODE right-hand sides, Jacobian
  // evaluations, preconditioner solves, etc.):
  //
  // int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
  // {
  //   auto* yv    = MyVec::extract(y);
  //   auto* ydotv = MyVec::extract(ydot);
  //   // ... work directly with MyVec members ...
  // }
  // ---------------------------------------------------------------------------

  // Extract a mutable Derived pointer from an N_Vector.
  static Derived* extract(N_Vector v)
  { return static_cast<Derived*>(v->content); }

  // Extract a const Derived pointer from an N_Vector.
  static const Derived* extractConst(N_Vector v)
  { return static_cast<const Derived*>(v->content); }

private:
  // The underlying C-level N_Vector struct.
  N_Vector nv_ = nullptr;

  // ===========================================================================
  // Vector array extraction helpers
  //
  // Convert N_Vector* (1D) and N_Vector** (2D) arrays into std::vectors of
  // typed Derived pointers. These are called inside the trampolines for fused
  // and vector-array operations. This assumes the allocation cost is negligible
  // compared to the work these operations perform on the vector data.
  // ===========================================================================

  // Extract an N_Vector* array into a vector of const Derived pointers.
  static std::vector<const Derived*> extractConstVec(N_Vector* nvs, int n)
  {
    std::vector<const Derived*> v(n);
    for (int i = 0; i < n; ++i) { v[i] = extractConst(nvs[i]); }
    return v;
  }

  // Extract an N_Vector* array into a vector of mutable Derived pointers.
  static std::vector<Derived*> extractMutVec(N_Vector* nvs, int n)
  {
    std::vector<Derived*> v(n);
    for (int i = 0; i < n; ++i) { v[i] = extract(nvs[i]); }
    return v;
  }

  // Extract an N_Vector** (2D) array into a vector of vectors of const Derived
  // pointers, nvs[j][i] -> result[j][i].
  static std::vector<std::vector<const Derived*>> extractConstVec2D(N_Vector** nvs,
                                                                    int nrows,
                                                                    int ncols)
  {
    std::vector<std::vector<const Derived*>> v(nrows);
    for (int j = 0; j < nrows; ++j)
    {
      v[j].resize(ncols);
      for (int i = 0; i < ncols; ++i) { v[j][i] = extractConst(nvs[j][i]); }
    }
    return v;
  }

  // Extract an N_Vector** (2D) array into a vector of vectors of mutable
  // Derived pointers, nvs[j][i] -> result[j][i].
  static std::vector<std::vector<Derived*>> extractMutVec2D(N_Vector** nvs,
                                                            int nrows, int ncols)
  {
    std::vector<std::vector<Derived*>> v(nrows);
    for (int j = 0; j < nrows; ++j)
    {
      v[j].resize(ncols);
      for (int i = 0; i < ncols; ++i) { v[j][i] = extract(nvs[j][i]); }
    }
    return v;
  }

  // ===========================================================================
  // Static trampoline functions
  //
  // These are the glue between SUNDIALS C code and the derived C++ class. Each
  // trampoline is a static function (so it has a plain C-compatible function
  // pointer) with the signature SUNDIALS expects.
  //
  // The full call chain for a typical operation (e.g., dot product):
  //
  // 1. User code:     sunrealtype r = N_VDotProd(x, y);
  // 2. SUNDIALS C:    x->ops->nvdotprod(x, y)
  // 3. Ops table:     nvdotprod points to &NVectorBase::dotprod_
  // 4. Trampoline:    dotprod_(x, y) does:
  //                     extractConst(x)->dotProd(*extractConst(y))
  //                   which is:
  //                     static_cast<const Derived*>(x->content)
  //                       ->dotProd(*static_cast<const Derived*>(y->content))
  // 5. Derived class: MyVec::dotProd(const MyVec& y) executes
  //
  // Conventions:
  //
  // - Input-only vectors -> const Derived& (via extractConst + dereference)
  // - Output vector      -> Derived*       (via extract, becomes "this")
  // ===========================================================================

  // Returns the vector type identifier.
  static N_Vector_ID getvectorid_(N_Vector v) { return SUNDIALS_NVEC_CUSTOM; }

  // Allocates a new Derived object by calling cloneNew() on the source. The
  // returned N_Vector is owned by SUNDIALS and will be freed via N_VDestroy ->
  // destroy_.
  static N_Vector clone_(N_Vector w)
  { return extractConst(w)->cloneNew(w->sunctx)->getNVector(); }

  // Frees a clone created by clone_. Detaches nv_ so the destructor does not
  // double-free, calls N_VFreeEmpty to release the C struct, then deletes the
  // C++ object.
  static void destroy_(N_Vector v)
  {
    Derived* self      = extract(v);
    self->nv_->content = nullptr;
    self->nv_          = nullptr;
    N_VFreeEmpty(v);
    delete self;
  }

  // Returns the global length of the vector.
  static sunindextype getlength_(N_Vector v)
  { return extractConst(v)->length(); }

  // Returns the local (per-process) length.
  static sunindextype getlocallength_(N_Vector v)
  { return extractConst(v)->localLength(); }

  // Returns the MPI communicator (or SUN_COMM_NULL).
  static SUNComm getcommunicator_(N_Vector v)
  { return extractConst(v)->communicator(); }

  // ---------------------------------------------------------------------------
  // Standard vector operation trampolines
  // ---------------------------------------------------------------------------

  // z = a*x + b*y  (dispatched on z)
  static void linearsum_(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                         N_Vector z)
  { extract(z)->linearSum(a, *extractConst(x), b, *extractConst(y)); }

  // z[i] = c for all i
  static void const_(sunrealtype c, N_Vector z) { extract(z)->setConst(c); }

  // z = x .* y  (element-wise product)
  static void prod_(N_Vector x, N_Vector y, N_Vector z)
  { extract(z)->prod(*extractConst(x), *extractConst(y)); }

  // z = x ./ y  (element-wise division)
  static void div_(N_Vector x, N_Vector y, N_Vector z)
  { extract(z)->div(*extractConst(x), *extractConst(y)); }

  // z = c * x
  static void scale_(sunrealtype c, N_Vector x, N_Vector z)
  { extract(z)->scale(c, *extractConst(x)); }

  // z = |x|
  static void abs_(N_Vector x, N_Vector z)
  { extract(z)->abs(*extractConst(x)); }

  // z = 1/x
  static void inv_(N_Vector x, N_Vector z)
  { extract(z)->inv(*extractConst(x)); }

  // z[i] = x[i] + b for all i
  static void addconst_(N_Vector x, sunrealtype b, N_Vector z)
  { extract(z)->addConst(*extractConst(x), b); }

  // Returns x . y  (dot product)
  static sunrealtype dotprod_(N_Vector x, N_Vector y)
  { return extractConst(x)->dotProd(*extractConst(y)); }

  // Returns max_i |x_i|  (infinity norm)
  static sunrealtype maxnorm_(N_Vector x) { return extractConst(x)->maxNorm(); }

  // Returns the weighted root-mean-square norm of x with weights w.
  static sunrealtype wrmsnorm_(N_Vector x, N_Vector w)
  { return extractConst(x)->wrmsNorm(*extractConst(w)); }

  // Returns the weighted RMS norm of x with weights w, masked by id.
  static sunrealtype wrmsnormmask_(N_Vector x, N_Vector w, N_Vector id)
  { return extractConst(x)->wrmsNormMask(*extractConst(w), *extractConst(id)); }

  // Returns min_i x_i
  static sunrealtype min_(N_Vector x) { return extractConst(x)->min(); }

  // Returns the weighted L2 norm of x with weights w.
  static sunrealtype wl2norm_(N_Vector x, N_Vector w)
  { return extractConst(x)->wl2Norm(*extractConst(w)); }

  // Returns the L1 norm of x.
  static sunrealtype l1norm_(N_Vector x) { return extractConst(x)->l1Norm(); }

  // z_i = |x_i| >= c ? 1 : 0
  static void compare_(sunrealtype c, N_Vector x, N_Vector z)
  { extract(z)->compare(c, *extractConst(x)); }

  // z = 1/x, returns false if any x_i == 0.
  static sunbooleantype invtest_(N_Vector x, N_Vector z)
  { return extract(z)->invTest(*extractConst(x)); }

  // Constraint mask check. m is the output mask vector.
  static sunbooleantype constrmask_(N_Vector c, N_Vector x, N_Vector m)
  { return extract(m)->constrMask(*extractConst(c), *extractConst(x)); }

  // Returns min_i (num_i / denom_i) where denom_i != 0.
  static sunrealtype minquotient_(N_Vector num, N_Vector denom)
  { return extractConst(num)->minQuotient(*extractConst(denom)); }

  // ---------------------------------------------------------------------------
  // Utility trampolines
  // ---------------------------------------------------------------------------

  // Returns a pointer to the host data array.
  static sunrealtype* getarraypointer_(N_Vector v)
  { return extract(v)->arrayPointer(); }

  // Returns a pointer to the device data array.
  static sunrealtype* getdevicearraypointer_(N_Vector v)
  { return extract(v)->deviceArrayPointer(); }

  // Sets the host data array pointer.
  static void setarraypointer_(sunrealtype* d, N_Vector v)
  { extract(v)->setArrayPointer(d); }

  // ---------------------------------------------------------------------------
  // Fused vector operation trampolines
  // ---------------------------------------------------------------------------

  // z = sum_j c[j]*X[j]  (linear combination, dispatched on z)
  static SUNErrCode linearcombination_(int nvec, sunrealtype* c, N_Vector* X,
                                       N_Vector z)
  {
    auto Xv = extractConstVec(X, nvec);
    return extract(z)->linearCombination(c, Xv);
  }

  // Z[j] = a[j]*x + Y[j]  (scale-add multi, dispatched on x)
  static SUNErrCode scaleaddmulti_(int nvec, sunrealtype* a, N_Vector x,
                                   N_Vector* Y, N_Vector* Z)
  {
    auto Yv = extractConstVec(Y, nvec);
    auto Zv = extractMutVec(Z, nvec);
    return extract(x)->scaleAddMulti(a, Yv, Zv);
  }

  // dotprods[j] = x . Y[j]  (dispatched on x)
  static SUNErrCode dotprodmulti_(int nvec, N_Vector x, N_Vector* Y,
                                  sunrealtype* dotprods)
  {
    auto Yv = extractConstVec(Y, nvec);
    return extractConst(x)->dotProdMulti(Yv, dotprods);
  }

  // ---------------------------------------------------------------------------
  // Vector array operation trampolines
  //
  // These have no single natural "this". The trampoline dispatches on an
  // arbitrary element (typically Z[0] for output arrays, X[0] for input-only
  // operations).
  // ---------------------------------------------------------------------------

  // Z[k] = a*X[k] + b*Y[k] for all k
  static SUNErrCode linearsumvectorarray_(int nvec, sunrealtype a, N_Vector* X,
                                          sunrealtype b, N_Vector* Y, N_Vector* Z)
  {
    auto Xv = extractConstVec(X, nvec);
    auto Yv = extractConstVec(Y, nvec);
    auto Zv = extractMutVec(Z, nvec);
    return extract(Z[0])->linearSumVectorArray(a, Xv, b, Yv, Zv);
  }

  // Z[k] = c[k]*X[k] for all k
  static SUNErrCode scalevectorarray_(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector* Z)
  {
    auto Xv = extractConstVec(X, nvec);
    auto Zv = extractMutVec(Z, nvec);
    return extract(Z[0])->scaleVectorArray(c, Xv, Zv);
  }

  // Z[k] = c for all k
  static SUNErrCode constvectorarray_(int nvec, sunrealtype c, N_Vector* Z)
  {
    auto Zv = extractMutVec(Z, nvec);
    return extract(Z[0])->constVectorArray(c, Zv);
  }

  // nrm[k] = wrmsNorm(X[k], W[k]) for all k
  static SUNErrCode wrmsnormvectorarray_(int nvec, N_Vector* X, N_Vector* W,
                                         sunrealtype* nrm)
  {
    auto Xv = extractConstVec(X, nvec);
    auto Wv = extractConstVec(W, nvec);
    return extractConst(X[0])->wrmsNormVectorArray(Xv, Wv, nrm);
  }

  // nrm[k] = wrmsNormMask(X[k], W[k], id) for all k
  static SUNErrCode wrmsnormmaskvectorarray_(int nvec, N_Vector* X, N_Vector* W,
                                             N_Vector id, sunrealtype* nrm)
  {
    auto Xv = extractConstVec(X, nvec);
    auto Wv = extractConstVec(W, nvec);
    return extractConst(X[0])->wrmsNormMaskVectorArray(Xv, Wv,
                                                       *extractConst(id), nrm);
  }

  // Z[j][i] = a[j]*X[i] + Y[j][i]
  static SUNErrCode scaleaddmultivectorarray_(int nvec, int nsum,
                                              sunrealtype* a, N_Vector* X,
                                              N_Vector** Y, N_Vector** Z)
  {
    auto Xv = extractConstVec(X, nvec);
    auto Yv = extractConstVec2D(Y, nsum, nvec);
    auto Zv = extractMutVec2D(Z, nsum, nvec);
    return extract(X[0])->scaleAddMultiVectorArray(a, Xv, Yv, Zv);
  }

  // Z[i] = sum_j c[j]*X[j][i]
  static SUNErrCode linearcombinationvectorarray_(int nvec, int nsum,
                                                  sunrealtype* c, N_Vector** X,
                                                  N_Vector* Z)
  {
    auto Xv = extractConstVec2D(X, nsum, nvec);
    auto Zv = extractMutVec(Z, nvec);
    return extract(Z[0])->linearCombinationVectorArray(c, Xv, Zv);
  }

  // ---------------------------------------------------------------------------
  // Local reduction kernel trampolines
  // ---------------------------------------------------------------------------

  static sunrealtype dotprodlocal_(N_Vector x, N_Vector y)
  { return extractConst(x)->dotProdLocal(*extractConst(y)); }

  static sunrealtype maxnormlocal_(N_Vector x)
  { return extractConst(x)->maxNormLocal(); }

  static sunrealtype minlocal_(N_Vector x)
  { return extractConst(x)->minLocal(); }

  static sunrealtype l1normlocal_(N_Vector x)
  { return extractConst(x)->l1NormLocal(); }

  static sunbooleantype invtestlocal_(N_Vector x, N_Vector z)
  { return extract(z)->invTestLocal(*extractConst(x)); }

  static sunbooleantype constrmasklocal_(N_Vector c, N_Vector x, N_Vector m)
  { return extract(m)->constrMaskLocal(*extractConst(c), *extractConst(x)); }

  static sunrealtype minquotientlocal_(N_Vector num, N_Vector denom)
  { return extractConst(num)->minQuotientLocal(*extractConst(denom)); }

  static sunrealtype wsqrsumlocal_(N_Vector x, N_Vector w)
  { return extractConst(x)->wsqrSumLocal(*extractConst(w)); }

  static sunrealtype wsqrsummasklocal_(N_Vector x, N_Vector w, N_Vector id)
  {
    return extractConst(x)->wsqrSumMaskLocal(*extractConst(w), *extractConst(id));
  }

  // ---------------------------------------------------------------------------
  // Single buffer reduction trampolines
  // ---------------------------------------------------------------------------

  // dotprods[j] = x . Y[j] (local, no communication)
  static SUNErrCode dotprodmultilocal_(int nvec, N_Vector x, N_Vector* Y,
                                       sunrealtype* dotprods)
  {
    auto Yv = extractConstVec(Y, nvec);
    return extractConst(x)->dotProdMultiLocal(Yv, dotprods);
  }

  // All-reduce the partial sum array. x provides the communicator.
  static SUNErrCode dotprodmultiallreduce_(int nvec_total, N_Vector x,
                                           sunrealtype* sum)
  { return extractConst(x)->dotProdMultiAllReduce(nvec_total, sum); }

  // ---------------------------------------------------------------------------
  // XBraid interface trampolines
  // ---------------------------------------------------------------------------

  static SUNErrCode bufsize_(N_Vector x, sunindextype* s)
  { return extractConst(x)->bufSize(s); }

  static SUNErrCode bufpack_(N_Vector x, void* b)
  { return extractConst(x)->bufPack(b); }

  static SUNErrCode bufunpack_(N_Vector x, void* b)
  { return extract(x)->bufUnpack(b); }

  // ---------------------------------------------------------------------------
  // Print trampolines
  // ---------------------------------------------------------------------------

  static void print_(N_Vector v) { extractConst(v)->print(); }

  static void printfile_(N_Vector v, FILE* f) { extractConst(v)->printFile(f); }

  // ===========================================================================
  // installOps
  //
  // Populates the C ops table. Setting an operation is guarded by if constexpr
  // with the corresponding SFINAE trait. If the derived class does not define a
  // method, the trampoline is never instantiated and the ops pointer stays at
  // the NULL value set by N_VNewEmpty. This allows the SUNDIALS to detect
  // missing operations at runtime and issue appropriate error messages.
  // ===========================================================================

  void installOps()
  {
    auto* ops = nv_->ops;

    ops->nvgetvectorid = &getvectorid_;

    if constexpr (impl::has_cloneNew<Derived>::value)
    {
      ops->nvclone   = &clone_;
      ops->nvdestroy = &destroy_;
    }

    if constexpr (impl::has_length<Derived>::value)
    {
      ops->nvgetlength = &getlength_;
    }

    if constexpr (impl::has_localLength<Derived>::value)
    {
      ops->nvgetlocallength = &getlocallength_;
    }

    if constexpr (impl::has_communicator<Derived>::value)
    {
      ops->nvgetcommunicator = &getcommunicator_;
    }

    // -- Standard vector operations -----------------------------------------

    if constexpr (impl::has_linearSum<Derived>::value)
    {
      ops->nvlinearsum = &linearsum_;
    }

    if constexpr (impl::has_setConst<Derived>::value)
    {
      ops->nvconst = &const_;
    }

    if constexpr (impl::has_prod<Derived>::value) { ops->nvprod = &prod_; }

    if constexpr (impl::has_div<Derived>::value) { ops->nvdiv = &div_; }

    if constexpr (impl::has_scale<Derived>::value) { ops->nvscale = &scale_; }

    if constexpr (impl::has_abs<Derived>::value) { ops->nvabs = &abs_; }

    if constexpr (impl::has_inv<Derived>::value) { ops->nvinv = &inv_; }

    if constexpr (impl::has_addConst<Derived>::value)
    {
      ops->nvaddconst = &addconst_;
    }

    if constexpr (impl::has_dotProd<Derived>::value)
    {
      ops->nvdotprod = &dotprod_;
    }

    if constexpr (impl::has_maxNorm<Derived>::value)
    {
      ops->nvmaxnorm = &maxnorm_;
    }

    if constexpr (impl::has_wrmsNorm<Derived>::value)
    {
      ops->nvwrmsnorm = &wrmsnorm_;
    }

    if constexpr (impl::has_wrmsNormMask<Derived>::value)
    {
      ops->nvwrmsnormmask = &wrmsnormmask_;
    }

    if constexpr (impl::has_min<Derived>::value) { ops->nvmin = &min_; }

    if constexpr (impl::has_wl2Norm<Derived>::value)
    {
      ops->nvwl2norm = &wl2norm_;
    }

    if constexpr (impl::has_l1Norm<Derived>::value)
    {
      ops->nvl1norm = &l1norm_;
    }

    if constexpr (impl::has_compare<Derived>::value)
    {
      ops->nvcompare = &compare_;
    }

    if constexpr (impl::has_invTest<Derived>::value)
    {
      ops->nvinvtest = &invtest_;
    }

    if constexpr (impl::has_constrMask<Derived>::value)
    {
      ops->nvconstrmask = &constrmask_;
    }

    if constexpr (impl::has_minQuotient<Derived>::value)
    {
      ops->nvminquotient = &minquotient_;
    }

    // -- Utility ------------------------------------------------------------

    if constexpr (impl::has_arrayPointer<Derived>::value)
    {
      ops->nvgetarraypointer = &getarraypointer_;
    }

    if constexpr (impl::has_deviceArrayPointer<Derived>::value)
    {
      ops->nvgetdevicearraypointer = &getdevicearraypointer_;
    }

    if constexpr (impl::has_setArrayPointer<Derived>::value)
    {
      ops->nvsetarraypointer = &setarraypointer_;
    }

    // -- Fused vector operations --------------------------------------------

    if constexpr (impl::has_linearCombination<Derived>::value)
    {
      ops->nvlinearcombination = &linearcombination_;
    }

    if constexpr (impl::has_scaleAddMulti<Derived>::value)
    {
      ops->nvscaleaddmulti = &scaleaddmulti_;
    }

    if constexpr (impl::has_dotProdMulti<Derived>::value)
    {
      ops->nvdotprodmulti = &dotprodmulti_;
    }

    // -- Vector array operations --------------------------------------------

    if constexpr (impl::has_linearSumVectorArray<Derived>::value)
    {
      ops->nvlinearsumvectorarray = &linearsumvectorarray_;
    }

    if constexpr (impl::has_scaleVectorArray<Derived>::value)
    {
      ops->nvscalevectorarray = &scalevectorarray_;
    }

    if constexpr (impl::has_constVectorArray<Derived>::value)
    {
      ops->nvconstvectorarray = &constvectorarray_;
    }

    if constexpr (impl::has_wrmsNormVectorArray<Derived>::value)
    {
      ops->nvwrmsnormvectorarray = &wrmsnormvectorarray_;
    }

    if constexpr (impl::has_wrmsNormMaskVectorArray<Derived>::value)
    {
      ops->nvwrmsnormmaskvectorarray = &wrmsnormmaskvectorarray_;
    }

    if constexpr (impl::has_scaleAddMultiVectorArray<Derived>::value)
    {
      ops->nvscaleaddmultivectorarray = &scaleaddmultivectorarray_;
    }

    if constexpr (impl::has_linearCombinationVectorArray<Derived>::value)
    {
      ops->nvlinearcombinationvectorarray = &linearcombinationvectorarray_;
    }

    // -- Local reduction kernels --------------------------------------------

    if constexpr (impl::has_dotProdLocal<Derived>::value)
    {
      ops->nvdotprodlocal = &dotprodlocal_;
    }

    if constexpr (impl::has_maxNormLocal<Derived>::value)
    {
      ops->nvmaxnormlocal = &maxnormlocal_;
    }

    if constexpr (impl::has_minLocal<Derived>::value)
    {
      ops->nvminlocal = &minlocal_;
    }

    if constexpr (impl::has_l1NormLocal<Derived>::value)
    {
      ops->nvl1normlocal = &l1normlocal_;
    }

    if constexpr (impl::has_invTestLocal<Derived>::value)
    {
      ops->nvinvtestlocal = &invtestlocal_;
    }

    if constexpr (impl::has_constrMaskLocal<Derived>::value)
    {
      ops->nvconstrmasklocal = &constrmasklocal_;
    }

    if constexpr (impl::has_minQuotientLocal<Derived>::value)
    {
      ops->nvminquotientlocal = &minquotientlocal_;
    }

    if constexpr (impl::has_wsqrSumLocal<Derived>::value)
    {
      ops->nvwsqrsumlocal = &wsqrsumlocal_;
    }

    if constexpr (impl::has_wsqrSumMaskLocal<Derived>::value)
    {
      ops->nvwsqrsummasklocal = &wsqrsummasklocal_;
    }

    // -- Single buffer reduction --------------------------------------------

    if constexpr (impl::has_dotProdMultiLocal<Derived>::value)
    {
      ops->nvdotprodmultilocal = &dotprodmultilocal_;
    }

    if constexpr (impl::has_dotProdMultiAllReduce<Derived>::value)
    {
      ops->nvdotprodmultiallreduce = &dotprodmultiallreduce_;
    }

    // -- XBraid interface ---------------------------------------------------

    if constexpr (impl::has_bufSize<Derived>::value)
    {
      ops->nvbufsize = &bufsize_;
    }

    if constexpr (impl::has_bufPack<Derived>::value)
    {
      ops->nvbufpack = &bufpack_;
    }

    if constexpr (impl::has_bufUnpack<Derived>::value)
    {
      ops->nvbufunpack = &bufunpack_;
    }

    // -- Debug --------------------------------------------------------------

    if constexpr (impl::has_print<Derived>::value) { ops->nvprint = &print_; }

    if constexpr (impl::has_printFile<Derived>::value)
    {
      ops->nvprintfile = &printfile_;
    }
  }
};

} // namespace sundials

#endif // NVECTOR_BASE_HPP
