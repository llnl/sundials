/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_TPETRA_VECTOR_KERNELS_HPP_
#define _SUNDIALS_TPETRA_VECTOR_KERNELS_HPP_

#include <Kokkos_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Trilinos_version.h>
#include <nvector/trilinos/SundialsTpetraVectorInterface.hpp>

#if TRILINOS_MAJOR_VERSION > 13
#include <Tpetra_Access.hpp>
#endif

namespace sundials {
namespace trilinos {

/**
   * The namespace contains custom Kokkos-based kernels needed by SUNDIALS
   *
   * Kernels are inlined in case this file is included in more than one
   * translation unit.
   */
namespace nvector_tpetra {
using Teuchos::outArg;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::REDUCE_SUM;
using Teuchos::reduceAll;

typedef sundials::trilinos::nvector_tpetra::TpetraVectorInterface::vector_type vector_type;
typedef vector_type::scalar_type scalar_type;
typedef vector_type::mag_type mag_type;
typedef vector_type::global_ordinal_type global_ordinal_type;
typedef vector_type::local_ordinal_type local_ordinal_type;
typedef vector_type::node_type::memory_space memory_space;
typedef vector_type::execution_space execution_space;

/*----------------------------------------------------------------
 *  Streaming vector kernels
 *---------------------------------------------------------------*/

/// Divide: z(i) = x(i)/y(i) forall i
inline void elementWiseDivide(const vector_type& x, const vector_type& y,
                              vector_type& z)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (y.need_sync<memory_space>())
    const_cast<vector_type&>(y).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto y_2d = y.getLocalView<memory_space>();
  auto y_1d = Kokkos::subview(y_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>();
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);

  z.modify<memory_space>();
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto y_2d = y.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto y_1d = Kokkos::subview(y_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);
#endif

  Kokkos::parallel_for(
    "elementWiseDivide", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i) { z_1d(i) = x_1d(i) / y_1d(i); });
}

/// Add constant to all vector elements: z(i) = x(i) + b
inline void addConst(const vector_type& x, scalar_type b, vector_type& z)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (z.need_sync<memory_space>())
    const_cast<vector_type&>(z).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>();
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);

  z.modify<memory_space>();
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);
#endif

  Kokkos::parallel_for(
    "addConst", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i) { z_1d(i) = x_1d(i) + b; });
}

/// Compare vector elements to c: z(i) = |x(i)| >= c ? 1 : 0
inline void compare(scalar_type c, const vector_type& x, vector_type& z)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (z.need_sync<memory_space>())
    const_cast<vector_type&>(z).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>();
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);

  z.modify<memory_space>();
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);
#endif

  Kokkos::parallel_for(
    "compare", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i) {
      z_1d(i) = std::abs(x_1d(i)) >= c ? static_cast<scalar_type>(1.0)
                                       : static_cast<scalar_type>(0.0);
    });
}

/*----------------------------------------------------------------
     *  Reduction vector kernels
     *---------------------------------------------------------------*/

/// Weighted root-mean-square norm
inline mag_type normWrms(const vector_type& x, const vector_type& w)
{
  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = x.getMap()->getComm();
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());
  const global_ordinal_type Nglob =
    static_cast<global_ordinal_type>(x.getGlobalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (w.need_sync<memory_space>())
    const_cast<vector_type&>(w).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d = w.getLocalView<memory_space>();
  auto w_1d = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d = w.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto w_1d = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "normWrms", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      local_sum += x_1d(i) * w_1d(i) * (x_1d(i) * w_1d(i));
    },
    sum);

  mag_type globalSum = static_cast<scalar_type>(0.0);
  reduceAll<int, mag_type>(*comm, REDUCE_SUM, sum, outArg(globalSum));
  return std::sqrt(globalSum / static_cast<mag_type>(Nglob));
}

/// Weighted root-mean-square norm with mask
inline mag_type normWrmsMask(const vector_type& x, const vector_type& w,
                             const vector_type& id)
{
  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = x.getMap()->getComm();
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());
  const global_ordinal_type Nglob =
    static_cast<global_ordinal_type>(x.getGlobalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (w.need_sync<memory_space>())
    const_cast<vector_type&>(w).sync<memory_space>();
  if (id.need_sync<memory_space>())
    const_cast<vector_type&>(id).sync<memory_space>();

  auto x_2d  = x.getLocalView<memory_space>();
  auto x_1d  = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d  = w.getLocalView<memory_space>();
  auto w_1d  = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
  auto id_2d = id.getLocalView<memory_space>();
  auto id_1d = Kokkos::subview(id_2d, Kokkos::ALL(), 0);
#else
  auto x_2d  = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d  = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d  = w.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto w_1d  = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
  auto id_2d = id.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto id_1d = Kokkos::subview(id_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "normWrmsMask", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      if (id_1d(i) > static_cast<scalar_type>(0.0))
        local_sum += x_1d(i) * w_1d(i) * (x_1d(i) * w_1d(i));
    },
    sum);

  mag_type globalSum = static_cast<scalar_type>(0.0);
  reduceAll<int, mag_type>(*comm, REDUCE_SUM, sum, outArg(globalSum));
  return std::sqrt(globalSum / static_cast<mag_type>(Nglob));
}

/// Find minimum element value in the vector
inline scalar_type minElement(const vector_type& x)
{
  using namespace Kokkos;

  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = x.getMap()->getComm();
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#endif

  scalar_type minimum;
  Min<scalar_type> min_reducer(minimum);

  Kokkos::parallel_reduce(
    "minElement", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_min) {
      min_reducer.join(local_min, x_1d(i));
    },
    min_reducer);

  scalar_type globalMin;
  reduceAll<int, scalar_type>(*comm, REDUCE_MIN, minimum, outArg(globalMin));
  return globalMin;
}

/// Weighted L2 norm
inline mag_type normWL2(const vector_type& x, const vector_type& w)
{
  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = x.getMap()->getComm();
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (w.need_sync<memory_space>())
    const_cast<vector_type&>(w).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d = w.getLocalView<memory_space>();
  auto w_1d = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d = w.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto w_1d = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "normWL2", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      local_sum += x_1d(i) * w_1d(i) * (x_1d(i) * w_1d(i));
    },
    sum);

  mag_type globalSum = static_cast<scalar_type>(0.0);
  reduceAll<int, mag_type>(*comm, REDUCE_SUM, sum, outArg(globalSum));
  return std::sqrt(globalSum);
}

/// Elementwise inverse, return false if any denominator is zero.
inline bool invTest(const vector_type& x, vector_type& z)
{
  using namespace Kokkos;

  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = x.getMap()->getComm();
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>();
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);

  z.modify<memory_space>();
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);
#endif

  scalar_type minimum;
  Min<scalar_type> min_reducer(minimum);

  Kokkos::parallel_reduce(
    "invTest", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_min) {
      if (x_1d(i) == static_cast<scalar_type>(0.0))
      {
        min_reducer.join(local_min, static_cast<scalar_type>(0.0));
      }
      else { z_1d(i) = static_cast<scalar_type>(1.0) / x_1d(i); }
    },
    min_reducer);

  scalar_type globalMin;
  reduceAll<int, scalar_type>(*comm, REDUCE_MIN, minimum, outArg(globalMin));
  return (globalMin > static_cast<scalar_type>(0.5));
}

/// Find constraint violations
inline bool constraintMask(const vector_type& c, const vector_type& x,
                           vector_type& m)
{
  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = x.getMap()->getComm();
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (c.need_sync<memory_space>())
    const_cast<vector_type&>(c).sync<memory_space>();
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto c_2d = c.getLocalView<memory_space>();
  auto c_1d = Kokkos::subview(c_2d, Kokkos::ALL(), 0);
  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto m_2d = m.getLocalView<memory_space>();
  auto m_1d = Kokkos::subview(m_2d, Kokkos::ALL(), 0);

  m.modify<memory_space>();
#else
  auto c_2d = c.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto c_1d = Kokkos::subview(c_2d, Kokkos::ALL(), 0);
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto m_2d = m.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto m_1d = Kokkos::subview(m_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "constraintMask", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      const bool test = (std::abs(c_1d(i)) > static_cast<scalar_type>(1.5) &&
                         c_1d(i) * x_1d(i) <= static_cast<scalar_type>(0.0)) ||
                        (std::abs(c_1d(i)) > static_cast<scalar_type>(0.5) &&
                         c_1d(i) * x_1d(i) < static_cast<scalar_type>(0.0));
      m_1d(i) = test ? static_cast<scalar_type>(1.0)
                     : static_cast<scalar_type>(0.0);
      local_sum += m_1d(i);
    },
    sum);

  mag_type globalSum = static_cast<scalar_type>(0.0);
  reduceAll<int, mag_type>(*comm, REDUCE_SUM, sum, outArg(globalSum));
  return (globalSum < static_cast<scalar_type>(0.5));
}

/// Minimum quotient: min_i(num(i)/den(i))
inline scalar_type minQuotient(const vector_type& num, const vector_type& den)
{
  using namespace Kokkos;

  const Teuchos::RCP<const Teuchos::Comm<int>>& comm = num.getMap()->getComm();
  const local_ordinal_type N =
    static_cast<local_ordinal_type>(num.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (num.need_sync<memory_space>())
    const_cast<vector_type&>(num).sync<memory_space>();
  if (den.need_sync<memory_space>())
    const_cast<vector_type&>(den).sync<memory_space>();

  auto num_2d = num.getLocalView<memory_space>();
  auto num_1d = Kokkos::subview(num_2d, Kokkos::ALL(), 0);
  auto den_2d = den.getLocalView<memory_space>();
  auto den_1d = Kokkos::subview(den_2d, Kokkos::ALL(), 0);
#else
  auto num_2d = num.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto num_1d = Kokkos::subview(num_2d, Kokkos::ALL(), 0);
  auto den_2d = den.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto den_1d = Kokkos::subview(den_2d, Kokkos::ALL(), 0);
#endif

  scalar_type minimum;
  Min<scalar_type> min_reducer(minimum);

  Kokkos::parallel_reduce(
    "minQuotient", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_min) {
      if (den_1d(i) != static_cast<scalar_type>(0.0))
        min_reducer.join(local_min, num_1d(i) / den_1d(i));
    },
    min_reducer);

  scalar_type globalMin;
  reduceAll<int, scalar_type>(*comm, REDUCE_MIN, minimum, outArg(globalMin));
  return globalMin;
}

/// MPI task-local dot-product
inline scalar_type dotProdLocal(const vector_type& x, const vector_type& y)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (y.need_sync<memory_space>())
    const_cast<vector_type&>(y).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto y_2d = y.getLocalView<memory_space>();
  auto y_1d = Kokkos::subview(y_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto y_2d = y.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto y_1d = Kokkos::subview(y_2d, Kokkos::ALL(), 0);
#endif

  scalar_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "dotProdLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_sum) {
      local_sum += x_1d(i) * y_1d(i);
    },
    sum);

  return sum;
}

/// MPI task-local maximum norm of a vector
inline mag_type maxNormLocal(const vector_type& x)
{
  using namespace Kokkos;

  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#endif

  mag_type maximum;
  Max<mag_type> max_reducer(maximum);

  Kokkos::parallel_reduce(
    "maxNormLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_max) {
      max_reducer.join(local_max, std::abs(x_1d(i)));
    },
    max_reducer);

  return maximum;
}

/// MPI task-local minimum element in the vector
inline scalar_type minLocal(const vector_type& x)
{
  using namespace Kokkos;

  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#endif

  scalar_type minimum;
  Min<scalar_type> min_reducer(minimum);

  Kokkos::parallel_reduce(
    "minElement", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_min) {
      min_reducer.join(local_min, x_1d(i));
    },
    min_reducer);

  return minimum;
}

/// MPI task-local L1 norm of a vector
inline mag_type L1NormLocal(const vector_type& x)
{
  using namespace Kokkos;

  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "L1NormLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      local_sum += std::abs(x_1d(i));
    },
    sum);

  return sum;
}

/// MPI task-local weighted squared sum
inline mag_type WSqrSumLocal(const vector_type& x, const vector_type& w)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (w.need_sync<memory_space>())
    const_cast<vector_type&>(w).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d = w.getLocalView<memory_space>();
  auto w_1d = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d = w.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto w_1d = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "WSqrSumLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      local_sum += x_1d(i) * w_1d(i) * (x_1d(i) * w_1d(i));
    },
    sum);

  return sum;
}

/// MPI task-local weighted squared masked sum
inline mag_type WSqrSumMaskLocal(const vector_type& x, const vector_type& w,
                                 const vector_type& id)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();
  if (w.need_sync<memory_space>())
    const_cast<vector_type&>(w).sync<memory_space>();
  if (id.need_sync<memory_space>())
    const_cast<vector_type&>(id).sync<memory_space>();

  auto x_2d  = x.getLocalView<memory_space>();
  auto x_1d  = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d  = w.getLocalView<memory_space>();
  auto w_1d  = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
  auto id_2d = id.getLocalView<memory_space>();
  auto id_1d = Kokkos::subview(id_2d, Kokkos::ALL(), 0);
#else
  auto x_2d  = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d  = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto w_2d  = w.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto w_1d  = Kokkos::subview(w_2d, Kokkos::ALL(), 0);
  auto id_2d = id.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto id_1d = Kokkos::subview(id_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "WSqrSumMaskLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      if (id_1d(i) > static_cast<scalar_type>(0.0))
        local_sum += x_1d(i) * w_1d(i) * (x_1d(i) * w_1d(i));
    },
    sum);

  return sum;
}

/// MPI task-local elementwise inverse, return false if any denominator is zero.
inline bool invTestLocal(const vector_type& x, vector_type& z)
{
  using namespace Kokkos;

  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>();
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);

  z.modify<memory_space>();
#else
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto z_2d = z.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto z_1d = Kokkos::subview(z_2d, Kokkos::ALL(), 0);
#endif

  scalar_type minimum;
  Min<scalar_type> min_reducer(minimum);

  Kokkos::parallel_reduce(
    "invTestLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_min) {
      if (x_1d(i) == static_cast<scalar_type>(0.0))
      {
        min_reducer.join(local_min, static_cast<scalar_type>(0.0));
      }
      else { z_1d(i) = static_cast<scalar_type>(1.0) / x_1d(i); }
    },
    min_reducer);

  return (minimum > static_cast<scalar_type>(0.5));
}

/// MPI task-local constraint violation check
inline bool constraintMaskLocal(const vector_type& c, const vector_type& x,
                                vector_type& m)
{
  const local_ordinal_type N = static_cast<local_ordinal_type>(x.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (c.need_sync<memory_space>())
    const_cast<vector_type&>(c).sync<memory_space>();
  if (x.need_sync<memory_space>())
    const_cast<vector_type&>(x).sync<memory_space>();

  auto c_2d = c.getLocalView<memory_space>();
  auto c_1d = Kokkos::subview(c_2d, Kokkos::ALL(), 0);
  auto x_2d = x.getLocalView<memory_space>();
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto m_2d = m.getLocalView<memory_space>();
  auto m_1d = Kokkos::subview(m_2d, Kokkos::ALL(), 0);

  m.modify<memory_space>();
#else
  auto c_2d = c.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto c_1d = Kokkos::subview(c_2d, Kokkos::ALL(), 0);
  auto x_2d = x.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
  auto m_2d = m.getLocalView<memory_space>(Tpetra::Access::ReadWrite);
  auto m_1d = Kokkos::subview(m_2d, Kokkos::ALL(), 0);
#endif

  mag_type sum = static_cast<scalar_type>(0.0);
  Kokkos::parallel_reduce(
    "constraintMaskLocal", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, mag_type& local_sum) {
      const bool test = (std::abs(c_1d(i)) > static_cast<scalar_type>(1.5) &&
                         c_1d(i) * x_1d(i) <= static_cast<scalar_type>(0.0)) ||
                        (std::abs(c_1d(i)) > static_cast<scalar_type>(0.5) &&
                         c_1d(i) * x_1d(i) < static_cast<scalar_type>(0.0));
      m_1d(i) = test ? static_cast<scalar_type>(1.0)
                     : static_cast<scalar_type>(0.0);
      local_sum += m_1d(i);
    },
    sum);

  return (sum < static_cast<scalar_type>(0.5));
}

/// MPI task-local minimum quotient: min_i(num(i)/den(i))
inline scalar_type minQuotientLocal(const vector_type& num, const vector_type& den)
{
  using namespace Kokkos;

  const local_ordinal_type N =
    static_cast<local_ordinal_type>(num.getLocalLength());

#if TRILINOS_MAJOR_VERSION < 14
  if (num.need_sync<memory_space>())
    const_cast<vector_type&>(num).sync<memory_space>();
  if (den.need_sync<memory_space>())
    const_cast<vector_type&>(den).sync<memory_space>();

  auto num_2d = num.getLocalView<memory_space>();
  auto num_1d = Kokkos::subview(num_2d, Kokkos::ALL(), 0);
  auto den_2d = den.getLocalView<memory_space>();
  auto den_1d = Kokkos::subview(den_2d, Kokkos::ALL(), 0);
#else
  auto num_2d = num.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto num_1d = Kokkos::subview(num_2d, Kokkos::ALL(), 0);
  auto den_2d = den.getLocalView<memory_space>(Tpetra::Access::ReadOnly);
  auto den_1d = Kokkos::subview(den_2d, Kokkos::ALL(), 0);
#endif

  scalar_type minimum;
  Min<scalar_type> min_reducer(minimum);

  Kokkos::parallel_reduce(
    "minQuotient", Kokkos::RangePolicy<execution_space>(0, N),
    KOKKOS_LAMBDA(const local_ordinal_type& i, scalar_type& local_min) {
      if (den_1d(i) != static_cast<scalar_type>(0.0))
        min_reducer.join(local_min, num_1d(i) / den_1d(i));
    },
    min_reducer);

  return minimum;
}

} // namespace nvector_tpetra

} // namespace trilinos
} // namespace sundials

#endif // _TPETRA_SUNDIALS_VECTOR_KERNELS_HPP_
