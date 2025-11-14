/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel McGreer, Cody Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the header file an NVector implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#ifndef _NVECTOR_KOKKOS_HPP
#define _NVECTOR_KOKKOS_HPP

#include <Kokkos_Core.hpp>
#include <memory>
#include <sundials/sundials_nvector.hpp>

namespace sundials {
namespace kokkos {

// Forward declaration
template<class ExecutionSpace, class MemorySpace, class MemoryTraits>
class Vector;

// Get the Kokkos vector wrapped by an N_Vector
template<class VectorType>
inline VectorType* GetVec(N_Vector v)
{
  return static_cast<VectorType*>(v->content);
}

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

/* N_Vector API ops */

inline N_Vector_ID N_VGetVectorID_Kokkos(N_Vector v)
{
  return SUNDIALS_NVEC_KOKKOS;
}

template<class VectorType>
sunindextype N_VGetLength_Kokkos(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  return static_cast<sunindextype>(vec->Length());
}

template<class VectorType>
sunrealtype* N_VGetArrayPointer_Kokkos(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  return vec->HostView().data();
}

template<class VectorType>
sunrealtype* N_VGetDeviceArrayPointer_Kokkos(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  return vec->View().data();
}

template<class VectorType>
N_Vector N_VClone_Kokkos(N_Vector w)
{
  auto vec{GetVec<VectorType>(w)};
  auto new_vec = VectorType::Clone(*vec, w->sunctx);
  return new_vec->Convert();
}

template<class VectorType>
void N_VDestroy_Kokkos(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  delete vec;
  return;
}

template<class VectorType>
void N_VPrint_Kokkos(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
}

template<class VectorType>
void N_VPrintFile_Kokkos(N_Vector v, FILE* outfile)
{
  auto vec{GetVec<VectorType>(v)};
}

/* OPTIONAL local reduction kernels (no parallel communication) */

template<class VectorType>
sunrealtype N_VWSqrSumLocal_Kokkos(N_Vector x, N_Vector w)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto wvec{GetVec<VectorType>(w)};
  auto wdata{wvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
    "N_VWSqrSumLocal", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
    },
    gpu_result);

  return gpu_result;
}

template<class VectorType>
sunrealtype N_VWSqrSumMaskLocal_Kokkos(N_Vector x, N_Vector w, N_Vector id)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto wvec{GetVec<VectorType>(w)};
  auto wdata{wvec->View()};
  auto idvec{GetVec<VectorType>(id)};
  auto iddata{idvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
    "N_VWSqrSumMaskLocal", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      if (iddata(i) > sunrealtype{0.0})
        update += (xdata(i) * wdata(i) * xdata(i) * wdata(i));
    },
    gpu_result);

  return gpu_result;
}

/* standard vector operations */

template<class VectorType>
void N_VAbs_Kokkos(N_Vector x, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VAbs", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = std::abs(xdata(i)); });
}

template<class VectorType>
void N_VAddConst_Kokkos(N_Vector x, sunrealtype b, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VAddConst", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = xdata(i) + b; });
}

template<class VectorType>
void N_VCompare_Kokkos(sunrealtype c, N_Vector x, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VCompare", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) {
      zdata(i) = std::abs(xdata(i)) >= c ? sunrealtype{1.0} : sunrealtype{0.0};
    });
}

template<class VectorType>
void N_VConst_Kokkos(sunrealtype c, N_Vector z)
{
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VConst", typename VectorType::range_policy(0, zvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = c; });
}

template<class VectorType>
sunbooleantype N_VConstrMask_Kokkos(N_Vector c, N_Vector x, N_Vector m)
{
  auto cvec{GetVec<VectorType>(c)};
  auto cdata{cvec->View()};
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto mvec{GetVec<VectorType>(m)};
  auto mdata{mvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype sum{0.0};
  Kokkos::parallel_reduce(
    "N_VConstrMask", typename VectorType::range_policy(0, mvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      bool test = (std::abs(cdata(i)) > sunrealtype{1.5} &&
                   cdata(i) * xdata(i) <= sunrealtype{0.0}) ||
                  (std::abs(cdata(i)) > sunrealtype{0.5} &&
                   cdata(i) * xdata(i) < sunrealtype{0.0});
      mdata(i) = test ? sunrealtype{1.0} : sunrealtype{0.0};
      update += mdata(i);
    },
    sum);

  return (sum < sunrealtype{0.5});
}

template<class VectorType>
void N_VDiv_Kokkos(N_Vector x, N_Vector y, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto yvec{GetVec<VectorType>(y)};
  auto zvec{GetVec<VectorType>(z)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VDiv", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = xdata(i) / ydata(i); });
}

template<class VectorType>
sunrealtype N_VDotProd_Kokkos(N_Vector x, N_Vector y)
{
  auto xvec{GetVec<VectorType>(x)};
  auto yvec{GetVec<VectorType>(y)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
    "N_VDotProd", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      update += xdata(i) * ydata(i);
    },
    gpu_result);

  return gpu_result;
}

template<class VectorType>
void N_VInv_Kokkos(N_Vector x, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VInv", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = sunrealtype{1.0} / xdata(i); });
}

template<class VectorType>
sunbooleantype N_VInvTest_Kokkos(N_Vector x, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype minimum{0.0};
  Kokkos::parallel_reduce(
    "N_VInvTest", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      if (xdata(i) == sunrealtype{0.0}) { update += sunrealtype{1.0}; }
      else { zdata(i) = sunrealtype{1.0} / xdata(i); }
    },
    minimum);

  return (minimum < sunrealtype{0.5});
}

template<class VectorType>
sunrealtype N_VL1Norm_Kokkos(N_Vector x)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
    "N_VL1Norm", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      update += (std::abs(xdata(i)));
    },
    gpu_result);

  return gpu_result;
}

template<class VectorType>
void N_VLinearSum_Kokkos(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                         N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto yvec{GetVec<VectorType>(y)};
  auto zvec{GetVec<VectorType>(z)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VLinearSum", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = a * xdata(i) + b * ydata(i); });
}

template<class VectorType>
sunrealtype N_VMaxNorm_Kokkos(N_Vector x)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{0.0};
  Kokkos::parallel_reduce(
    "N_VMaxNorm", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      if (std::abs(xdata(i)) > update) update = std::abs(xdata(i));
    },
    Kokkos::Max<sunrealtype>(gpu_result));

  return gpu_result;
}

template<class VectorType>
sunrealtype N_VMin_Kokkos(N_Vector x)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{std::numeric_limits<sunrealtype>::max()};
  Kokkos::parallel_reduce(
    "N_VMin", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      if (xdata(i) < update) update = xdata(i);
    },
    Kokkos::Min<sunrealtype>(gpu_result));

  return gpu_result;
}

template<class VectorType>
sunrealtype N_VMinQuotient_Kokkos(N_Vector num, N_Vector denom)
{
  auto nvec{GetVec<VectorType>(num)};
  auto ndata{nvec->View()};
  auto dvec{GetVec<VectorType>(denom)};
  auto ddata{dvec->View()};

  using size_type = typename VectorType::size_type;

  sunrealtype gpu_result{std::numeric_limits<sunrealtype>::max()};
  Kokkos::parallel_reduce(
    "N_VMinQuotient", typename VectorType::range_policy(0, nvec->Length()),
    KOKKOS_LAMBDA(const size_type i, sunrealtype& update) {
      if (ddata(i) != sunrealtype{0.0})
      {
        if ((ndata(i) / ddata(i)) < update) update = ndata(i) / ddata(i);
      }
    },
    Kokkos::Min<sunrealtype>(gpu_result));

  return gpu_result;
}

template<class VectorType>
void N_VProd_Kokkos(N_Vector x, N_Vector y, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto yvec{GetVec<VectorType>(y)};
  auto zvec{GetVec<VectorType>(z)};
  auto xdata{xvec->View()};
  auto ydata{yvec->View()};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VProd", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = xdata(i) * ydata(i); });
}

template<class VectorType>
void N_VScale_Kokkos(sunrealtype c, N_Vector x, N_Vector z)
{
  auto xvec{GetVec<VectorType>(x)};
  auto xdata{xvec->View()};
  auto zvec{GetVec<VectorType>(z)};
  auto zdata{zvec->View()};

  using size_type = typename VectorType::size_type;

  Kokkos::parallel_for(
    "N_VScale", typename VectorType::range_policy(0, xvec->Length()),
    KOKKOS_LAMBDA(const size_type i) { zdata(i) = c * xdata(i); });
}

template<class VectorType>
sunrealtype N_VWL2Norm_Kokkos(N_Vector x, N_Vector w)
{
  return std::sqrt(impl::N_VWSqrSumLocal_Kokkos<VectorType>(x, w));
}

template<class VectorType>
sunrealtype N_VWrmsNorm_Kokkos(N_Vector x, N_Vector w)
{
  auto xvec{GetVec<VectorType>(x)};
  return std::sqrt(impl::N_VWSqrSumLocal_Kokkos<VectorType>(x, w) /
                   static_cast<sunrealtype>(xvec->Length()));
}

template<class VectorType>
sunrealtype N_VWrmsNormMask_Kokkos(N_Vector x, N_Vector w, N_Vector id)
{
  auto xvec{GetVec<VectorType>(x)};
  return std::sqrt(impl::N_VWSqrSumMaskLocal_Kokkos<VectorType>(x, w, id) /
                   static_cast<sunrealtype>(xvec->Length()));
}

} // namespace impl

// =============================================================================
// Public namespace
// =============================================================================

template<class ExecutionSpace = Kokkos::DefaultExecutionSpace,
         class MemorySpace    = typename ExecutionSpace::memory_space,
         class MemoryTraits   = Kokkos::MemoryManaged>
class Vector : public sundials::impl::BaseNVector,
               public sundials::ConvertibleTo<N_Vector>
{
public:
  using view_type      = Kokkos::View<sunrealtype*, MemorySpace, MemoryTraits>;
  using host_view_type = typename view_type::HostMirror;
  using memory_space   = MemorySpace;
  using memory_traits  = MemoryTraits;
  using exec_space     = typename MemorySpace::execution_space;
  using range_policy   = Kokkos::RangePolicy<exec_space>;
  using size_type      = typename view_type::size_type;

  // Default constructor
  Vector() = default;

  Vector(size_type length, SUNContext sunctx)
    : sundials::impl::BaseNVector(sunctx),
      view_("Vector device view", length),
      host_view_(Kokkos::create_mirror_view(view_))

  {
    initNvector();
  }

  Vector(view_type view, SUNContext sunctx)
    : sundials::impl::BaseNVector(sunctx),
      view_(view),
      host_view_(Kokkos::create_mirror_view(view_))

  {
    initNvector();
  }

  Vector(view_type view, host_view_type host_view, SUNContext sunctx)
    : sundials::impl::BaseNVector(sunctx), view_(view), host_view_(host_view)
  {
    initNvector();
  }

  // Move constructor
  Vector(Vector&& that_vector) noexcept
    : sundials::impl::BaseNVector(std::move(that_vector)),
      view_(std::move(that_vector.view_)),
      host_view_(std::move(that_vector.host_view_))

  {
    initNvector();
  }

  // Copy constructor
  Vector(const Vector& that_vector)
    : sundials::impl::BaseNVector(that_vector),
      view_("Vector device view", that_vector.Length()),
      host_view_(Kokkos::create_mirror_view(view_))
  {
    initNvector();
  }

  // Move assignment
  Vector& operator=(Vector&& rhs) noexcept
  {
    view_      = std::move(rhs.view_);
    host_view_ = std::move(rhs.host_view_);

    sundials::impl::BaseNVector::operator=(std::move(rhs));

    return *this;
  }

  // Copy assignment
  Vector& operator=(const Vector& rhs)
  {
    view_      = Kokkos::View<view_type>("Vector device view", rhs.Length());
    host_view_ = Kokkos::create_mirror_view(view_);

    sundials::impl::BaseNVector::operator=(rhs);

    return *this;
  }

  // Default destructor
  virtual ~Vector() = default;

  // Accessors

  size_type Length() const { return static_cast<size_type>(view_.extent(0)); }

  view_type View() { return view_; }

  host_view_type HostView() { return host_view_; }

  // Override ConvertibleTo operations

  operator N_Vector() override { return object_.get(); }

  operator N_Vector() const override { return object_.get(); }

  N_Vector Convert() override { return object_.get(); }

  N_Vector Convert() const override { return object_.get(); }

  // Static routines to create clones of the vector that are always managed

  static Vector<exec_space, memory_space, Kokkos::MemoryManaged>* Clone(
    const Vector<exec_space, memory_space, Kokkos::MemoryManaged>& that_vector,
    SUNContext sunctx)
  {
    return new Vector<exec_space, memory_space,
                      Kokkos::MemoryManaged>(that_vector.Length(), sunctx);
  }

  static Vector<exec_space, memory_space, Kokkos::MemoryManaged>* Clone(
    const Vector<exec_space, memory_space, Kokkos::MemoryUnmanaged>& that_vector,
    SUNContext sunctx)
  {
    return new Vector<exec_space, memory_space,
                      Kokkos::MemoryManaged>(that_vector.Length(), sunctx);
  }

private:
  view_type view_;
  host_view_type host_view_;

  void initNvector()
  {
    using this_type = Vector<ExecutionSpace, MemorySpace, MemoryTraits>;

    this->object_->content = this;

    /* constructors, destructors, and utility operations */
    this->object_->ops->nvclone   = impl::N_VClone_Kokkos<this_type>;
    this->object_->ops->nvdestroy = impl::N_VDestroy_Kokkos<this_type>;
    this->object_->ops->nvgetarraypointer =
      impl::N_VGetArrayPointer_Kokkos<this_type>;
    this->object_->ops->nvgetdevicearraypointer =
      impl::N_VGetDeviceArrayPointer_Kokkos<this_type>;
    this->object_->ops->nvgetlength   = impl::N_VGetLength_Kokkos<this_type>;
    this->object_->ops->nvgetvectorid = impl::N_VGetVectorID_Kokkos;

    /* standard vector operations */
    this->object_->ops->nvabs         = impl::N_VAbs_Kokkos<this_type>;
    this->object_->ops->nvaddconst    = impl::N_VAddConst_Kokkos<this_type>;
    this->object_->ops->nvcompare     = impl::N_VCompare_Kokkos<this_type>;
    this->object_->ops->nvconst       = impl::N_VConst_Kokkos<this_type>;
    this->object_->ops->nvconstrmask  = impl::N_VConstrMask_Kokkos<this_type>;
    this->object_->ops->nvdiv         = impl::N_VDiv_Kokkos<this_type>;
    this->object_->ops->nvdotprod     = impl::N_VDotProd_Kokkos<this_type>;
    this->object_->ops->nvinv         = impl::N_VInv_Kokkos<this_type>;
    this->object_->ops->nvinvtest     = impl::N_VInvTest_Kokkos<this_type>;
    this->object_->ops->nvl1norm      = impl::N_VL1Norm_Kokkos<this_type>;
    this->object_->ops->nvlinearsum   = impl::N_VLinearSum_Kokkos<this_type>;
    this->object_->ops->nvmaxnorm     = impl::N_VMaxNorm_Kokkos<this_type>;
    this->object_->ops->nvmin         = impl::N_VMin_Kokkos<this_type>;
    this->object_->ops->nvminquotient = impl::N_VMinQuotient_Kokkos<this_type>;
    this->object_->ops->nvprod        = impl::N_VProd_Kokkos<this_type>;
    this->object_->ops->nvscale       = impl::N_VScale_Kokkos<this_type>;
    this->object_->ops->nvwl2norm     = impl::N_VWL2Norm_Kokkos<this_type>;
    this->object_->ops->nvwrmsnorm    = impl::N_VWrmsNorm_Kokkos<this_type>;
    this->object_->ops->nvwrmsnormmask = impl::N_VWrmsNormMask_Kokkos<this_type>;

    /* local reduction operations */
    this->object_->ops->nvconstrmasklocal = impl::N_VConstrMask_Kokkos<this_type>;
    this->object_->ops->nvdotprodlocal = impl::N_VDotProd_Kokkos<this_type>;
    this->object_->ops->nvinvtestlocal = impl::N_VInvTest_Kokkos<this_type>;
    this->object_->ops->nvl1normlocal  = impl::N_VL1Norm_Kokkos<this_type>;
    this->object_->ops->nvmaxnormlocal = impl::N_VMaxNorm_Kokkos<this_type>;
    this->object_->ops->nvminlocal     = impl::N_VMin_Kokkos<this_type>;
    this->object_->ops->nvminquotientlocal =
      impl::N_VMinQuotient_Kokkos<this_type>;
    this->object_->ops->nvwsqrsumlocal = impl::N_VWSqrSumLocal_Kokkos<this_type>;
    this->object_->ops->nvwsqrsummasklocal =
      impl::N_VWSqrSumMaskLocal_Kokkos<this_type>;
  }
};

template<class VectorType>
void CopyToDevice(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  CopyToDevice(*vec);
}

template<class VectorType>
void CopyFromDevice(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  CopyFromDevice(*vec);
}

template<class VectorType>
void CopyToDevice(VectorType& v)
{
  Kokkos::deep_copy(v.View(), v.HostView());
}

template<class VectorType>
void CopyFromDevice(VectorType& v)
{
  Kokkos::deep_copy(v.HostView(), v.View());
}

template<class VectorType, class view_type>
view_type GetView(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  return vec->View();
}

template<class VectorType, class host_view_type>
host_view_type GetHostView(N_Vector v)
{
  auto vec{GetVec<VectorType>(v)};
  return vec->HostView();
}

} // namespace kokkos
} // namespace sundials

#endif
