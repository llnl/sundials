/* -----------------------------------------------------------------------------
 * A std::vector-based NVector using NVectorBase.
 *
 * Usage:
 *
 *   // 1. Allocating
 *   StdVector v(100, sunctx);
 *
 *   // 2. Move from existing data
 *   std::vector<sunrealtype> my_data(100, 1.0);
 *   StdVector v(std::move(my_data), sunctx);  // my_data is now empty
 *
 *   // 3. View existing data (no copy, no ownership transfer)
 *   std::vector<sunrealtype> my_data(100, 1.0);
 *   StdVector v(my_data, sunctx);  // v reads/writes my_data directly
 *   // my_data must outlive v!
 * ---------------------------------------------------------------------------*/

#ifndef NVECTOR_STDVECTOR_HPP
#define NVECTOR_STDVECTOR_HPP

#include "sundials/sundials_nvector_base.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace sundials {

class StdVector : public NVectorBase<StdVector>
{
public:
  // ---------------------------------------------------------------------------
  // Construction
  // ---------------------------------------------------------------------------

  // Allocating constructor. Creates and owns a vector.
  StdVector(sunindextype len, SUNContext sunctx)
    : NVectorBase<StdVector>(sunctx), owned_(len), ptr_(owned_.data()), len_(len)
  {}

  // Move constructor. Takes ownership of an existing std::vector. The source
  // is left in a moved-from state (typically empty).
  StdVector(std::vector<sunrealtype>&& v, SUNContext sunctx)
    : NVectorBase<StdVector>(sunctx),
      owned_(std::move(v)),
      ptr_(owned_.data()),
      len_(static_cast<sunindextype>(owned_.size()))
  {}

  // View constructor. Wraps an existing std::vector WITHOUT copying or taking
  // ownership. The caller must ensure the source outlives this StdVector.
  // Reads and writes go directly to the source vector's storage.
  StdVector(std::vector<sunrealtype>& v, SUNContext sunctx)
    : NVectorBase<StdVector>(sunctx),
      owned_(),
      ptr_(v.data()),
      len_(static_cast<sunindextype>(v.size()))
  {}

  // Clones always allocate fresh owning storage. SUNDIALS takes ownership of
  // the clone and frees it via N_VDestroy.
  StdVector* cloneNew(SUNContext sunctx) const
  { return new StdVector(len_, sunctx); }

  // ---------------------------------------------------------------------------
  // Utilities
  // ---------------------------------------------------------------------------

  sunindextype length() const { return len_; }

  // Returns a pointer to the underlying data. Always valid for both
  // owning and view modes.
  sunrealtype* data() { return ptr_; }

  const sunrealtype* data() const { return ptr_; }

  // Convenience operator for element access.
  sunrealtype& operator[](sunindextype i) { return ptr_[i]; }

  const sunrealtype& operator[](sunindextype i) const { return ptr_[i]; }

  // Returns true if this StdVector owns its data (allocated or moved-in).
  // Returns false for views.
  bool owning() const { return !owned_.empty(); }

  sunrealtype* arrayPointer() { return ptr_; }

  // ---------------------------------------------------------------------------
  // Standard vector operations
  // ---------------------------------------------------------------------------

  void linearSum(sunrealtype a, const StdVector& x, sunrealtype b,
                 const StdVector& y)
  {
    for (sunindextype i = 0; i < len_; ++i)
      ptr_[i] = a * x.ptr_[i] + b * y.ptr_[i];
  }

  void setConst(sunrealtype c) { std::fill(ptr_, ptr_ + len_, c); }

  void prod(const StdVector& x, const StdVector& y)
  {
    for (sunindextype i = 0; i < len_; ++i) ptr_[i] = x.ptr_[i] * y.ptr_[i];
  }

  void div(const StdVector& x, const StdVector& y)
  {
    for (sunindextype i = 0; i < len_; ++i) ptr_[i] = x.ptr_[i] / y.ptr_[i];
  }

  void scale(sunrealtype c, const StdVector& x)
  {
    for (sunindextype i = 0; i < len_; ++i) ptr_[i] = c * x.ptr_[i];
  }

  void abs(const StdVector& x)
  {
    for (sunindextype i = 0; i < len_; ++i) ptr_[i] = std::abs(x.ptr_[i]);
  }

  void inv(const StdVector& x)
  {
    for (sunindextype i = 0; i < len_; ++i)
      ptr_[i] = sunrealtype{1} / x.ptr_[i];
  }

  void addConst(const StdVector& x, sunrealtype b)
  {
    for (sunindextype i = 0; i < len_; ++i) ptr_[i] = x.ptr_[i] + b;
  }

  sunrealtype dotProd(const StdVector& y) const
  {
    sunrealtype s = 0;
    for (sunindextype i = 0; i < len_; ++i) s += ptr_[i] * y.ptr_[i];
    return s;
  }

  sunrealtype maxNorm() const
  {
    sunrealtype r = 0;
    for (sunindextype i = 0; i < len_; ++i) r = std::max(r, std::abs(ptr_[i]));
    return r;
  }

  sunrealtype wrmsNorm(const StdVector& w) const
  { return std::sqrt(wsqrSumLocal(w) / static_cast<sunrealtype>(len_)); }

  sunrealtype wrmsNormMask(const StdVector& w, const StdVector& id) const
  {
    return std::sqrt(wsqrSumMaskLocal(w, id) / static_cast<sunrealtype>(len_));
  }

  sunrealtype min() const { return *std::min_element(ptr_, ptr_ + len_); }

  sunrealtype wl2Norm(const StdVector& w) const
  { return std::sqrt(wsqrSumLocal(w)); }

  sunrealtype l1Norm() const
  {
    sunrealtype s = 0;
    for (sunindextype i = 0; i < len_; ++i) s += std::abs(ptr_[i]);
    return s;
  }

  void compare(sunrealtype c, const StdVector& x)
  {
    for (sunindextype i = 0; i < len_; ++i)
      ptr_[i] = std::abs(x.ptr_[i]) >= c ? 1.0 : 0.0;
  }

  sunbooleantype invTest(const StdVector& x)
  {
    for (sunindextype i = 0; i < len_; ++i)
    {
      if (x.ptr_[i] == 0.0) return SUNFALSE;
      ptr_[i] = 1.0 / x.ptr_[i];
    }
    return SUNTRUE;
  }

  sunbooleantype constrMask(const StdVector& c, const StdVector& x)
  {
    sunrealtype sum = 0;
    for (sunindextype i = 0; i < len_; ++i)
    {
      bool test = (std::abs(c.ptr_[i]) > 1.5 && c.ptr_[i] * x.ptr_[i] <= 0.0) ||
                  (std::abs(c.ptr_[i]) > 0.5 && c.ptr_[i] * x.ptr_[i] < 0.0);
      ptr_[i]   = test ? 1.0 : 0.0;
      sum += ptr_[i];
    }
    return (sum < 0.5);
  }

  sunrealtype minQuotient(const StdVector& denom) const
  {
    sunrealtype r = std::numeric_limits<sunrealtype>::max();
    for (sunindextype i = 0; i < len_; ++i)
      if (denom.ptr_[i] != 0.0) r = std::min(r, ptr_[i] / denom.ptr_[i]);
    return r;
  }

  // ---------------------------------------------------------------------------
  // Fused ops
  // ---------------------------------------------------------------------------

  SUNErrCode linearCombination(const sunrealtype* c,
                               const std::vector<const StdVector*>& X)
  {
    setConst(0.0);
    for (size_t j = 0; j < X.size(); ++j)
      for (sunindextype i = 0; i < len_; ++i) ptr_[i] += c[j] * X[j]->ptr_[i];
    return SUN_SUCCESS;
  }

  SUNErrCode scaleAddMulti(const sunrealtype* a,
                           const std::vector<const StdVector*>& Y,
                           const std::vector<StdVector*>& Z)
  {
    for (size_t j = 0; j < Y.size(); ++j)
      for (sunindextype i = 0; i < len_; ++i)
        Z[j]->ptr_[i] = a[j] * ptr_[i] + Y[j]->ptr_[i];
    return SUN_SUCCESS;
  }

  SUNErrCode dotProdMulti(const std::vector<const StdVector*>& Y,
                          sunrealtype* dotprods) const
  {
    for (size_t j = 0; j < Y.size(); ++j)
    {
      dotprods[j] = 0.0;
      for (sunindextype i = 0; i < len_; ++i)
        dotprods[j] += ptr_[i] * Y[j]->ptr_[i];
    }
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Vector array ops
  // ---------------------------------------------------------------------------

  SUNErrCode linearSumVectorArray(sunrealtype a,
                                  const std::vector<const StdVector*>& X,
                                  sunrealtype b,
                                  const std::vector<const StdVector*>& Y,
                                  const std::vector<StdVector*>& Z)
  {
    for (size_t k = 0; k < X.size(); ++k) Z[k]->linearSum(a, *X[k], b, *Y[k]);
    return SUN_SUCCESS;
  }

  SUNErrCode scaleVectorArray(const sunrealtype* c,
                              const std::vector<const StdVector*>& X,
                              const std::vector<StdVector*>& Z)
  {
    for (size_t k = 0; k < X.size(); ++k) Z[k]->scale(c[k], *X[k]);
    return SUN_SUCCESS;
  }

  SUNErrCode constVectorArray(sunrealtype c, const std::vector<StdVector*>& Z)
  {
    for (size_t k = 0; k < Z.size(); ++k) Z[k]->setConst(c);
    return SUN_SUCCESS;
  }

  SUNErrCode wrmsNormVectorArray(const std::vector<const StdVector*>& X,
                                 const std::vector<const StdVector*>& W,
                                 sunrealtype* nrm) const
  {
    for (size_t k = 0; k < X.size(); ++k) nrm[k] = X[k]->wrmsNorm(*W[k]);
    return SUN_SUCCESS;
  }

  SUNErrCode wrmsNormMaskVectorArray(const std::vector<const StdVector*>& X,
                                     const std::vector<const StdVector*>& W,
                                     const StdVector& id, sunrealtype* nrm) const
  {
    for (size_t k = 0; k < X.size(); ++k)
      nrm[k] = X[k]->wrmsNormMask(*W[k], id);
    return SUN_SUCCESS;
  }

  SUNErrCode scaleAddMultiVectorArray(
    const sunrealtype* a, const std::vector<const StdVector*>& X,
    const std::vector<std::vector<const StdVector*>>& Y,
    const std::vector<std::vector<StdVector*>>& Z)
  {
    for (size_t j = 0; j < Y.size(); ++j)
      for (size_t i = 0; i < X.size(); ++i)
        Z[j][i]->linearSum(a[j], *X[i], 1.0, *Y[j][i]);
    return SUN_SUCCESS;
  }

  SUNErrCode linearCombinationVectorArray(
    const sunrealtype* c, const std::vector<std::vector<const StdVector*>>& X,
    const std::vector<StdVector*>& Z)
  {
    for (size_t i = 0; i < Z.size(); ++i)
    {
      Z[i]->setConst(0.0);
      for (size_t j = 0; j < X.size(); ++j)
        for (sunindextype k = 0; k < Z[i]->len_; ++k)
          Z[i]->ptr_[k] += c[j] * X[j][i]->ptr_[k];
    }
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Local reduction kernels
  // ---------------------------------------------------------------------------

  sunrealtype wsqrSumLocal(const StdVector& w) const
  {
    sunrealtype s = 0;
    for (sunindextype i = 0; i < len_; ++i)
      s += ptr_[i] * w.ptr_[i] * ptr_[i] * w.ptr_[i];
    return s;
  }

  sunrealtype wsqrSumMaskLocal(const StdVector& w, const StdVector& id) const
  {
    sunrealtype s = 0;
    for (sunindextype i = 0; i < len_; ++i)
      if (id.ptr_[i] > 0.0) s += ptr_[i] * w.ptr_[i] * ptr_[i] * w.ptr_[i];
    return s;
  }

  // ---------------------------------------------------------------------------
  // Single buffer reduction
  // ---------------------------------------------------------------------------

  SUNErrCode dotProdMultiLocal(const std::vector<const StdVector*>& Y,
                               sunrealtype* dotprods) const
  {
    for (size_t j = 0; j < Y.size(); ++j)
    {
      dotprods[j] = 0.0;
      for (sunindextype i = 0; i < len_; ++i)
        dotprods[j] += ptr_[i] * Y[j]->ptr_[i];
    }
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Buffer Pack/Unpack
  // ---------------------------------------------------------------------------

  SUNErrCode bufSize(sunindextype* size) const
  {
    *size = len_ * static_cast<sunindextype>(sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  SUNErrCode bufPack(void* buf) const
  {
    std::memcpy(buf, ptr_, len_ * sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  SUNErrCode bufUnpack(void* buf)
  {
    std::memcpy(ptr_, buf, len_ * sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  // ---------------------------------------------------------------------------
  // Print
  // ---------------------------------------------------------------------------

  void print() const
  {
    for (sunindextype i = 0; i < len_; ++i)
      std::printf("  data[%ld] = %g\n", (long)i, ptr_[i]);
  }

  void printFile(FILE* f) const
  {
    for (sunindextype i = 0; i < len_; ++i)
      std::fprintf(f, "  data[%ld] = %g\n", (long)i, ptr_[i]);
  }

private:
  std::vector<sunrealtype> owned_; // Storage when owning (empty for views)
  sunrealtype* ptr_;               // Active data pointer (always valid)
  sunindextype len_;               // Number of elements
};

} // namespace sundials

#endif // NVECTOR_STDVECTOR_HPP
