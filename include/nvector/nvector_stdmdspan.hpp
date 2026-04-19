/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * This NVector wraps an arbitrary std::mdspan and supports any rank, extents,
 * layout, and accessor.
 *
 * ---------------------------------------------------------------------------
 * C++ TECHNIQUES USED IN THIS FILE
 * ---------------------------------------------------------------------------
 *
 * PARAMETER PACKS AND INDEX SEQUENCES (std::index_sequence)
 *
 *   Problem: we need to call span_[i0, i1, ..., iN-1] but the number of indices
 *   depends on the rank. Solution: std::make_index_sequence<N> generates the
 *   compile-time list {0, 1, ..., N-1}, a template captures it as a "parameter
 *   pack" (Is...), and idx[Is]... expands to idx[0], idx[1], ..., idx[N-1].
 *   Used in MultiFor::unpack and computeSize.
 *
 * FOLD EXPRESSIONS
 *
 *   (s.extent(Is) * ...) multiplies all expanded values together: s.extent(0) *
 *   s.extent(1) * ... * s.extent(N-1). Used to compute the total element count
 *   from the extents.
 *
 * IF CONSTEXPR
 *
 *   Evaluates at compile time. Only the matching branch is compiled. Used in
 *   forEach to select the right iteration strategy based on the layout type,
 *   and in the conditionally-available arrayPointer/operator[].
 *
 * SFINAE (std::enable_if_t)
 *
 *   Makes arrayPointer() and operator[] disappear at compile time for
 *   non-contiguous layouts. See the NVectorBase header for a detailed
 *   explanation.
 *
 * DECLTYPE WITH FUNCTION DECLARATIONS (no body)
 *
 *   The RightOrder and LeftOrder type aliases use a trick: a function is
 *   DECLARED (not defined) and decltype extracts the return type without ever
 *   calling it. This is a pure compile-time computation to build the
 *   MultiFor<DimOrder...> type with the right dimension ordering.
 *
 * ---------------------------------------------------------------------------
 * LAYOUT-AWARE ITERATION
 * ---------------------------------------------------------------------------
 *
 * All element-wise operations iterate over multi-indices using nested loops.
 * The NESTING ORDER of those loops is chosen to give stride-1 access in the
 * innermost loop, which is critical for cache performance:
 *
 *   layout_right (row-major):
 *     Outermost -> dim 0, innermost -> dim N-1.
 *     Last index varies fastest, matching the memory layout.
 *     Loop order determined at compile time.
 *
 *   layout_left (column-major):
 *     Outermost -> dim N-1, innermost -> dim 0.
 *     First index varies fastest, matching the memory layout.
 *     Loop order determined at compile time.
 *
 *   layout_stride (arbitrary strides):
 *     Dimensions sorted by increasing stride, so the dimension with the
 *     smallest stride (most contiguous) is the innermost loop.
 *     Order computed at construction from the mapping's runtime strides
 *     and cached for the lifetime of the object.
 *
 *   Unknown/custom layouts:
 *     Default to layout_right order (dim 0 outermost).
 *
 * The iteration order does NOT affect correctness - all N_Vector operations are
 * element-wise and order-independent. Only cache performance differs.
 *
 * In all cases, the lambda/callback receives indices in DIMENSION order (i0,
 * i1, ..., iN-1) regardless of which loop is outermost. This means
 * span_[is...] always accesses the correct element through the mapping.
 *
 * ---------------------------------------------------------------------------
 * CONTIGUOUS-ONLY FLAT ACCESS
 * ---------------------------------------------------------------------------
 *
 * For layouts that are both unique and exhaustive (layout_right, layout_left,
 * and equivalent custom layouts), arrayPointer() and flat operator[] are
 * conditionally available via SFINAE. For non-contiguous layouts these methods
 * do not exist and the corresponding SUNDIALS ops pointers stay NULL.
 *
 * ---------------------------------------------------------------------------
 * MDSPAN REQUIREMENTS
 * ---------------------------------------------------------------------------
 *
 *   MdSpanType::value_type               -> must be sunrealtype
 *   MdSpanType::rank()                   -> constexpr, number of dimensions
 *   MdSpanType::layout_type              -> the layout policy
 *   MdSpanType::mapping_type             -> the layout mapping type
 *   MdSpanType::accessor_type            -> the accessor type
 *   span.data_handle()                   -> pointer to storage
 *   span.extent(r)                       -> size of dimension r
 *   span.extents()                       -> the extents object
 *   span.mapping()                       -> the mapping object
 *   span.accessor()                      -> the accessor object
 *   span[i, j, ...]                      -> element access (C++23)
 *   MdSpanType(ptr, mapping, accessor)   -> three-argument constructor
 *   mapping.required_span_size()         -> total memory span needed
 *   mapping.stride(r)                    -> stride of dimension r
 *                                          (only for layout_stride dispatch)
 * ---------------------------------------------------------------------------*/

#ifndef NVECTOR_MDSPAN_HPP
#define NVECTOR_MDSPAN_HPP

#include "sundials/sundials_nvector_base.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>

namespace sundials {
namespace example {

// ============================================================================
// Iteration helpers
//
// The core problem: we need nested loops over N dimensions where N is a
// compile-time constant, and the loop nesting order varies by layout. These
// helpers generate the loops at compile time and call a user-provided function
// f(i0, i1, ..., iN-1) for every multi-index.
//
// The indices are always passed to f in DIMENSION order (i0 = dim 0, i1 = dim
// 1, etc.) regardless of which dimension is the outermost loop. This is
// achieved by storing indices in an array at their correct position (idx[dim] =
// i) and then unpacking the array into the function call using the
// index_sequence technique.
// ============================================================================

namespace impl {

// --------------------------------------------------------------------------
// Compile-time ordered iteration (MultiFor)
//
// DimOrder... specifies the loop nesting from outermost to innermost.
// Example: MultiFor<2, 1, 0> generates:
//   for (idx[2] = 0..ext(2))     <- outermost
//     for (idx[1] = 0..ext(1))
//       for (idx[0] = 0..ext(0)) <- innermost
//         f(idx[0], idx[1], idx[2])  <- always dimension order
//
// The nest() method is recursive: each call handles one loop level and recurses
// for the next. "Level" is a compile-time counter that advances through the
// DimOrder list. When Level reaches Rank, the recursion ends and unpack()
// calls f with the accumulated indices.
// --------------------------------------------------------------------------

template<size_t... DimOrder>
struct MultiFor
{
  static constexpr size_t Rank = sizeof...(DimOrder);

  template<class Extents, class Func>
  static void apply(const Extents& ext, Func&& f)
  {
    std::array<size_t, Rank> idx{};
    nest<0>(ext, f, idx);
  }

private:
  // The dimension ordering stored as a compile-time array.
  // order[0] = outermost dimension, order[Rank-1] = innermost.
  static constexpr size_t order[] = {DimOrder...};

  // Recursive loop generator. At each level, loops over the dimension specified
  // by order[Level] and stores the index at idx[dim].
  template<size_t Level, class Extents, class Func>
  static void nest(const Extents& ext, Func& f, std::array<size_t, Rank>& idx)
  {
    constexpr size_t dim = order[Level]; // which dimension this level iterates
    for (idx[dim] = 0; idx[dim] < ext.extent(dim); ++idx[dim])
    {
      if constexpr (Level + 1 < Rank)
        nest<Level + 1>(ext, f, idx); // recurse to next loop level
      else
        unpack(f, idx, std::make_index_sequence<Rank>{}); // innermost: call f
    }
  }

  // Convert the index array into function arguments.
  // std::make_index_sequence<Rank>{} generates {0, 1, ..., Rank-1}.
  // The template captures this as Is..., then idx[Is]... expands to
  // idx[0], idx[1], ..., idx[Rank-1] - passing the indices to f in
  // dimension order regardless of the loop nesting.
  template<class Func, size_t... Is>
  static void unpack(Func& f, const std::array<size_t, Rank>& idx,
                     std::index_sequence<Is...>)
  { f(idx[Is]...); }
};

// --------------------------------------------------------------------------
// Type aliases for standard dimension orderings
//
// A function is DECLARED (not defined) whose return type is the desired
// MultiFor<...> specialization. decltype extracts the return type without
// calling the function. This is a pure compile-time mechanism to transform an
// index_sequence into a MultiFor template parameter list.
//
// RightOrder<3> = MultiFor<0, 1, 2> (dim 0 outermost -> stride-1 for row-major)
// LeftOrder<3>  = MultiFor<2, 1, 0> (dim 2 outermost -> stride-1 for col-major)
// --------------------------------------------------------------------------

// layout_right order: 0, 1, ..., Rank-1
// make_index_sequence<3> gives index_sequence<0, 1, 2>,
// the function declaration captures Is = {0, 1, 2},
// the return type is MultiFor<0, 1, 2>.
template<size_t... Is>
MultiFor<Is...> makeRightOrder(std::index_sequence<Is...>);

template<size_t Rank>
using RightOrder = decltype(makeRightOrder(std::make_index_sequence<Rank>{}));

// layout_left order: Rank-1, ..., 1, 0
// (N - 1 - Is)... reverses the sequence: for N=3 and Is={0,1,2},
// this produces {2, 1, 0}.
template<size_t N, size_t... Is>
MultiFor<(N - 1 - Is)...> makeLeftOrder(std::index_sequence<Is...>);

template<size_t Rank>
using LeftOrder = decltype(makeLeftOrder<Rank>(std::make_index_sequence<Rank>{}));

// --------------------------------------------------------------------------
// Runtime ordered iteration (for layout_stride and similar)
//
// Same semantics as the compile-time version, but the dimension order is stored
// in a runtime array. The recursive call uses a runtime level counter rather
// than a template parameter, so the compiler cannot fully unroll.
// --------------------------------------------------------------------------

template<size_t Rank>
struct MultiForRuntime
{
  template<class Extents, class Func>
  static void apply(const Extents& ext, Func&& f,
                    const std::array<size_t, Rank>& dim_order)
  {
    std::array<size_t, Rank> idx{};
    nest(ext, f, idx, dim_order, 0);
  }

private:
  template<class Extents, class Func>
  static void nest(const Extents& ext, Func& f, std::array<size_t, Rank>& idx,
                   const std::array<size_t, Rank>& dim_order, size_t level)
  {
    const size_t dim = dim_order[level];
    for (idx[dim] = 0; idx[dim] < ext.extent(dim); ++idx[dim])
    {
      if (level + 1 < Rank) nest(ext, f, idx, dim_order, level + 1);
      else unpack(f, idx, std::make_index_sequence<Rank>{});
    }
  }

  template<class Func, size_t... Is>
  static void unpack(Func& f, const std::array<size_t, Rank>& idx,
                     std::index_sequence<Is...>)
  { f(idx[Is]...); }
};

// --------------------------------------------------------------------------
// Compute dimension order sorted by increasing stride (smallest stride
// -> innermost loop -> best cache behavior).
// --------------------------------------------------------------------------

template<size_t Rank, class MappingType>
std::array<size_t, Rank> strideOrder(const MappingType& m)
{
  std::array<size_t, Rank> order;
  std::iota(order.begin(), order.end(), size_t{0});
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b) { return m.stride(a) > m.stride(b); });
  return order;
  // order[0] = largest stride (outermost), order[Rank-1] = smallest (innermost)
}

// --------------------------------------------------------------------------
// Conditional storage for the cached iteration order. Empty for layouts where
// the order is known at compile time; an array for layout_stride where it
// depends on runtime strides.
// --------------------------------------------------------------------------

struct NoOrder
{};

template<class LayoutType, size_t Rank>
struct OrderStorage
{
  using type = NoOrder;
};

} // namespace impl

// ============================================================================
// MdSpanVector
// ============================================================================

template<class MdSpanType>
class MdSpanVector : public NVectorBase<MdSpanVector<MdSpanType>>
{
  using Base = NVectorBase<MdSpanVector<MdSpanType>>;

public:
  using mdspan_type   = MdSpanType;
  using mapping_type  = typename mdspan_type::mapping_type;
  using accessor_type = typename mdspan_type::accessor_type;
  using layout_type   = typename mdspan_type::layout_type;

  static_assert(
    std::is_same_v<typename mdspan_type::value_type, sunrealtype>,
    "MdSpanVector requires an mdspan whose value_type is sunrealtype");

  static constexpr size_t Rank = mdspan_type::rank();

  // --------------------------------------------------------------------------
  // Construction
  // --------------------------------------------------------------------------

  /// View constructor. Wraps the user's mdspan. The caller must keep the
  /// underlying storage alive for the lifetime of this object.
  MdSpanVector(mdspan_type span, SUNContext sunctx)
    : Base(sunctx),
      span_(span),
      len_(computeSize(span)),
      owned_(),
      stride_order_(computeStrideOrder(span))
  {}

  // --------------------------------------------------------------------------
  // Lifecycle
  // --------------------------------------------------------------------------

  /// Clones allocate fresh storage and construct a new mdspan over it.
  ///
  /// required_span_size() is used instead of the product of extents because
  /// non-contiguous layouts (e.g., layout_stride) may need more memory than the
  /// logical element count - there can be gaps between elements.
  ///
  /// The mapping and accessor are copied from the original so the clone has the
  /// same shape, layout, and element access behavior.
  ///
  /// std::make_unique<sunrealtype[]>(n) allocates n elements on the heap and
  /// wraps the pointer in a unique_ptr that automatically frees the memory when
  /// the MdSpanVector is destroyed.
  ///
  /// std::move(storage) transfers ownership of the allocation into the new
  /// MdSpanVector without copying - the source unique_ptr becomes null.
  MdSpanVector* cloneNew(SUNContext sunctx) const
  {
    auto storage =
      std::make_unique<sunrealtype[]>(span_.mapping().required_span_size());
    mdspan_type new_span(storage.get(), span_.mapping(), span_.accessor());
    return new MdSpanVector(std::move(storage), new_span, sunctx);
  }

  /// Returns the number of logical elements (product of all extents).
  sunindextype length() const { return len_; }

  // --------------------------------------------------------------------------
  // Data access
  // --------------------------------------------------------------------------

  /// Access the underlying mdspan for multi-dimensional element access.
  mdspan_type& span() { return span_; }

  const mdspan_type& span() const { return span_; }

  /// Returns true if this object owns its backing storage (clone).
  bool owning() const { return owned_ != nullptr; }

  // -- Contiguous-only access (conditionally available) ---------------------
  //
  // These methods only exist when the layout is contiguous (both unique and
  // exhaustive). The mechanism:
  //
  //   template<class M = mapping_type,
  //            std::enable_if_t<M::is_always_exhaustive() &&
  //                             M::is_always_unique(), int> = 0>
  //   sunrealtype* arrayPointer() { ... }
  //
  // This is a template function with a default parameter M = mapping_type. The
  // enable_if_t check evaluates M::is_always_exhaustive() and
  // M::is_always_unique() at compile time:
  //
  //   - Both true (layout_right, layout_left): enable_if_t produces int,
  //     the function exists, NVectorBase detects it and installs the op.
  //
  //   - Either false (layout_stride): enable_if_t fails (SFINAE), the
  //     function silently disappears, the ops pointer stays NULL.
  //
  // The "template<class M = mapping_type>" wrapper is necessary because
  // enable_if_t must depend on a template parameter of the FUNCTION (not just
  // the class) for SFINAE to apply. The default M = mapping_type means callers
  // don't need to specify it: just call arrayPointer().

  /// Flat pointer to the data. Only available when the mapping is both unique
  /// and exhaustive (layout_right, layout_left, or equivalent).
  template<class M = mapping_type,
           std::enable_if_t<M::is_always_exhaustive() && M::is_always_unique(), int> = 0>
  sunrealtype* arrayPointer()
  { return span_.data_handle(); }

  /// Flat element access by linear index. Same contiguity requirement.
  template<class M = mapping_type,
           std::enable_if_t<M::is_always_exhaustive() && M::is_always_unique(), int> = 0>
  sunrealtype& operator[](sunindextype i)
  { return span_.data_handle()[i]; }

  template<class M = mapping_type,
           std::enable_if_t<M::is_always_exhaustive() && M::is_always_unique(), int> = 0>
  const sunrealtype& operator[](sunindextype i) const
  { return span_.data_handle()[i]; }

  // --------------------------------------------------------------------------
  // Standard vector operations
  //
  // Each operation calls forEach() with a lambda that receives the
  // multi-indices as a parameter pack (auto... is). The lambda body uses
  // span_[is...] which expands to span_[i0, i1, ..., iN-1] - the correct mdspan
  // element access for any rank.
  //
  // [&] captures all local variables by reference. This gives the lambda
  // access to span_ (for writing), x.span_ (for reading), and any scalar
  // parameters (a, b, c).
  //
  // Example expansion for a rank-2 mdspan:
  //   forEach calls the lambda with (i0, i1)
  //   span_[is...] becomes span_[i0, i1]
  //   The mdspan's operator[] maps (i0, i1) to the correct memory offset
  // --------------------------------------------------------------------------

  void linearSum(sunrealtype a, const MdSpanVector& x, sunrealtype b,
                 const MdSpanVector& y)
  {
    forEach([&](auto... is)
            { span_[is...] = a * x.span_[is...] + b * y.span_[is...]; });
  }

  void setConst(sunrealtype c)
  {
    forEach([&](auto... is) { span_[is...] = c; });
  }

  void prod(const MdSpanVector& x, const MdSpanVector& y)
  {
    forEach([&](auto... is) { span_[is...] = x.span_[is...] * y.span_[is...]; });
  }

  void div(const MdSpanVector& x, const MdSpanVector& y)
  {
    forEach([&](auto... is) { span_[is...] = x.span_[is...] / y.span_[is...]; });
  }

  void scale(sunrealtype c, const MdSpanVector& x)
  {
    forEach([&](auto... is) { span_[is...] = c * x.span_[is...]; });
  }

  void abs(const MdSpanVector& x)
  {
    forEach([&](auto... is) { span_[is...] = std::abs(x.span_[is...]); });
  }

  void inv(const MdSpanVector& x)
  {
    forEach([&](auto... is) { span_[is...] = sunrealtype{1} / x.span_[is...]; });
  }

  void addConst(const MdSpanVector& x, sunrealtype b)
  {
    forEach([&](auto... is) { span_[is...] = x.span_[is...] + b; });
  }

  sunrealtype dotProd(const MdSpanVector& y) const
  {
    sunrealtype s = 0;
    forEach([&](auto... is) { s += span_[is...] * y.span_[is...]; });
    return s;
  }

  sunrealtype maxNorm() const
  {
    sunrealtype r = 0;
    forEach([&](auto... is) { r = std::max(r, std::abs(span_[is...])); });
    return r;
  }

  sunrealtype wrmsNorm(const MdSpanVector& w) const
  { return std::sqrt(wsqrSumLocal(w) / static_cast<sunrealtype>(len_)); }

  sunrealtype wrmsNormMask(const MdSpanVector& w, const MdSpanVector& id) const
  {
    return std::sqrt(wsqrSumMaskLocal(w, id) / static_cast<sunrealtype>(len_));
  }

  sunrealtype min() const
  {
    sunrealtype r = std::numeric_limits<sunrealtype>::max();
    forEach([&](auto... is) { r = std::min(r, span_[is...]); });
    return r;
  }

  sunrealtype wl2Norm(const MdSpanVector& w) const
  { return std::sqrt(wsqrSumLocal(w)); }

  sunrealtype l1Norm() const
  {
    sunrealtype s = 0;
    forEach([&](auto... is) { s += std::abs(span_[is...]); });
    return s;
  }

  void compare(sunrealtype c, const MdSpanVector& x)
  {
    forEach([&](auto... is)
            { span_[is...] = std::abs(x.span_[is...]) >= c ? 1.0 : 0.0; });
  }

  sunbooleantype invTest(const MdSpanVector& x)
  {
    sunbooleantype result = SUNTRUE;
    forEach(
      [&](auto... is)
      {
        if (result == SUNFALSE) return;
        if (x.span_[is...] == 0.0)
        {
          result = SUNFALSE;
          return;
        }
        span_[is...] = 1.0 / x.span_[is...];
      });
    return result;
  }

  sunbooleantype constrMask(const MdSpanVector& c, const MdSpanVector& x)
  {
    sunrealtype sum = 0;
    forEach(
      [&](auto... is)
      {
        bool test    = (std::abs(c.span_[is...]) > 1.5 &&
                        c.span_[is...] * x.span_[is...] <= 0.0) ||
                       (std::abs(c.span_[is...]) > 0.5 &&
                        c.span_[is...] * x.span_[is...] < 0.0);
        span_[is...] = test ? 1.0 : 0.0;
        sum += span_[is...];
      });
    return (sum < 0.5);
  }

  sunrealtype minQuotient(const MdSpanVector& denom) const
  {
    sunrealtype r = std::numeric_limits<sunrealtype>::max();
    forEach(
      [&](auto... is)
      {
        if (denom.span_[is...] != 0.0)
          r = std::min(r, span_[is...] / denom.span_[is...]);
      });
    return r;
  }

  // --------------------------------------------------------------------------
  // Local reduction kernels
  // --------------------------------------------------------------------------

  sunrealtype wsqrSumLocal(const MdSpanVector& w) const
  {
    sunrealtype s = 0;
    forEach(
      [&](auto... is)
      { s += span_[is...] * w.span_[is...] * span_[is...] * w.span_[is...]; });
    return s;
  }

  sunrealtype wsqrSumMaskLocal(const MdSpanVector& w, const MdSpanVector& id) const
  {
    sunrealtype s = 0;
    forEach(
      [&](auto... is)
      {
        if (id.span_[is...] > 0.0)
          s += span_[is...] * w.span_[is...] * span_[is...] * w.span_[is...];
      });
    return s;
  }

  // --------------------------------------------------------------------------
  // XBraid interface
  //
  // Pack/unpack serialize logical elements into a dense buffer. The buffer
  // size is the number of logical elements, NOT required_span_size(). The
  // iteration order ensures a deterministic packing regardless of layout.
  // --------------------------------------------------------------------------

  SUNErrCode bufSize(sunindextype* size) const
  {
    *size = len_ * static_cast<sunindextype>(sizeof(sunrealtype));
    return SUN_SUCCESS;
  }

  SUNErrCode bufPack(void* buf) const
  {
    sunrealtype* b = static_cast<sunrealtype*>(buf);
    sunindextype k = 0;
    // The lambda captures k by reference and increments it for each element
    // visited. forEach iterates in a deterministic order (determined by the
    // layout), so pack/unpack visit elements in the same sequence - the buffer
    // is always consistent.
    forEach([&](auto... is) { b[k++] = span_[is...]; });
    return SUN_SUCCESS;
  }

  SUNErrCode bufUnpack(void* buf)
  {
    const sunrealtype* b = static_cast<const sunrealtype*>(buf);
    sunindextype k       = 0;
    forEach([&](auto... is) { span_[is...] = b[k++]; });
    return SUN_SUCCESS;
  }

  // --------------------------------------------------------------------------
  // Debug
  // --------------------------------------------------------------------------

  void print() const
  {
    sunindextype k = 0;
    forEach([&](auto... is)
            { std::printf("  data[%ld] = %g\n", (long)k++, span_[is...]); });
  }

  void printFile(FILE* f) const
  {
    sunindextype k = 0;
    forEach([&](auto... is)
            { std::fprintf(f, "  data[%ld] = %g\n", (long)k++, span_[is...]); });
  }

private:
  mdspan_type span_; ///< The mdspan (view or over owned storage).
  sunindextype len_; ///< Product of all extents (logical element count).
  std::unique_ptr<sunrealtype[]> owned_; ///< Backing storage for clones (null for views).

  /// Cached dimension iteration order for layout_stride. For layout_right and
  /// layout_left the iteration order is known at compile time, so this array is
  /// unused (but must exist because the member is always declared). For
  /// layout_stride, dimensions are sorted by decreasing stride at construction
  /// time so the innermost loop has the smallest stride.
  std::array<size_t, Rank> stride_order_;

  // --------------------------------------------------------------------------
  // Private clone constructor
  // --------------------------------------------------------------------------

  MdSpanVector(std::unique_ptr<sunrealtype[]> storage, mdspan_type span,
               SUNContext sunctx)
    : Base(sunctx),
      span_(span),
      len_(computeSize(span)),
      owned_(std::move(storage)),
      stride_order_(computeStrideOrder(span))
  {}

  // --------------------------------------------------------------------------
  // forEach - layout-aware dispatch
  //
  // Selects the iteration strategy based on the layout type. The decision is
  // made at compile time via if constexpr - only the matching branch is
  // compiled, so there is no runtime overhead from the other branches.
  //
  // std::forward<Func>(f) is "perfect forwarding" - it passes the lambda to the
  // iteration helper with the same value category (lvalue or rvalue) it arrived
  // with, avoiding unnecessary copies.
  //
  // forEach is marked const because it only reads span_.extents(). The lambda
  // itself determines whether span_ elements are read or written - it captures
  // span_ from the calling method's scope, where the const qualification of
  // "this" controls mutability.
  // --------------------------------------------------------------------------

  template<class Func>
  void forEach(Func&& f) const
  {
    if constexpr (Rank == 0) { f(); }
    else if constexpr (std::is_same_v<layout_type, std::layout_left>)
    {
      // Column-major: dim 0 innermost (smallest stride)
      impl::LeftOrder<Rank>::apply(span_.extents(), std::forward<Func>(f));
    }
    else if constexpr (std::is_same_v<layout_type, std::layout_stride>)
    {
      // Arbitrary strides: use cached stride-sorted order
      impl::MultiForRuntime<Rank>::apply(span_.extents(),
                                           std::forward<Func>(f), stride_order_);
    }
    else
    {
      // layout_right and unknown layouts: dim N-1 innermost
      impl::RightOrder<Rank>::apply(span_.extents(), std::forward<Func>(f));
    }
  }

  // --------------------------------------------------------------------------
  // Compute stride-sorted iteration order. For layout_stride, sorts dimensions
  // by decreasing stride (so the smallest stride is the innermost loop). For
  // other layouts, returns a dummy value (unused but must compile).
  // --------------------------------------------------------------------------

  static std::array<size_t, Rank> computeStrideOrder(const mdspan_type& span)
  {
    std::array<size_t, Rank> order{};
    if constexpr (Rank == 0) { return order; }
    else if constexpr (std::is_same_v<layout_type, std::layout_stride>)
    {
      return impl::strideOrder<Rank>(span.mapping());
    }
    else
    {
      // Unused for compile-time dispatch, but must be initialized.
      std::iota(order.begin(), order.end(), size_t{0});
      return order;
    }
  }

  // --------------------------------------------------------------------------
  // Product of all extents -> logical element count.
  //
  // Uses the index_sequence technique to expand s.extent(0), s.extent(1), ...,
  // s.extent(Rank-1) and multiplies them with a fold expression: (s.extent(Is)
  // * ...) -> s.extent(0) * s.extent(1) * ... * s.extent(N-1)
  //
  // The two-function pattern (computeSize calls computeSizeImpl with
  // make_index_sequence) exists because you can't "open" an index_sequence
  // directly - the template deduction step in computeSizeImpl converts the type
  // index_sequence<0,1,2> into the usable pack Is = {0, 1, 2}.
  // --------------------------------------------------------------------------

  template<size_t... Is>
  static sunindextype computeSizeImpl(const mdspan_type& s,
                                      std::index_sequence<Is...>)
  { return static_cast<sunindextype>((s.extent(Is) * ...)); }

  static sunindextype computeSize(const mdspan_type& s)
  {
    if constexpr (Rank == 0) return 1;
    else return computeSizeImpl(s, std::make_index_sequence<Rank>{});
  }
};

} // namespace example
} // namespace sundials

#endif // NVECTOR_MDSPAN_HPP
