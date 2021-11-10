/* General setup, data init, and cleanup for any C++/Kokkos implementation.
 */

#ifndef INCLUDE_CKE_IMPL_HPP
#define INCLUDE_CKE_IMPL_HPP

#include <memory>

#include "Kokkos_Core.hpp"

#include "ekat/ekat_session.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

namespace cke {

struct Data {
  typedef std::shared_ptr<Data> Ptr; // convenience: Data::Ptr

  enum : int { packn = CKE_PACK_SIZE }; // pack size
  typedef ekat::Pack<Real,packn> Pr; // shorthand for a pack of reals

  typedef Kokkos::DefaultExecutionSpace ExeSpace;
  typedef Kokkos::LayoutRight Layout; // data layout; can switch to experiment

  template <typename Data>
  using View = Kokkos::View<Data,Layout,ExeSpace,Kokkos::MemoryTraits<Kokkos::Restrict>>;

  // Some handy view aliases.
  typedef View<Int*> Ai1;
  typedef View<const Int*> Aci1;
  typedef View<Int**> Ai2;
  typedef View<const Int**> Aci2;

  typedef View<Real*> Ar1;
  typedef View<const Real*> Acr1;
  typedef View<Real**> Ar2;
  typedef View<const Real**> Acr2;

  typedef View<Pr*> Apr1;
  typedef View<const Pr*> Acpr1;
  typedef View<Pr**> Apr2;
  typedef View<const Pr**> Acpr2;

  // Read-only data from F90.
  Int nIters, nEdges, nCells, nVertLevels, nAdv;
  Real coef3rdOrder;
  Aci1 nAdvCellsForEdge, minLevelCell, maxLevelCell;
  Aci2 advCellsForEdge;
  Acr2 advCoefs, advCoefs3rd;
  Acpr2 tracerCur, cellMask, normalThicknessFlux, advMaskHighOrder;

  // Output.
  Apr2 highOrderFlx;

  void init(
    const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
    const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
    const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
    const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
    const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder);  
};

Data::Ptr get_Data_singleton();

template <typename Fn, typename ExeSpace = Kokkos::DefaultExecutionSpace>
void parfor_iEdge_kPack (const Data& d, const Fn& f) {
  const auto nvlpk = ekat::PackInfo<Data::packn>::num_packs(d.nVertLevels);
  if (ekat::OnGpu<ExeSpace>::value) {
    const auto p = Kokkos::RangePolicy<ExeSpace>(0, d.nEdges*nvlpk);
    const auto g = KOKKOS_LAMBDA(const int idx) {
      const int iEdge = idx / nvlpk, k = idx % nvlpk;
      f(iEdge, k);
    };
    Kokkos::parallel_for(p, g);
  } else {
#   pragma omp parallel for
    for (int iEdge = 0; iEdge < d.nEdges; ++iEdge)
      for (int k = 0; k < nvlpk; ++k)
        f(iEdge, k);
  }
}

} // namespace cke

#endif
