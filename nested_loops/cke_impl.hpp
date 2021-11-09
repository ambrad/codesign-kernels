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

  typedef Kokkos::LayoutRight Layout; // data layout; can switch to experiment

  template <typename Data>
  using View = Kokkos::View<Data,Layout,Kokkos::MemoryTraits<Kokkos::Restrict>>;

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

} // namespace cke

#endif
