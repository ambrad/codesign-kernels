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

  // Some handy view aliases.
  typedef Kokkos::View<Int*,Layout> Ai1;
  typedef Kokkos::View<const Int*,Layout> Aci1;
  typedef Kokkos::View<Int**,Layout> Ai2;
  typedef Kokkos::View<const Int**,Layout> Aci2;

  typedef Kokkos::View<Real*,Layout> Ar1;
  typedef Kokkos::View<const Real*,Layout> Acr1;
  typedef Kokkos::View<Real**,Layout> Ar2;
  typedef Kokkos::View<const Real**,Layout> Acr2;

  typedef Kokkos::View<Pr*,Layout> Apr1;
  typedef Kokkos::View<const Pr*,Layout> Acpr1;
  typedef Kokkos::View<Pr**,Layout> Apr2;
  typedef Kokkos::View<const Pr**,Layout> Acpr2;

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
