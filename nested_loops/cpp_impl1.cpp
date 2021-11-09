#include "cpp_impl1.hpp"
#include "Kokkos_Core.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

struct Data {
  enum : int { packn = CPP_IMPL1_PACK_SIZE };

  typedef Kokkos::LayoutRight Layout;

  typedef Kokkos::View<Int*,Layout> Ai1;
  typedef Kokkos::View<const Int*,Layout> Aci1;
  typedef Kokkos::View<Int**,Layout> Ai2;
  typedef Kokkos::View<const Int**,Layout> Aci2;

  typedef Kokkos::View<Real*,Layout> Ar1;
  typedef Kokkos::View<const Real*,Layout> Acr1;
  typedef Kokkos::View<Real**,Layout> Ar2;
  typedef Kokkos::View<const Real**,Layout> Acr2;

  typedef ekat::Pack<Real,packn> Pr;
  typedef ekat::Pack<const Real,packn> Pcr;

  typedef Kokkos::View<Pr*,Layout> Apr1;
  typedef Kokkos::View<Pcr*,Layout> Acpr1;
  typedef Kokkos::View<Pr**,Layout> Apr2;
  typedef Kokkos::View<Pcr**,Layout> Acpr2;

  Int nIters, nEdges, nCells, nVertLevels, nAdv;
  Real coef3rdOrder;
  Aci1 nAdvCellsForEdge, minLevelCell, maxLevelCell;
  Aci2 advCellsForEdge;
  Acpr2 tracerCur, cellMask, normalThicknessFlux, advMaskHighOrder, advCoefs, advCoefs3rd;

  void init (
    const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
    const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
    const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
    const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
    const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder)
  {
    
  }
  
};

void cpp_impl1_init (
  const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
  const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
  const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
  const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
  const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder)
{
  
}

void cpp_impl1_run () {
  
}

void cpp_impl1_get_results (
  const Int nEdges, const Int nVertLevels, Real* highOrderFlx)
{
  
}
