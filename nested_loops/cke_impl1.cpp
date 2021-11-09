#include "cke.hpp"
#include "cke_impl.hpp"

//sec C++/Kokkos/EKAT demo implementation 1.

namespace cke {

void qad (const Data& d) {
  const auto tracerCur = scalarize(d.tracerCur);
  const auto cellMask = scalarize(d.cellMask);
  const auto normalThicknessFlux = scalarize(d.normalThicknessFlux);
  const auto advMaskHighOrder = scalarize(d.advMaskHighOrder);
  const auto highOrderFlx = scalarize(d.highOrderFlx);
  for (int iEdge = 0; iEdge < d.nEdges; ++iEdge)
    for (int k = 0; k < d.nVertLevels; ++k) {
      const auto coef2 = normalThicknessFlux(iEdge,k)*advMaskHighOrder(iEdge,k);
      const auto csgn = normalThicknessFlux(iEdge,k) < 0 ? -1 : 1;
      Real edgeFlx = 0;
      for (int i = 0; i < d.nAdvCellsForEdge(iEdge); ++i) {
        const auto iCell = d.advCellsForEdge(iEdge,i);
        // not necessary: if (k < d.minLevelCell(iCell) || k > d.maxLevelCell(iCell)) continue;
        const auto coef1 = d.advCoefs(iEdge,i);
        const auto coef3 = d.advCoefs3rd(iEdge,i)*d.coef3rdOrder;
        edgeFlx += tracerCur(iCell,k)*cellMask(iCell,k)*coef2*(coef1 + coef3*csgn);
      }
      highOrderFlx(iEdge,k) = edgeFlx;
    }
}

void run (const Data& d) {
}

} // namespace cke

void cke_impl1_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  cke::qad(*d);
  for (int iter = 0; iter < d->nIters; ++iter)
    cke::run(*d);
}
