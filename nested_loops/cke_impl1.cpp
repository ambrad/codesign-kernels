#include "cke.hpp"
#include "cke_impl.hpp"

//sec C++/Kokkos/EKAT demo implementation 1.

namespace cke {

void run (const Data& d) {
  const auto nvlpk = d.highOrderFlx.extent_int(1);
  for (int iEdge = 0; iEdge < d.nEdges; ++iEdge)
    for (int k = 0; k < nvlpk; ++k) {
      const auto coef2 = d.normalThicknessFlux(iEdge,k)*d.advMaskHighOrder(iEdge,k);
      Data::Pr csgn(1);
      csgn.set(d.normalThicknessFlux(iEdge,k) < 0, -1);
      Data::Pr edgeFlx = 0;
      for (int i = 0; i < d.nAdvCellsForEdge(iEdge); ++i) {
        const auto iCell = d.advCellsForEdge(iEdge,i);
        const auto coef1 = d.advCoefs(iEdge,i);
        const auto coef3 = d.advCoefs3rd(iEdge,i)*d.coef3rdOrder;
        edgeFlx += d.tracerCur(iCell,k)*d.cellMask(iCell,k)*coef2*(coef1 + coef3*csgn);
      }
      d.highOrderFlx(iEdge,k) = edgeFlx;
    }
}

} // namespace cke

void cke_impl1_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  for (int iter = 0; iter < d->nIters; ++iter)
    cke::run(*d);
}
