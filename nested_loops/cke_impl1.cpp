#include "cke.hpp"
#include "cke_impl.hpp"

//sec C++/Kokkos/EKAT demo implementation 1.

namespace cke {

void run (const Data& d) {
  const auto nvlpk = d.highOrderFlx.extent_int(1);
#if 1
  const auto normalThicknessFlux_s = scalarize(d.normalThicknessFlux);
  for (int iEdge = 0; iEdge < d.nEdges; ++iEdge) {
    const auto iend = d.nAdvCellsForEdge(iEdge);
    for (int k = 0; k < nvlpk; ++k) {
      const auto coef2 = d.normalThicknessFlux(iEdge,k)*d.advMaskHighOrder(iEdge,k);
# if 1
      Data::Pr csgn;
      vector_simd for (int s = 0; s < Data::packn; ++s)
        csgn[s] = normalThicknessFlux_s(iEdge,k*Data::packn+s) < 0 ? -1 : 1;
# else
      Data::Pr csgn(1);
      csgn.set(d.normalThicknessFlux(iEdge,k) < 0, -1);
# endif
      Data::Pr edgeFlx = 0;
      for (int i = 0; i < iend; ++i) {
        const auto coef1 = d.advCoefs(iEdge,i);
        const auto coef3 = d.advCoefs3rd(iEdge,i)*d.coef3rdOrder;
        const auto iCell = d.advCellsForEdge(iEdge,i);
        edgeFlx += d.tracerCur(iCell,k)*d.cellMask(iCell,k)*coef2*(coef1 + coef3*csgn);
      }
      d.highOrderFlx(iEdge,k) = edgeFlx;
    }
  }
#else
  const auto normalThicknessFlux_s = scalarize(d.normalThicknessFlux);
  Data::Pr edgeFlx[128/Data::packn];
  for (int iEdge = 0; iEdge < d.nEdges; ++iEdge) {
    Data::Pr coef2[128/Data::packn];
    for (int k = 0; k < nvlpk; ++k)
      coef2[k] = d.normalThicknessFlux(iEdge,k)*d.advMaskHighOrder(iEdge,k);
    Data::Pr csgn[128/Data::packn];
    Real* scsgn = reinterpret_cast<Real*>(&csgn[0]);
    for (int k = 0; k < d.nVertLevels; ++k)
      scsgn[k] = normalThicknessFlux_s(iEdge,k) < 0 ? -1 : 1;
    for (int k = 0; k < nvlpk; ++k)
      edgeFlx[k] = 0;
    const auto iend = d.nAdvCellsForEdge(iEdge);
    for (int i = 0; i < iend; ++i) {
      const auto coef1 = d.advCoefs(iEdge,i);
      const auto coef3 = d.advCoefs3rd(iEdge,i)*d.coef3rdOrder;
      const auto iCell = d.advCellsForEdge(iEdge,i);
      const auto kbeg = d.minLevelCell(iCell)/Data::packn;
      const auto kend = d.maxLevelCell(iCell)/Data::packn;
      for (int k = kbeg; k <= kend; ++k)
        edgeFlx[k] += d.tracerCur(iCell,k)*d.cellMask(iCell,k)*coef2[k]*(coef1 + coef3*csgn[k]);
    }
    for (int k = 0; k < nvlpk; ++k)
      d.highOrderFlx(iEdge,k) = edgeFlx[k];
  }
#endif
}

} // namespace cke

void cke_impl1_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  for (int iter = 0; iter < d->nIters; ++iter)
    cke::run(*d);
}
