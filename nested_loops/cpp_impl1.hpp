#ifndef INCLUDE_CPP_IMPL1_HPP
#define INCLUDE_CPP_IMPL1_HPP

typedef int Int;
typedef double Real;

extern "C" {
  void kokkos_init();
  void kokkos_finalize();

  void cpp_impl1_init(
    const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
    const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
    const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
    const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
    const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder);
  void cpp_impl1_run();
  void cpp_impl1_get_results(const Int nEdges, const Int nVertLevels, Real* highOrderFlx);
  void cpp_impl1_cleanup();
}

#endif
