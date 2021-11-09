/* This file implements one approach to using C++/Kokkos/EKAT to implemented the
   nested loop kernel.

   See comments beginning with "//sec" to see each section of code: data
   structures, test administrative details, and the actual demo.
 */

#include "cpp_impl1.hpp"

#include <memory>

#include "Kokkos_Core.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

//sec Kokkos init/finalize

static bool in_charge_of_kokkos = false;

void kokkos_init () {
  if (Kokkos::is_initialized()) return;
  in_charge_of_kokkos = true;
  std::vector<char*> args;
#ifdef KOKKOS_ENABLE_CUDA
  // Round-robin assignment of ranks to GPUs.
  int nd;
  const auto ret = cudaGetDeviceCount(&nd);
  if (ret != cudaSuccess) {
    // It isn't a big deal if we can't get the device count.
    nd = 1;
  }
  std::stringstream ss;
  ss << "--kokkos-ndevices=" << nd;
  const auto key = ss.str();
  std::vector<char> str(key.size()+1);
  std::copy(key.begin(), key.end(), str.begin());
  str.back() = 0;
  args.push_back(const_cast<char*>(str.data()));
#endif
  const char* silence = "--kokkos-disable-warnings";
  args.push_back(const_cast<char*>(silence));
  int narg = args.size();
  Kokkos::initialize(narg, args.data());
}

void kokkos_finalize () {
  if (in_charge_of_kokkos && Kokkos::is_initialized())
    Kokkos::finalize();
}

//sec Data structures

struct Data {
  typedef std::shared_ptr<Data> Ptr;

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

  void init(
    const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
    const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
    const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
    const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
    const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder);  
};

//sec Test administrative details

Data::Ptr g_data;

template <typename Scalar, typename V> static
void init1d (const Scalar* raw, const int d1, const std::string& name, V& v,
             const Scalar delta = 0) {
  const auto vnc = typename V::non_const_type(name, d1);
  const auto h = Kokkos::create_mirror_view(vnc);
  for (int i = 0; i < d1; ++i) h(i) = raw[i] + delta;
  Kokkos::deep_copy(vnc, h);
  v = vnc;
}

void Data::init (
  const Int nIters_, const Int nEdges_, const Int nCells_, const Int nVertLevels_,
  const Int nAdv_, const Int* nAdvCellsForEdge_, const Int* minLevelCell_,
  const Int* maxLevelCell_, const Int* advCellsForEdge_, const Real* tracerCur_,
  const Real* normalThicknessFlux_, const Real* advMaskHighOrder_, const Real* cellMask_,
  const Real* advCoefs_, const Real* advCoefs3rd_, const Real coef3rdOrder_)
{
  nIters = nIters_; nEdges = nEdges_; nCells = nCells_; nVertLevels = nVertLevels_;
  nAdv = nAdv_; coef3rdOrder = coef3rdOrder_;

  const int npack = ekat::PackInfo<Data::packn>::num_packs(Data::packn);
  init1d(nAdvCellsForEdge_, nEdges, "nAdvCellsForEdge", nAdvCellsForEdge);
}

void cpp_impl1_init (
  const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
  const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
  const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
  const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
  const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder)
{
  g_data = std::make_shared<Data>();
  g_data->init(nIters, nEdges, nCells, nVertLevels, nAdv,
               nAdvCellsForEdge, minLevelCell, maxLevelCell, advCellsForEdge,
               tracerCur, normalThicknessFlux, advMaskHighOrder, cellMask,
               advCoefs, advCoefs3rd, coef3rdOrder);
}

static void run(const Data& d);

void cpp_impl1_run () {
  assert(g_data);
  run(*g_data);
}

void cpp_impl1_get_results (const Int nEdges, const Int nVertLevels,
                            Real* highOrderFlx) {
  
}

void cpp_impl1_cleanup () { g_data = nullptr; }

//sec C++/Kokkos/EKAT demo implementation 1.

void run (const Data& d) {
  
}
