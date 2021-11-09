/* This file implements one or more approaches to using C++/Kokkos/EKAT to
   implemented the nested loop kernel.

   See comments beginning with "//sec" to see each section of code: data
   structures, test administrative details, and the actual demo.
 */

#include "cpp_impl.hpp"

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
  typedef std::shared_ptr<Data> Ptr; // convenience: Data::Ptr

  enum : int { packn = CPP_IMPL_PACK_SIZE }; // pack size
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
  Acpr2 tracerCur, cellMask, normalThicknessFlux, advMaskHighOrder, advCoefs, advCoefs3rd;

  // Output.
  Apr2 highOrderFlx;

  void init(
    const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
    const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
    const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
    const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
    const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder);  
};

//sec Administrative details

// Initialize a View<const Pack<Real,packn>**> from raw(1:d1,1:d2), where dim 2
// has the fast index.
template <typename Scalar, typename V> static
void initvpk (const Scalar* raw, const int d1, const int d2, const int packn,
              const std::string& name, V& v) {
  // Get the number of packs that cover the scalar length.
  const int d2pk = ekat::PackInfo<Data::packn>::num_packs(d2);
  // Allocate the view as writeable.
  const auto vnc = typename V::non_const_type(name, d1, d2pk);
  // For convenience, take a scalar view of the original Pack<Scalar> view.
  const auto svnc = scalarize(vnc);
  // Copy the F90 data to a host mirror of the scalar view.
  const auto h = Kokkos::create_mirror_view(svnc);
  for (int i = 0; i < d1; ++i)
    for (int j = 0; j < d2; ++j)
      h(i,j) = raw[d2*i+j];
  // Copy the data to device.
  Kokkos::deep_copy(svnc, h);
  // Set the possibly read-only view with this data.
  v = vnc;
}

// Simpler forms of the above: scalar 1D and 2D views.

template <typename Scalar, typename V> static
void initv (const Scalar* raw, const int d1, const std::string& name, V& v,
            const Scalar delta = 0) {
  const auto vnc = typename V::non_const_type(name, d1);
  const auto h = Kokkos::create_mirror_view(vnc);
  for (int i = 0; i < d1; ++i) h(i) = raw[i] + delta;
  Kokkos::deep_copy(vnc, h);
  v = vnc;
}

template <typename Scalar, typename V> static
void initv (const Scalar* raw, const int d1, const int d2, const std::string& name,
            V& v, const Scalar delta = 0) {
  const auto vnc = typename V::non_const_type(name, d1, d2);
  const auto h = Kokkos::create_mirror_view(vnc);
  for (int i = 0; i < d1; ++i)
    for (int j = 0; j < d2; ++j)
      h(i,j) = raw[d2*i+j] + delta;
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

  initv(nAdvCellsForEdge_, nEdges, "nAdvCellsForEdge", nAdvCellsForEdge);
  initv(minLevelCell_, nCells, "minLevelCell", minLevelCell, -1);
  initv(maxLevelCell_, nCells, "maxLevelCell", maxLevelCell, -1);
  initv(advCellsForEdge_, nEdges, nAdv, "advCellsForEdge", advCellsForEdge, -1);
  initvpk(tracerCur_, nCells, nVertLevels, packn, "tracerCur", tracerCur);
  initvpk(cellMask_, nCells, nVertLevels, packn, "cellMask", cellMask);
  initvpk(normalThicknessFlux_, nEdges, nVertLevels, packn,
          "normalThicknessFlux", normalThicknessFlux);
  initvpk(advMaskHighOrder_, nEdges, nVertLevels, packn,
          "advMaskHighOrder", advMaskHighOrder);
  initv(advCoefs_, nEdges, nAdv, "advCoefs", advCoefs);
  initv(advCoefs3rd_, nEdges, nAdv, "advCoefs3rd", advCoefs3rd);

  const int npack = ekat::PackInfo<packn>::num_packs(nVertLevels);
  highOrderFlx = Apr2("highOrderFlx", nEdges, npack);
}

//sec API

Data::Ptr g_data;

void cpp_impl_init (
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

void cpp_impl_get_results (const Int nEdges, const Int nVertLevels,
                           Real* highOrderFlx) {
  
}

void cpp_impl_cleanup () { g_data = nullptr; }

//sec C++/Kokkos/EKAT demo implementation 1.

void run_impl1 (const Data& d) {
  
}
void cpp_impl1_run () {
  assert(g_data);
  run_impl1(*g_data);
}
