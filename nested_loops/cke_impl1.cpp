#include "cke.hpp"
#include "cke_impl.hpp"

//sec C++/Kokkos/EKAT demo implementation 1.

namespace cke {

void run (const Data& d) {
  
}

} // namespace cke

void cke_impl1_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  cke::run(*d);
}

