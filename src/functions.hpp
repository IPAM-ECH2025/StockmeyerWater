#pragma once

#include "parameters.hpp"
#include <Kokkos_Core.hpp>

template <typename RealType, typename ViewType>
KOKKOS_FUNCTION RealType compute_integral(const ViewType &f) {
  // NOTE: This only works for 3D views
  // NOTE: This only works for uniform grids
  using P = Parameters<RealType>;

  RealType sum = 0.0;

  Kokkos::parallel_reduce(
      "compute_integral",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {P::N_y, P::N_theta, P::N_omega}),
      KOKKOS_LAMBDA(const int i, const int j, const int k,
                    RealType &local_sum) { local_sum += f(i, j, k); },
      sum);

  return sum * P::dy * P::dtheta * P::domega;
}
