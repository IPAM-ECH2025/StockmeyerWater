#pragma once

#include "parameters.hpp"
#include <Kokkos_Core.hpp>
#include <array>

// TODO: Use this
template <typename RealType>
KOKKOS_INLINE_FUNCTION std::array<RealType, 2> angle_to_vec(RealType theta) {
  using Kokkos::cos;
  using Kokkos::sin;
  return {cos(theta), sin(theta)};
}

// TODO: Use this
template <typename RealType>
KOKKOS_INLINE_FUNCTION std::array<RealType, 2> angle_to_perp(RealType theta) {
  using Kokkos::cos;
  using Kokkos::sin;
  return {-sin(theta), cos(theta)};
}

template <typename RealType>
KOKKOS_INLINE_FUNCTION RealType japanese_bracket(RealType x) {
  return Kokkos::sqrt(x * x + Parameters<RealType>::delta *
                                  Parameters<RealType>::delta);
}

template <typename RealType, typename ViewType>
KOKKOS_FUNCTION RealType compute_integral(const ViewType &f) {
  // NOTE: This only works for 3D views
  // NOTE: This only works for uniform grids
  RealType sum = 0.0;

  Kokkos::parallel_reduce(
      "compute_integral",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          {0, 0, 0}, {Parameters<RealType>::N_x, Parameters<RealType>::N_y,
                      Parameters<RealType>::N_z}),
      KOKKOS_LAMBDA(const int &i, const int &j, const int &k,
                    RealType &local_sum) { local_sum += f(i, j, k); },
      sum);

  return sum * Parameters<RealType>::dx * Parameters<RealType>::dy *
         Parameters<RealType>::dz;
}
