#pragma once

#include <Kokkos_Core.hpp>
#include <array>

template <typename RealType>
KOKKOS_INLINE_FUNCTION std::array<RealType, 2> angle_to_vec(RealType theta) {
  return {Kokkos::cos(theta), Kokkos::sin(theta)};
}

template <typename RealType>
KOKKOS_INLINE_FUNCTION std::array<RealType, 2> angle_to_perp(RealType theta) {
  return {-Kokkos::sin(theta), Kokkos::cos(theta)};
}

template <typename RealType>
KOKKOS_INLINE_FUNCTION RealType japanese_bracket(RealType x, RealType delta) {
  return Kokkos::sqrt(x * x + delta * delta);
}
