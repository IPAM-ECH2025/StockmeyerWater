#pragma once

#include <cmath>
#include <numbers>

template <typename RealType> struct Parameters {
  /*
   * Spatial discretization parameters
   */
  static constexpr int N_x = 128; // Corresponds to y
  static constexpr int N_y = 128; // Corresponds to \theta
  static constexpr int N_z = 128; // Corresponds to \omega
  static constexpr RealType L_x_lower = -std::numbers::pi;
  static constexpr RealType L_y_lower = -std::numbers::pi;
  static constexpr RealType L_z_lower = -std::numbers::pi;
  static constexpr RealType L_x_upper = std::numbers::pi;
  static constexpr RealType L_y_upper = std::numbers::pi;
  static constexpr RealType L_z_upper = std::numbers::pi;

  static constexpr RealType dx = (L_x_upper - L_x_lower) / RealType(N_x);
  static constexpr RealType dy = (L_y_upper - L_y_lower) / RealType(N_y);
  static constexpr RealType dz = (L_z_upper - L_z_lower) / RealType(N_z);

  static constexpr int n_points = N_x * N_y * N_z;
  static constexpr int n_cells = (N_x - 1) * (N_y - 1) * (N_z - 1);

  /**
   * Temporal discretization parameters
   */
  static constexpr RealType dt = 0.005;
  static constexpr RealType t_final = 100 * dt; // 10.0;

  static constexpr auto n_steps = (int)std::ceil(t_final / dt);

  /**
   * Initial condition parameters
   */
  static constexpr RealType k = 2.0;
  static constexpr RealType sigma_star = 0.3;

  /**
   * Model parameters
   */
  static constexpr RealType gamma_3 = 22.0;
  static constexpr RealType gamma_4 = 28.0;
  static constexpr RealType E_x = -1.0;
  static constexpr RealType E_y = 0.0;
  static constexpr RealType delta = 1.0; // Japanese bracket regularization
};
