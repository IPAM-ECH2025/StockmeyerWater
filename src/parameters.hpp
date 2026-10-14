#pragma once

#include <cmath>
#include <numbers>

template <typename RealType> struct Parameters {
  /*
   * Spatial discretization parameters
   */
  static constexpr int N_y = 64;
  static constexpr int N_theta = 128;
  static constexpr int N_omega = 128;
  static constexpr RealType L_y_lower = -M_PI;
  static constexpr RealType L_theta_lower = -M_PI;
  static constexpr RealType L_omega_lower = -5.0 * M_PI;
  static constexpr RealType L_y_upper = M_PI;
  static constexpr RealType L_theta_upper = M_PI;
  static constexpr RealType L_omega_upper = 5.0 * M_PI;

  static constexpr RealType dy = (L_y_upper - L_y_lower) / RealType(N_y);
  static constexpr RealType dtheta =
      (L_theta_upper - L_theta_lower) / RealType(N_theta);
  static constexpr RealType domega =
      (L_omega_upper - L_omega_lower) / RealType(N_omega);

  static constexpr int n_points = N_y * N_theta * N_omega;
  static constexpr int n_cells = (N_y - 1) * (N_theta - 1) * (N_omega - 1);

  /**
   * Temporal discretization parameters
   */
  static constexpr RealType dt = 0.005;
  static constexpr RealType t_final = 10.0;

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
  static constexpr RealType E_magnitude = 0.2;
  static constexpr RealType E_angle = M_PI / 4.0;
  static constexpr RealType delta = 1.0; // Japanese bracket regularization

  static constexpr RealType E_x = E_magnitude * std::cos(E_angle);
  static constexpr RealType E_y = E_magnitude * std::sin(E_angle);

  /**
   * Output parameters
   */
  static constexpr int output_stride = 10;
};
