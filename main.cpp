#include "Kokkos_Array.hpp"
#include "impl/Kokkos_Profiling.hpp"
#include "src/functions.hpp"
#include "src/output.hpp"
#include "src/parameters.hpp"
#include <KokkosFFT.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>
#include <sstream>

class WaterProblem {
public:
  using Exec = Kokkos::DefaultExecutionSpace;
  using Real = double;
  using Complex = Kokkos::complex<Real>;

  using R1 = Kokkos::View<Real *, Kokkos::LayoutRight, Exec>;
  using R2 = Kokkos::View<Real **, Kokkos::LayoutRight, Exec>;
  using R3 = Kokkos::View<Real ***, Kokkos::LayoutRight, Exec>;
  using Z1 = Kokkos::View<Complex *, Kokkos::LayoutRight, Exec>;
  using Z2 = Kokkos::View<Complex **, Kokkos::LayoutRight, Exec>;
  using Z3 = Kokkos::View<Complex ***, Kokkos::LayoutRight, Exec>;

  using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
  using Policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
  using TeamPolicy = Kokkos::TeamPolicy<>;
  using MemberType = TeamPolicy::member_type;

  using P = Parameters<Real>;

  struct Vec2D {
    Real x;
    Real y;
  };

  WaterProblem()
      : f("f", P::N_y, P::N_theta, P::N_omega), f_bar("f_bar", P::N_y, 2),
        G("G", P::N_y, 2),

        y("y", P::N_y), theta("theta", P::N_theta), omega("omega", P::N_omega),
        f_hat_theta("f_hat_theta", P::N_y, P::N_theta / 2 + 1, P::N_omega),
        f_hat_omega("f_hat_omega", P::N_y, P::N_theta, P::N_omega / 2 + 1),
        k_theta(KokkosFFT::rfftfreq(Exec{}, P::N_theta, P::dtheta)),
        k_omega(KokkosFFT::rfftfreq(Exec{}, P::N_omega, P::domega)),
        acceleration("acceleration", P::N_y, P::N_theta) {
    Kokkos::parallel_for(
        "initial_condition_y", P::N_y,
        KOKKOS_CLASS_LAMBDA(const int i) { y(i) = P::L_y_lower + i * P::dy; });
    Kokkos::parallel_for(
        "initial_condition_theta", P::N_theta,
        KOKKOS_CLASS_LAMBDA(const int i) {
          theta(i) = P::L_theta_lower + i * P::dtheta;
        });
    Kokkos::parallel_for(
        "initial_condition_omega", P::N_omega,
        KOKKOS_CLASS_LAMBDA(const int i) {
          omega(i) = P::L_omega_lower + i * P::domega;
        });

    Kokkos::parallel_for(
        "scale_k_theta", P::N_theta / 2 + 1,
        KOKKOS_CLASS_LAMBDA(const int i) { k_theta(i) *= 2.0 * M_PI; });
    Kokkos::parallel_for(
        "scale_k_omega", P::N_omega / 2 + 1,
        KOKKOS_CLASS_LAMBDA(const int i) { k_omega(i) *= 2.0 * M_PI; });

    // TODO: Should add check for k-space and f_hat to have same dimensions
  }

  void run() {
    initial_condition();

    Real time = 0.0;

    for (int step = 0; step <= P::n_steps; ++step) {
      if (step % P::output_stride == 0) {
        print_solution(step, time);
      }

      if (step == P::n_steps) {
        break;
      }

      const Real _dt = std::min(P::dt, P::t_final - time);

      increment(_dt);

      time += _dt;
    }
  }

private:
  /**
   * For the initial conditions we use the following function:
   *
   * TODO: Write function
   *
   * Importantly, we normalize the total mass to 1.
   */
  KOKKOS_INLINE_FUNCTION Real f_0(const Real _y, const Real _theta,
                                  const Real _omega) const {
    using Kokkos::exp;
    using Kokkos::sin;

    /*
        const Real sigma = P::sigma_star * (1.0 + 0.5 * sin(P::k * _y));

        return 1.0 / (sqrt(2.0 * M_PI) * sigma) *
               exp(-_omega * _omega / (2.0 * sigma * sigma));
    */
    return exp(-8.0 * _omega * _omega) * exp(-8.0 * _theta * _theta);
  }

  void initial_condition() {
    // Apply the provided f_0 function
    Kokkos::parallel_for(
        "initial_condition_initial",
        Policy({0, 0, 0}, {P::N_y, P::N_theta, P ::N_omega}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k) {
          f(i, j, k) = f_0(y(i), theta(j), omega(k));
        });

    // Compute the integral of the density and renormalize
    const Real integral_f = compute_integral<Real>(f);
    Kokkos::parallel_for(
        "initial_condition_renormalize",
        Policy({0, 0, 0}, {P::N_y, P::N_theta, P ::N_omega}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k) {
          f(i, j, k) /= integral_f;
        });
  }

  /**
   * Solve the update for a single step.
   *
   * The update term has two advective components: one in theta and one in
   * omega. To make things easier, we'll split the update into two semigroups.
   */
  void compute_f_bar() {
    using Kokkos::cos;
    using Kokkos::sin;

    Kokkos::parallel_for(
        "compute_f_bar", TeamPolicy(P::N_y, Kokkos::AUTO),
        KOKKOS_CLASS_LAMBDA(const MemberType &team) {
          int i = team.league_rank();

          Vec2D layer_sum{0.0, 0.0};

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, P::N_theta),
              [=](const int j, Vec2D &local_sum) {
                Real row_sum = 0.0;

                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(team, P::N_omega),
                    [=](const int k, Real &vector_sum) {
                      vector_sum += f(i, j, k);
                    },
                    row_sum);

                const Real _theta = theta(j);

                local_sum.x += row_sum * cos(_theta);
                local_sum.y += row_sum * sin(_theta);
              },
              layer_sum);

          f_bar(i, 0) = layer_sum.x * P::dtheta * P::domega;
          f_bar(i, 1) = layer_sum.y * P::dtheta * P::domega;
        });
    Kokkos::fence();
  }
  KOKKOS_INLINE_FUNCTION Real japanese_bracket(const Real x) const {
    using Kokkos::sqrt;

    return sqrt(x * x + P::delta * P::delta);
  }
  KOKKOS_INLINE_FUNCTION Real periodic_distance(const Real x1,
                                                const Real x2) const {
    using Kokkos::abs;
    using Kokkos::copysign;

    const Real s = x1 - x2;
    static constexpr Real L_half = 0.5 * (P::L_y_upper - P::L_y_lower);

    const Real d = L_half - abs(s);
    return copysign(L_half - abs(d), d);
  }
  KOKKOS_INLINE_FUNCTION Real compute_h_11(const Real r) const {
    static constexpr Real L = P::L_y_upper - P::L_y_lower;
    static constexpr Real L2 = L * L;
    const Real r_b = japanese_bracket(r);
    const Real r_b2 = r_b * r_b;

    return 4.0 * L / (L2 + 4.0 * r_b2);
  }
  KOKKOS_INLINE_FUNCTION Real compute_h_22(const Real r) const {
    using Kokkos::atan;

    static constexpr Real L = P::L_y_upper - P::L_y_lower;
    static constexpr Real L2 = L * L;
    static constexpr Real delta2 = P::delta * P::delta;
    const Real r2 = r * r;
    const Real r_b = japanese_bracket(r);
    const Real r_b2 = r_b * r_b;

    const Real term_1 =
        2.0 * ((delta2 - r2) / r_b2 - 1.0) * L / (L2 + 4.0 * r_b2);
    const Real term_2 =
        ((delta2 - r2) / r_b2 + 1.0) * atan(L / (2.0 * r_b)) / r_b;

    return term_1 + term_2;
  }
  void compute_G() {
    Kokkos::parallel_for(
        "compute_G", P::N_y, KOKKOS_CLASS_LAMBDA(const int i) {
          Real sum_x = 0.0;
          Real sum_y = 0.0;

          for (int j = 0; j < P::N_y; ++j) {
            const Real r = periodic_distance(y(i), y(j));

            sum_x += compute_h_11(r) * f_bar(j, 0);
            sum_y += compute_h_22(r) * f_bar(j, 1);
          }

          G(i, 0) = P::dy * sum_x;
          G(i, 1) = P::dy * sum_y;
        });
    Kokkos::fence();
  }
  void compute_acceleration() {
    Kokkos::parallel_for(
        "compute_acceleration", Policy2D({0, 0}, {P::N_y, P::N_theta}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j) {
          using Kokkos::cos;
          using Kokkos::sin;

          const Real _theta = theta(j);

          const Real G_dot_n_perp =
              G(i, 0) * sin(_theta) + G(i, 1) * cos(_theta);

          const Real E_dot_n_perp =
              -P::E_x * sin(_theta) + P::E_y * cos(_theta);

          acceleration(i, j) =
              P::gamma_3 * G_dot_n_perp - P::gamma_4 * E_dot_n_perp;
        });
    Kokkos::fence();
  }
  void increment(const Real dt) {
    // Compute the prefactors
    compute_f_bar();
    compute_G();
    compute_acceleration();

    // Theta advection
    KokkosFFT::rfft(Exec{}, f, f_hat_theta, KokkosFFT::Normalization::none, 1);
    Kokkos::parallel_for(
        "apply_theta_semigroup",
        Policy({0, 0, 0}, {P::N_y, P::N_theta / 2 + 1, P::N_omega}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k) {
          using Kokkos::cos;
          using Kokkos::sin;

          const Real phase = -k_theta(j) * omega(k) * dt;
          f_hat_theta(i, j, k) *= Complex(cos(phase), sin(phase));
        });
    KokkosFFT::irfft(Exec{}, f_hat_theta, f, KokkosFFT::Normalization::backward,
                     1);

    // Omega advection
    KokkosFFT::rfft(Exec{}, f, f_hat_omega, KokkosFFT::Normalization::none, 2);
    Kokkos::parallel_for(
        "apply_omega_semigroup",
        Policy({0, 0, 0}, {P::N_y, P::N_theta, P::N_omega / 2 + 1}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k) {
          using Kokkos::cos;
          using Kokkos::sin;

          const Real phase = -k_omega(k) * acceleration(i, j) * dt;
          f_hat_omega(i, j, k) *= Complex(cos(phase), sin(phase));
        });
    KokkosFFT::irfft(Exec{}, f_hat_omega, f, KokkosFFT::Normalization::backward,
                     2);
  }

  /**
   * Output solution and print information to terminal
   */
  void print_solution(const int step, const Real time) const {
    std::cout << "Step " << step << " time " << time
              << " integrated value of f " << compute_integral<Real>(f)
              << std::endl;

    std::ostringstream filename;
    filename << "solution-" << std::setw(6) << std::setfill('0') << step
             << ".vtr";

    write_vtu<Real>(f, y, theta, omega, filename.str(), time, step);
  }

  /**
   * Data objects
   */
  R3 f;     // Density [N_y, N_theta, N_omega]
  R2 f_bar; // n(theta)-weighted marginal [N_y, 2]
  R2 G;     // Convolution term [N_y, 2]
  R1 y;     // Distance [N_y]
  R1 theta; // Angle [N_theta]
  R1 omega; // Angular velocity [N_omega]

  Z3 f_hat_theta;
  Z3 f_hat_omega;
  R1 k_theta;
  R1 k_omega;
  R2 acceleration;
};

KOKKOS_INLINE_FUNCTION
WaterProblem::Vec2D &operator+=(WaterProblem::Vec2D &a,
                                const WaterProblem::Vec2D &b) {
  a.x += b.x;
  a.y += b.y;
  return a;
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    WaterProblem problem;
    problem.run();
  }
  Kokkos::finalize();
}
