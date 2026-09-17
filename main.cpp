#include "KokkosFFT_Common_Types.hpp"
#include "src/functions.hpp"
#include "src/output.hpp"
#include "src/parameters.hpp"
#include <KokkosFFT.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>
#include <filesystem>
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

  using P = Parameters<Real>;

  WaterProblem()
      : f("f", P::N_x, P::N_y, P::N_z), x("x", P::N_x), y("y", P::N_y),
        z("z", P::N_z),
        f_hat_theta("f_hat_theta", P::N_x, P::N_y / 2 + 1, P::N_z),
        f_hat_omega("f_hat_omega", P::N_x, P::N_y, P::N_z / 2 + 1),
        k_theta(KokkosFFT::rfftfreq(Exec{}, std::size_t{P::N_y}, Real{P::dy})),
        k_omega(KokkosFFT::rfftfreq(Exec{}, std::size_t{P::N_z}, Real{P::dz})),
        acceleration("acceleration", P::N_x, P::N_y) {
    Kokkos::parallel_for(
        "initial_condition_x", P::N_x,
        KOKKOS_LAMBDA(const int i) { x(i) = P::L_x_lower + i * P::dx; });
    Kokkos::parallel_for(
        "initial_condition_y", P::N_y,
        KOKKOS_LAMBDA(const int i) { y(i) = P::L_y_lower + i * P::dy; });
    Kokkos::parallel_for(
        "initial_condition_z", P::N_z,
        KOKKOS_LAMBDA(const int i) { z(i) = P::L_z_lower + i * P::dz; });

    Kokkos::parallel_for(
        "scale_k_theta", P::N_y / 2 + 1,
        KOKKOS_LAMBDA(const int i) { k_theta(i) *= 2.0 * M_PI; });
    Kokkos::parallel_for(
        "scale_k_omega", P::N_z / 2 + 1,
        KOKKOS_LAMBDA(const int i) { k_omega(i) *= 2.0 * M_PI; });

    // TODO: Should add check for k-space and f_hat to have same dimensions
  }

  void run() {
    initial_condition();

    Real time = 0.0;

    for (int step = 0; step <= P::n_steps; ++step) {
      print_solution(step, time);

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
   *
   * @note Remember x->y y->\theta z->\omega
   */
  KOKKOS_INLINE_FUNCTION Real f_0(Real x, Real y, Real z) const {
    // TODO: Change x, y, and z here to y, theta, omega
    using Kokkos::exp;
    using Kokkos::sin;

    const auto _y = x;
    const auto _theta = y;
    const auto _omega = z;
    /*
        const auto sigma = P::sigma_star * (1.0 + 0.5 * sin(P::k * _y));

        return 1.0 / (sqrt(2.0 * M_PI) * sigma) *
               exp(-_omega * _omega / (2.0 * sigma * sigma));
    */
    return exp(-8.0 * _omega * _omega) *
           exp(-8.0 * (_theta - M_PI / 4.0) * (_theta - M_PI / 4.0));
  }

  void initial_condition() {
    // Apply the provided f_0 function
    Kokkos::parallel_for(
        "initial_condition_initial",
        Policy({0, 0, 0}, {P::N_x, P::N_y, P ::N_z}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          const auto x_ = x(i);
          const auto y_ = y(j);
          const auto z_ = z(k);

          f(i, j, k) = f_0(x_, y_, z_);
        });

    // Compute the integral of the density and renormalize
    const auto integral_f = compute_integral<Real>(f);
    Kokkos::parallel_for(
        "initial_condition_renormalize",
        Policy({0, 0, 0}, {P::N_x, P::N_y, P ::N_z}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          f(i, j, k) /= integral_f;
        });
  }

  /**
   * Solve the update for a single step.
   *
   * The update term has two advective components: one in theta and one in
   * omega. To make things easier, we'll split the update into two semigroups.
   */
  void compute_acceleration() {
    using Policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;

    Kokkos::parallel_for(
        "compute_acceleration", Policy2D({0, 0}, {P::N_x, P::N_y}),
        KOKKOS_LAMBDA(const int i, const int j) {
          using Kokkos::cos;
          using Kokkos::sin;

          const Real theta = y(j);

          const Real E_dot_n_perp = -P::E_x * sin(theta) + P::E_y * cos(theta);

          acceleration(i, j) = -P::gamma_4 * E_dot_n_perp;
        });
  }
  void increment(Real dt) {
    // Compute the prefactors
    compute_acceleration();

    // TODO: I think the indexing is off

    // Theta advection
    KokkosFFT::rfft(Exec{}, f, f_hat_theta, KokkosFFT::Normalization::none, 1);
    Kokkos::parallel_for(
        "apply_theta_semigroup",
        Policy({0, 0, 0}, {P::N_x, P::N_y / 2 + 1, P::N_z}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          using Kokkos::cos;
          using Kokkos::sin;

          const Real phase = -k_theta(j) * z(k) * dt;
          f_hat_theta(i, j, k) *= Complex(cos(phase), sin(phase));
        });
    KokkosFFT::irfft(Exec{}, f_hat_theta, f, KokkosFFT::Normalization::backward,
                     1);

    // Omega advection
    KokkosFFT::rfft(Exec{}, f, f_hat_omega, KokkosFFT::Normalization::none, 2);
    Kokkos::parallel_for(
        "apply_omega_semigroup",
        Policy({0, 0, 0}, {P::N_x, P::N_y, P::N_z / 2 + 1}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
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
  void print_solution(int step, Real time) const {

    std::cout << "Integrated value of f " << compute_integral<Real>(f)
              << std::endl;

    std::ostringstream filename;
    filename << "solution-" << std::setw(6) << std::setfill('0') << step
             << ".vtr";

    write_vtu<Real>(f, x, y, z, filename.str(), time, step);
  }

  /**
   * Data objects
   */
  R3 f; // probability density
  R1 x; // x-coordinates
  R1 y; // y-coordinates
  R1 z; // z-coordinates

  Z3 f_hat_theta;
  Z3 f_hat_omega;
  R1 k_theta;
  R1 k_omega;
  R2 acceleration;
};

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    WaterProblem problem;
    problem.run();
  }
  Kokkos::finalize();
}
