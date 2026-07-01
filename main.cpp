#include <KokkosFFT.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>
#include <string>

#include "src/functions.hpp"
#include "src/output.hpp"

using Exec = Kokkos::DefaultExecutionSpace;
using Real = double;
using Complex = Kokkos::complex<Real>;

using R1 = Kokkos::View<Real *, Kokkos::LayoutRight, Exec>;
using R2 = Kokkos::View<Real **, Kokkos::LayoutRight, Exec>;
using Z1 = Kokkos::View<Complex *, Kokkos::LayoutRight, Exec>;
using Z2 = Kokkos::View<Complex **, Kokkos::LayoutRight, Exec>;

// Parameters
constexpr int N_x = 128;       // Grid points in x direction
constexpr int N_y = 128;       // Grid points in y direction
constexpr Real dt = 0.001;     // Timestep
constexpr Real t_final = 30.0; // Final time

constexpr Real gamma_3 = 1.0;
constexpr Real gamma_4 = 1.0;
constexpr Real E_x_0 = 0.0;
constexpr Real E_y_0 = 0.0;
constexpr Real ef_angle = 0.0;
constexpr Real sigma_star = 0.3;
constexpr Real k = 2.0;
constexpr Real sine_coefficient = 0.0;
constexpr Real delta = 1.0; // Japanese bracket regularization

constexpr Real L_x_lower = -M_PI; // Grid lower bound in x direction
constexpr Real L_x_upper = M_PI;  // Grid upper bound in x direction
constexpr Real L_y_lower = -M_PI; // Grid lower bound in y direction
constexpr Real L_y_upper = M_PI;  // Grid upper bound in y direction

// Computed values
constexpr Real dx = (L_x_upper - L_x_lower) / Real(N_x);
constexpr Real dy = (L_y_upper - L_y_lower) / Real(N_y);

constexpr int n_points = N_x * N_y;
constexpr int n_cells = (N_x - 1) * (N_y - 1);

// TODO: Document
KOKKOS_INLINE_FUNCTION Real fftfreq(int i, int n, Real h) {
  return Real(i < n / 2 ? i : i - n) / (Real(n) * h);
}

// TODO: Document
KOKKOS_INLINE_FUNCTION Real B_func(Real z) { return -z; }

// Compute the x-semigroup multiplier: exp(i * dt * freq_x * B(y))
// sg_x[i, j] = exp(i * dt * k_x(i) * B(y(j)))
// Note: In MATLAB this is sg_x = exp(1i*dt*freq_x.*By')
//       freq_x is row vector, By is column vector -> outer product
void compute_sg_x(Z2 sg_x, R1 k_x, R1 y) {
  Kokkos::parallel_for(
      "compute_sg_x",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
      KOKKOS_LAMBDA(int i, int j) {
        // B(y(j)) = -y(j)
        Real By_j = B_func(y(j));
        Real phase = dt * k_x(i) * By_j;
        sg_x(i, j) = Complex(Kokkos::cos(phase), Kokkos::sin(phase));
      });
}

// Compute the y-semigroup multiplier: exp(i * dt * A(x,ef,f,c) * freq_y)
// sg_y[i, j] = exp(i * dt * Ax(i) * k_y(j))
// A(x) is computed per row (x-index), freq_y is per column (y-index)
void compute_sg_y(Z2 sg_y, R1 Ax, R1 k_y) {
  Kokkos::parallel_for(
      "compute_sg_y",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
      KOKKOS_LAMBDA(int i, int j) {
        Real phase = dt * Ax(i) * k_y(j);
        sg_y(i, j) = Complex(Kokkos::cos(phase), Kokkos::sin(phase));
      });
}

// Computes Ax(i) for each x-index i
// Matches: return_data = -real_interaction/max(real_interaction) + vec
void compute_Ax(R1 Ax, R1 x, R2 f, R1 c, Real ef_x,
                Real ef_y, // components of ef unit vector
                Real t) {

  const int Nx = N_x;
  const int Ny = N_y;

  // Step 1: no_vel(i) = trapz(f(i,:)) over y  (integral over omega)
  R1 no_vel("no_vel", Nx);
  Kokkos::parallel_for(
      "trapz_y", Nx, KOKKOS_LAMBDA(int i) {
        Real sum = 0.0;
        for (int j = 0; j < Ny - 1; j++) {
          sum += 0.5 * (f(i, j) + f(i, j + 1));
        }
        no_vel(i) = sum * dy;
      });
  Kokkos::fence();

  // Step 2: real_interaction(i) = sum_k c(k) * no_vel(i-k)  (same-mode conv)
  // "same" convolution: output length = max(Nx, len(c)) = Nx here
  R1 real_interaction("real_interaction", Nx);
  Kokkos::parallel_for(
      "conv_same", Nx, KOKKOS_LAMBDA(int i) {
        Real val = 0.0;
        // Full convolution index: out[i] = sum_m a[m] * b[i - m]
        // "same" centers this: offset by (Nx-1)/2
        const int offset = (Nx - 1) / 2;
        for (int m = 0; m < Nx; m++) {
          int idx = i - m + offset;
          if (idx >= 0 && idx < Nx) {
            val += c(m) * no_vel(idx);
          }
        }
        real_interaction(i) = val;
      });
  Kokkos::fence();

  // Step 3: max(real_interaction) - need a reduction
  Real max_val = 0.0;
  Kokkos::Max<Real> reducer(max_val);
  Kokkos::parallel_reduce(
      "max_ri", Nx,
      KOKKOS_LAMBDA(int i, Real &lmax) {
        if (real_interaction(i) > lmax)
          lmax = real_interaction(i);
      },
      reducer);
  Kokkos::fence();

  // Step 4: vec(i) = dot((0 - ef), angle_to_perp(x(i)))
  // angle_to_perp(theta) = (-sin(theta), cos(theta))
  // ef = (ef_x, ef_y)
  // vec = sum((0-ef) .* angle_to_perp) = (-ef_x)*(-sin) + (-ef_y)*(cos)
  Kokkos::parallel_for(
      "compute_Ax", Nx, KOKKOS_LAMBDA(int i) {
        Real theta = x(i);
        Real perp_x = -Kokkos::sin(theta);
        Real perp_y = Kokkos::cos(theta);
        Real vec_i = (-ef_x) * perp_x + (-ef_y) * perp_y;

        Real ri = (max_val > 0.0) ? -real_interaction(i) / max_val : 0.0;
        Ax(i) = ri + vec_i;
      });
}

// Computes the ff_influence kernel value at a single point x
// shape = exp(-5*(x + pi/16)^2) - exp(-5*(x - pi/16)^2)
// return_data = shape / max(shape)
//
// Note: Normalization requires a two-pass approach in parallel:
//   Pass 1: compute raw shape values and find the max
//   Pass 2: normalize by the max

KOKKOS_INLINE_FUNCTION Real ff_influence_shape(Real x) {
  const Real shift = M_PI / 16.0;
  const Real left = x + shift;
  const Real right = x - shift;
  return Kokkos::exp(-5.0 * left * left) - Kokkos::exp(-5.0 * right * right);
}

// Computes c(i) = ff_influence(x(i)) for all i
// Matches MATLAB: shape/max(shape)
void compute_ff_influence(R1 c, R1 x) {

  const int Nx = x.extent(0);

  // --- Pass 1: compute raw shape and find maximum ---
  Real max_val = 0.0;
  Kokkos::parallel_reduce(
      "ff_influence_max", Nx,
      KOKKOS_LAMBDA(int i, Real &local_max) {
        Real val = ff_influence_shape(x(i));
        if (val > local_max)
          local_max = val;
      },
      Kokkos::Max<Real>(max_val));
  Kokkos::fence();

  // --- Pass 2: normalize and store ---
  Kokkos::parallel_for(
      "ff_influence_fill", Nx,
      KOKKOS_LAMBDA(int i) { c(i) = ff_influence_shape(x(i)) / max_val; });
  Kokkos::fence();
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    // Set up data structures
    R2 f("f", N_x, N_y); // probability density
    R1 x("x", N_x);      // x-coordinates
    R1 y("y", N_y);      // y-coordinates

    R1 k_x("k_x", N_x); // x-direction frequencies
    R1 k_y("k_y", N_y); // y-direction frequencies

    Z2 F("F", N_x, N_y);         // complex version of f
    Z2 F_hat("F_hat", N_x, N_y); // Fourier space buffer
    Z2 sg_x("sg_x", N_x, N_y);   // x semigroup multiplier
    Z2 sg_y("sg_y", N_x, N_y);   // y semigroup multiplier
    R1 Ax("Ax", N_x);            // A(x) values
    R1 c("c", N_x);              // ff_influence(x)

    // ef vector from angle
    Real ef_x = Kokkos::cos(ef_angle);
    Real ef_y = Kokkos::sin(ef_angle);

    Kokkos::parallel_for(
        "init_grid_x", N_x, KOKKOS_LAMBDA(int i) {
          x(i) = L_x_lower + i * dx;
          k_x(i) = fftfreq(i, N_x, dx);
        });
    Kokkos::parallel_for(
        "init_grid_y", N_y, KOKKOS_LAMBDA(int j) {
          y(j) = L_y_lower + j * dy;
          k_y(j) = fftfreq(j, N_y, dy);
        });
    Kokkos::fence();

    // Pass 1: fill unnormalized f
    Kokkos::parallel_for(
        "init_f", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
        KOKKOS_LAMBDA(int i, int j) {
          const Real theta = x(i);
          const Real omega = y(j);

          const Real x_side = Kokkos::exp(-10.0 * theta * theta);
          const Real y_side =
              Kokkos::exp(-10.0 * omega * omega) +
              Kokkos::exp(-10.0 * (omega - 0.2 * M_PI) * (omega - 0.2 * M_PI));
          f(i, j) = x_side * y_side;
        });
    Kokkos::fence();

    // Pass 2: normalize f so it integrates to 1
    // trapz in x for each j, then trapz in y
    R1 integral_x("integral_x", N_y);
    Kokkos::parallel_for(
        "trapz_x_norm", N_y, KOKKOS_LAMBDA(int j) {
          Real sum = 0.0;
          for (int i = 0; i < N_x - 1; i++) {
            sum += 0.5 * (f(i, j) + f(i + 1, j));
          }
          integral_x(j) = sum * dx;
        });
    Kokkos::fence();

    Real total_integral = 0.0;
    Kokkos::parallel_reduce(
        "trapz_y_norm", N_y,
        KOKKOS_LAMBDA(int j, Real &lsum) {
          if (j < N_y - 1) {
            lsum += 0.5 * (integral_x(j) + integral_x(j + 1)) * dy;
          }
        },
        total_integral);
    Kokkos::fence();

    Kokkos::parallel_for(
        "normalize_f",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
        KOKKOS_LAMBDA(int i, int j) { f(i, j) /= total_integral; });
    Kokkos::fence();

    // Initialize c = ff_influence(x)  -- fill with your actual function
    compute_ff_influence(c, x);

    // Precompute sg_x (constant if B only depends on y)
    compute_sg_x(sg_x, k_x, y);

    // KokkosFFT plan axes
    // Axis 0 = x-direction (second index in row-major = columns)
    // Axis 1 = y-direction (first index in row-major = rows)
    // NOTE: KokkosFFT axis convention matches array layout

    write_vtu<Real>("solution_0.vti", x, y, f);

    const int n_iterations = std::ceil(t_final / dt);

    for (int iteration = 1; iteration <= n_iterations; iteration++) {

      // --- Cast f (real) to complex F ---
      Kokkos::parallel_for(
          "real_to_complex",
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
          KOKKOS_LAMBDA(int i, int j) { F(i, j) = Complex(f(i, j), 0.0); });

      // =========================================================
      // STEP 1: FFT in x (axis=1), apply sg_x, IFFT in x
      // MATLAB: X_n0 = fft(U,[],2)  -- fft along columns
      // =========================================================
      KokkosFFT::fft(Exec(), F, F_hat, KokkosFFT::Normalization::none,
                     /*axis=*/0);

      // Apply x-semigroup multiplier
      Kokkos::parallel_for(
          "apply_sg_x",
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
          KOKKOS_LAMBDA(int i, int j) {
            F_hat(i, j) = F_hat(i, j) * sg_x(i, j);
          });

      // IFFT in x, normalize by 1/N_x
      KokkosFFT::ifft(Exec(), F_hat, F, KokkosFFT::Normalization::backward,
                      /*axis=*/0);

      // Take real part back to f
      Kokkos::parallel_for(
          "take_real_x",
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
          KOKKOS_LAMBDA(int i, int j) { f(i, j) = F(i, j).real(); });
      Kokkos::fence();

      // =========================================================
      // STEP 2: Compute A(x, ef, f, c) -> Ax
      // =========================================================
      compute_Ax(Ax, x, f, c, ef_x, ef_y, iteration * dt);

      // Compute sg_y from updated Ax
      compute_sg_y(sg_y, Ax, k_y);

      // =========================================================
      // STEP 3: FFT in y (axis=0), apply sg_y, IFFT in y
      // MATLAB: Y_n12 = fft(U,[],1)  -- fft along rows
      // =========================================================

      // Re-cast f to complex
      Kokkos::parallel_for(
          "real_to_complex_2",
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
          KOKKOS_LAMBDA(int i, int j) { F(i, j) = Complex(f(i, j), 0.0); });

      KokkosFFT::fft(Exec(), F, F_hat, KokkosFFT::Normalization::none,
                     /*axis=*/1);

      // Apply y-semigroup multiplier
      Kokkos::parallel_for(
          "apply_sg_y",
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
          KOKKOS_LAMBDA(int i, int j) {
            F_hat(i, j) = F_hat(i, j) * sg_y(i, j);
          });

      // IFFT in y
      KokkosFFT::ifft(Exec(), F_hat, F, KokkosFFT::Normalization::backward,
                      /*axis=*/1);

      // Take real part
      Kokkos::parallel_for(
          "take_real_y",
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
          KOKKOS_LAMBDA(int i, int j) { f(i, j) = F(i, j).real(); });
      Kokkos::fence();

      if (iteration % 1000 == 0) {
        write_vtu<Real>("solution_" + std::to_string(iteration) + ".vti", x, y,
                        f);
      }
    }
  }
  Kokkos::finalize();
}
