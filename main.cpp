#include "Kokkos_Macros.hpp"
#include <KokkosFFT.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

// Parameters
constexpr int N_x = 128;         // Grid points in x direction
constexpr int N_y = 128;         // Grid points in y direction
constexpr double dt = 0.001;     // Timestep
constexpr double t_final = 30.0; // Final time

constexpr double gamma_3 = 1.0;
constexpr double gamma_4 = 1.0;
constexpr double E_x_0 = 0.0;
constexpr double E_y_0 = 0.0;
constexpr double sigma_star = 0.3;
constexpr double k = 2.0;
constexpr double sine_coefficient = 0.0;
constexpr double delta = 1.0; // Japanese bracket regularization

constexpr double L_x_lower = -M_PI; // Grid lower bound in x direction
constexpr double L_x_upper = M_PI;  // Grid upper bound in x direction
constexpr double L_y_lower = -M_PI; // Grid lower bound in y direction
constexpr double L_y_upper = M_PI;  // Grid upper bound in y direction

// Computed values
constexpr double dx = (L_x_upper - L_x_lower) / double(N_x);
constexpr double dy = (L_y_upper - L_y_lower) / double(N_y);

constexpr int n_points = N_x * N_y;
constexpr int n_cells = (N_x - 1) * (N_y - 1);

// Utility functions
KOKKOS_INLINE_FUNCTION std::array<double, 2> angle_to_vec(double theta) {
  return {Kokkos::cos(theta), Kokkos::sin(theta)};
}
KOKKOS_INLINE_FUNCTION std::array<double, 2> angle_to_perp(double theta) {
  return {-Kokkos::sin(theta), Kokkos::cos(theta)};
}
KOKKOS_INLINE_FUNCTION double japanese_bracket(double x) {
  return Kokkos::sqrt(x * x + delta * delta);
}

// VTU writer
void write_vtu(const std::string &filename, const Kokkos::View<double *> &x,
               const Kokkos::View<double *> &y,
               const Kokkos::View<double **> &f) {
  std::ofstream out(filename);
  out << std::setprecision(16);

  int _N_x = x.extent(0);
  int _N_y = y.extent(0);

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n";

  out << "<ImageData WholeExtent=\"0 " << _N_x - 1 << " 0 " << _N_y - 1
      << " 0 0\" "
      << "Origin=\"" << x(0) << " " << y(0) << " 0\" "
      << "Spacing=\"" << dx << " " << dy << " 1\">\n";

  out << "    <Piece Extent=\"0 " << _N_x - 1 << " 0 " << _N_y - 1
      << " 0 0\">\n";

  out << "      <PointData Scalars=\"f\">\n";

  out << "        <DataArray type=\"Float64\" Name=\"f\" format=\"ascii\">\n";

  for (int j = 0; j < _N_y; ++j) {
    for (int i = 0; i < _N_x; ++i) {
      out << f(i, j) << "\n";
    }
  }

  out << "        </DataArray>\n";
  out << "      </PointData>\n";

  out << "    </Piece>\n";
  out << "  </ImageData>\n";
  out << "</VTKFile>\n";
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    // Set up grid and coordinates
    Kokkos::View<double **> f("f", N_x, N_y);
    Kokkos::View<double *> x("x", N_x);
    Kokkos::View<double *> y("y", N_y);

    Kokkos::parallel_for(
        "init_grid_x", N_x,
        KOKKOS_LAMBDA(int i) { x(i) = L_x_lower + i * dx; });

    Kokkos::parallel_for(
        "init_grid_y", N_y,
        KOKKOS_LAMBDA(int j) { y(j) = L_y_lower + j * dy; });
    Kokkos::fence();

    Kokkos::parallel_for(
        "init_f", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N_x, N_y}),
        KOKKOS_LAMBDA(int i, int j) {
          using Kokkos::exp;
          using Kokkos::sqrt;

          const double theta = x(i);
          const double omega = y(j);

          const double exponential_term =
              exp(-omega * omega / (2.0 * sigma_star * sigma_star));
          const double exponential_prefactor =
              1.0 / (sqrt(2.0 * M_PI) * sigma_star);

          const double some_extra_terms_idk =
              1.0 / (4.0 * M_PI * M_PI) /
              (1.0 / (sqrt(2.0 * M_PI) * 0.1 * 2.0 * M_PI)) *
              exp(-theta * theta / (0.02 * 4.0 * M_PI * M_PI));

          f(i, j) =
              exponential_prefactor * exponential_term * some_extra_terms_idk;
        });
    Kokkos::fence();

    write_vtu("solution_0.vti", x, y, f);

    const int n_iterations = std::ceil(t_final / dt);

    for (int iteration = 0; iteration < n_iterations; iteration++) {

      if (iteration % 1000 == 0) {
        write_vtu("solution_" + std::to_string(iteration) + ".vti", x, y, f);
      }
    }
  }
  Kokkos::finalize();
}
