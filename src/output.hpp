#pragma once

#include <Kokkos_Core.hpp>
#include <fstream>
#include <iomanip>
#include <string>

template <typename RealType>
void write_vtu(const std::string &filename, const Kokkos::View<RealType *> &x,
               const Kokkos::View<RealType *> &y,
               const Kokkos::View<RealType **> &f) {
  std::ofstream out(filename);
  out << std::setprecision(16);

  auto _N_x = x.extent(0);
  auto _N_y = y.extent(0);
  auto _dx = x(1) - x(0);
  auto _dy = y(1) - y(0);

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n";

  out << "<ImageData WholeExtent=\"0 " << _N_x - 1 << " 0 " << _N_y - 1
      << " 0 0\" "
      << "Origin=\"" << x(0) << " " << y(0) << " 0\" "
      << "Spacing=\"" << _dx << " " << _dy << " 1\">\n";

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
