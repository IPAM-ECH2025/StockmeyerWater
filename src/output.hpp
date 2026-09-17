#pragma once

#include "parameters.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

template <typename RealType, typename ViewType3D, typename ViewType1D>
void write_vtu(const ViewType3D &f, const ViewType1D &x, const ViewType1D &y,
               const ViewType1D &z, const std::string &filename,
               RealType time = 0.0, int increment = 0,
               const std::string &output_dir = "outputs/") {
  // NOTE: This only works for 3D views
  // NOTE: Casts to float32
  using P = Parameters<RealType>;

  std::filesystem::create_directory(output_dir);
  std::string filepath = output_dir + filename;
  std::ofstream fstream(filepath);
  if (!fstream) {
    throw std::runtime_error("Failed to open file: " + filepath);
  }
  fstream << std::setprecision(std::numeric_limits<float>::max_digits10);

  const auto _N_x = P::N_x;
  const auto _N_y = P::N_y;
  const auto _N_z = P::N_z;

  fstream << "<?xml version=\"1.0\"?>\n";
  fstream << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
             "byte_order=\"LittleEndian\">\n";

  fstream << "  <RectilinearGrid WholeExtent=\""
          << "0 " << _N_x - 1 << " "
          << "0 " << _N_y - 1 << " "
          << "0 " << _N_z - 1 << "\">\n";

  fstream << "    <FieldData>\n";
  fstream << "      <DataArray type=\"Float64\" Name=\"Time\" "
             "NumberOfTuples=\"1\" format=\"ascii\">\n";
  fstream << "        " << time << "\n";
  fstream << "      </DataArray>\n";
  fstream << "      <DataArray type=\"Int32\" Name=\"Increment\" "
             "NumberOfTuples=\"1\" format=\"ascii\">\n";
  fstream << "        " << increment << "\n";
  fstream << "      </DataArray>\n";
  fstream << "    </FieldData>\n";

  fstream << "    <Piece Extent=\""
          << "0 " << _N_x - 1 << " "
          << "0 " << _N_y - 1 << " "
          << "0 " << _N_z - 1 << "\">\n";

  fstream << "      <PointData Scalars=\"f\">\n";
  fstream << "        <DataArray type=\"Float32\" Name=\"f\" "
             "NumberOfComponents=\"1\" format=\"ascii\">\n";
  for (int k = 0; k < _N_z; ++k) {
    for (int j = 0; j < _N_y; ++j) {
      for (int i = 0; i < _N_x; ++i) {
        fstream << "          " << (float)f(i, j, k) << "\n";
      }
    }
  }
  fstream << "        </DataArray>\n";
  fstream << "      </PointData>\n";

  fstream << "      <Coordinates>\n";
  fstream << "        <DataArray type=\"Float32\" Name=\"XCoordinates\" "
             "format=\"ascii\">\n";
  for (int i = 0; i < _N_x; ++i) {
    fstream << "          " << (float)x(i) << "\n";
  }
  fstream << "        </DataArray>\n";
  fstream << "        <DataArray type=\"Float32\" Name=\"YCoordinates\" "
             "format=\"ascii\">\n";
  for (int j = 0; j < _N_y; ++j) {
    fstream << "          " << (float)y(j) << "\n";
  }
  fstream << "        </DataArray>\n";
  fstream << "        <DataArray type=\"Float32\" Name=\"ZCoordinates\" "
             "format=\"ascii\">\n";
  for (int k = 0; k < _N_z; ++k) {
    fstream << "          " << (float)z(k) << "\n";
  }
  fstream << "        </DataArray>\n";
  fstream << "      </Coordinates>\n";

  fstream << "    </Piece>\n";
  fstream << "  </RectilinearGrid>\n";
  fstream << "</VTKFile>\n";
}
