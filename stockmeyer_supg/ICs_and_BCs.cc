// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  scalar_value = 0.0;

  if (index == 0)
    {
      // Coordinate x = y
      double y     = point[0];
      double theta = point[1];
      double omega = point[2] * omega_max / pi - omega_max;
      double sigma = sigma_s;
      // double sigma = sigma_s * (1.0 + sine_coeff * sin(k * y));
      double f0 = (1.0 / (std::sqrt(2.0 * pi) * sigma)) *
                  std::exp(-omega * omega / (2.0 * sigma * sigma));
      f0 = f0 / (4.0 * pi * pi);
      // f0           = f0 * (2.0 + sin(theta)) / 2.0;
      f0 = f0 / (1.0 / (std::sqrt(2.0 * pi) * 0.1 * 2.0 * pi)) *
           std::exp(-(theta - pi) * (theta - pi) / (0.02 * 4.0 * pi * pi));
      scalar_value = f0;
    }
  else
    {
      scalar_value = 0.0;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE