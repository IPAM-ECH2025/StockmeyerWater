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
void CustomPDE<dim, degree, number>::set_initial_condition(
    [[maybe_unused]] const unsigned int &index,
    [[maybe_unused]] const unsigned int &component,
    [[maybe_unused]] const dealii::Point<dim> &point,
    [[maybe_unused]] number &scalar_value,
    [[maybe_unused]] number &vector_component_value) const {

  if (index == 0) {
    if constexpr (dim == 3) {
      auto sigma = [&](auto _y) {
        return sigma_star * (1.0 + sine_coefficient * std::sin(k * _y));
      };

      auto y = point[0];
      auto theta = point[1];
      auto omega = point[2];

      auto exponential_term =
          std::exp(-omega * omega / (2.0 * sigma(y) * sigma(y)));
      auto exponential_prefactor = 1.0 / (std::sqrt(2.0 * M_PI) * sigma(y));

      auto some_extra_terms_idk =
          1.0 / (4.0 * M_PI * M_PI) /
          (1.0 / (std::sqrt(2.0 * M_PI) * 0.1 * 2.0 * M_PI)) *
          std::exp(-theta * theta / (0.02 * 4.0 * M_PI * M_PI));

      scalar_value =
          exponential_prefactor * exponential_term * some_extra_terms_idk;
    } else if constexpr (dim == 2) {

      auto theta = point[0];
      auto omega = point[1];

      auto exponential_term =
          std::exp(-omega * omega / (2.0 * sigma_star * sigma_star));
      auto exponential_prefactor = 1.0 / (std::sqrt(2.0 * M_PI) * sigma_star);

      auto some_extra_terms_idk =
          1.0 / (4.0 * M_PI * M_PI) /
          (1.0 / (std::sqrt(2.0 * M_PI) * 0.1 * 2.0 * M_PI)) *
          std::exp(-theta * theta / (0.02 * 4.0 * M_PI * M_PI));

      scalar_value =
          exponential_prefactor * exponential_term * some_extra_terms_idk;
    }
  }
}

template <unsigned int dim, unsigned int degree, typename number>
void CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
    [[maybe_unused]] const unsigned int &index,
    [[maybe_unused]] const unsigned int &boundary_id,
    [[maybe_unused]] const unsigned int &component,
    [[maybe_unused]] const dealii::Point<dim> &point,
    [[maybe_unused]] number &scalar_value,
    [[maybe_unused]] number &vector_component_value) const {}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
