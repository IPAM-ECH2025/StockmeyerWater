// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void CustomAttributeLoader::load_variable_attributes() {
  set_variable_name(0, "f");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "f");
  set_dependencies_gradient_term_rhs(0, "grad(f)");
}

template <unsigned int dim, unsigned int degree, typename number>
void CustomPDE<dim, degree, number>::compute_explicit_rhs(
    [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
    [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>>
        &q_point_loc,
    [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
    [[maybe_unused]] Types::Index solve_block) const {
  // Direction indices
  //  0 = y direction
  ScalarValue y = q_point_loc[0];
  // 1 = theta direction
  ScalarValue theta = q_point_loc[1]-pi;
  // 2 = omega_r direction (p[2] = pi*(omega+omega_max)/omega_max)
  ScalarValue omega = (omega_max / pi) * q_point_loc[2] - omega_max;

  // Obtaining f and its gradient with respect to y, theta and
  ScalarValue f = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad fx = variable_list.template get_gradient<ScalarGrad>(0);

  // The components of the vector perpendiclar to the normal vector
  // nn(\theta)=n(\theta+\pi/2)
  ScalarValue nnx = -std::sin(theta);
  ScalarValue nny = std::cos(theta);

  // Time-dependent componenents of the electric field
  ScalarValue Ex = Ex0 + 0.0 * std::sin(get_timestep());
  ScalarValue Ey = Ey0 + 0.0 * std::cos(get_timestep());

  // E dot nn
  ScalarValue Edotnn = Ex * nnx + Ey * nny;

  // Add calculation of G term here, add it to E term

  ScalarGrad u;
  u[0] = 0.0;
  u[1] = -get_timestep() * omega;
  u[2] = get_timestep() * nu4 * Edotnn * (pi / omega_max);

  variable_list.set_value_term(0, f);
  variable_list.set_gradient_term(0, -f * u);
}

template <unsigned int dim, unsigned int degree, typename number>
void CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
    [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
    [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>>
        &q_point_loc,
    [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
    [[maybe_unused]] Types::Index solve_block,
    [[maybe_unused]] Types::Index current_index) const {}

template <unsigned int dim, unsigned int degree, typename number>
void CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
    [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
    [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>>
        &q_point_loc,
    [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
    [[maybe_unused]] Types::Index solve_block,
    [[maybe_unused]] Types::Index current_index) const {}

template <unsigned int dim, unsigned int degree, typename number>
void CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
    [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
    [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>>
        &q_point_loc,
    [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
    [[maybe_unused]] Types::Index solve_block) const {}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE