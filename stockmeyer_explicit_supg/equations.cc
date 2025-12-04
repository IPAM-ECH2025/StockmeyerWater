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
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "f");
  set_dependencies_gradient_term_rhs(0, "f");
}

template <unsigned int dim, unsigned int degree, typename number>
void CustomPDE<dim, degree, number>::compute_explicit_rhs(
    [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
    [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>>
        &q_point_loc,
    [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
    [[maybe_unused]] Types::Index solve_block) const {
  if constexpr (dim == 3) {
    ScalarValue y = q_point_loc[0];
    ScalarValue theta = q_point_loc[1];
    ScalarValue omega = q_point_loc[2];

    ScalarValue f = variable_list.template get_value<ScalarValue>(0);

    // The components of the vector perpendicular to the normal vector
    ScalarGrad orthonormal;
    orthonormal[0] = -std::sin(theta);
    orthonormal[1] = std::cos(theta);
    orthonormal[2] = 0.0;

    // Electric field
    ScalarGrad E;
    E[0] = E_x_0 + 0.0 * std::sin(dt);
    E[1] = E_y_0 + 0.0 * std::cos(dt);
    E[2] = 0.0;

    // TODO: Add G term here
    ScalarGrad G;
    G = 0.0 * G;

    // Create a vector for the velocity
    ScalarGrad velocity;
    velocity[0] = 0.0;
    velocity[1] = -omega;
    velocity[2] = -(nu_3 * G - nu_4 * E) * orthonormal;

    // u[2] = nu4 * Edotnn * (pi / omega_max);

    variable_list.set_value_term(0, f);
    variable_list.set_gradient_term(0, -dt * velocity * f);
  }
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
