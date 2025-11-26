// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  set_variable_name(0, "f");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ImplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "f, old_1(f), grad(f)");
  set_dependencies_gradient_term_rhs(0, "f, old_1(f), grad(f)");
  set_dependencies_value_term_lhs(0, "change(f), grad(change(f))");
  set_dependencies_gradient_term_lhs(0, "change(f), grad(change(f))");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{
  // Direction indices
  //  0 = y direction
  ScalarValue X = q_point_loc[0];
  // 1 = theta direction
  ScalarValue Y     = q_point_loc[1];
  ScalarValue theta = Y;
  // 2 = omega_r direction (p[2] = pi*(omega+omega_max)/omega_max)
  ScalarValue Z     = q_point_loc[2];
  ScalarValue omega = (omega_max / pi) * Z - omega_max;

  // Obtaining f, its gradient and previous time step value
  ScalarValue f     = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  fX    = variable_list.template get_gradient<ScalarGrad>(0);
  ScalarValue f_old = variable_list.template get_value<ScalarValue>(0, OldOne);

  // The components of the vector perpendiclar to the normal vector
  // nn(\theta)=n(\theta+\pi/2)
  ScalarValue nnx = -std::sin(theta);
  ScalarValue nny = std::cos(theta);

  // Time-dependent componenents of the electric field
  ScalarValue Ex = Ex0 + 0.0 * std::sin(get_timestep());
  ScalarValue Ey = Ey0 + 0.0 * std::cos(get_timestep());

  //"Velocity term"
  ScalarGrad u;
  u[0] = 0.0;
  u[1] = omega;
  u[2] = 0.0;

  // Stabilization parameter
  ScalarValue stabilization_parameter =
    compute_stabilization_parameter_2(u, element_volume);

  // Advection term u dot grad(f)
  ScalarValue advection_term = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      advection_term = advection_term + u[i] * fX[i];
    }

  ScalarValue residual = -(f - f_old) - get_timestep() * advection_term;
  ScalarValue eq_f     = residual;
  ScalarGrad  eqx_f    = residual * stabilization_parameter * u;
  variable_list.set_value_term(0, eq_f);
  variable_list.set_gradient_term(0, eqx_f);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{
  // Direction indices
  // 0 = y direction
  ScalarValue X = q_point_loc[0];
  // 1 = theta direction
  ScalarValue Y     = q_point_loc[1];
  ScalarValue theta = Y;
  // 2 = omega_r direction (p[2] = pi*(omega+omega_max)/omega_max)
  ScalarValue Z     = q_point_loc[2];
  ScalarValue omega = (omega_max / pi) * Z - omega_max;

  // Obtaining change_f and its gradient
  ScalarValue change_f  = variable_list.template get_value<ScalarValue>(0, Change);
  ScalarGrad  change_fX = variable_list.template get_gradient<ScalarGrad>(0, Change);

  // The components of the vector perpendiclar to the normal vector
  // nn(\theta)=n(\theta+\pi/2)
  ScalarValue nnx = -std::sin(theta);
  ScalarValue nny = std::cos(theta);

  // Time-dependent componenents of the electric field
  ScalarValue Ex = Ex0 + 0.0 * std::sin(get_timestep());
  ScalarValue Ey = Ey0 + 0.0 * std::cos(get_timestep());

  //"Velocity term"
  ScalarGrad u;
  u[0] = 0.0;
  u[1] = omega;
  u[2] = 0.0;

  // Stabilization parameter
  ScalarValue stabilization_parameter =
    compute_stabilization_parameter_2(u, element_volume);

  // Advection term u dot grad(f)
  ScalarValue advection_term = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      advection_term = advection_term + u[i] * change_fX[i];
    }

  ScalarValue residual = change_f + get_timestep() * advection_term;
  ScalarValue eq_f     = residual;
  ScalarGrad  eqx_f    = residual * stabilization_parameter * u;
  variable_list.set_value_term(0, eq_f, Change);
  variable_list.set_gradient_term(0, eqx_f, Change);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE