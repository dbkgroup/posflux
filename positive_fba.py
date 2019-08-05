'''
(c) University of Liverpool 2019

All rights reserved.
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-locals
import lp_analyse
import lp_solver
import positive_constraint_based_modelling
import sbml_read


def main():
    '''main method.'''
    model = 'data/yeast_5.21_MCISB.xml'

    s_matrix, _, num_reactions, reaction_names, lbs, ubs, \
        _ = sbml_read.retrieve_information_from_sbml(model)

    # Force glucose into the system (note the exchange reaction is written in
    # the above model as glucose moving from inside the cell to outside, so
    # the flux is negative in the SBML).
    glucose_index = reaction_names.index('D-glucose exchange')
    lbs[glucose_index] = -16.5
    ubs[glucose_index] = -16.5

    # Set the the growth reaction as the objective
    growth_index = reaction_names.index('growth')
    objective = {}
    objective[growth_index] = -1

    lhs, equalities, rhs, variable_types, lbs, ubs, objective = \
        positive_constraint_based_modelling.build_positive_fba(
            s_matrix, lbs, ubs, objective)  # Build matrix problem

    # Build solver problem
    lp = lp_solver.create_gurobi_lp(
        lhs, equalities, rhs, variable_types, lbs, ubs, objective)

    lp_solved = lp_solver.run_gurobi_lp(lp)  # Solve LP

    flux = lp_analyse.extract_flux_patterns_from_positive_fba(
        lp_solved, num_reactions)  # Extract flux pattern

    for i, val in enumerate(flux):
        print(round(val, 3), '\t', reaction_names[i])
