'''
(c) University of Liverpool 2019

All rights reserved.
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-locals
# pylint: disable=wrong-import-order
import lp_analyse
import lp_solver
import matplotlib.pyplot as plt
import positive_constraint_based_modelling
import sbml_data_mapper


def main():
    '''main method.'''
    model = 'data/yeast_5.21_MCISB.xml'
    data = 'data/genedata_75.txt'

    s_matrix, _, num_reactions, reaction_names, lbs, ubs, \
        rxn_mean_expression, rxn_stdev_expression = \
        sbml_data_mapper.read_and_map_sbml_and_expression_data(
            model, data, 'D-glucose exchange', 16.5)

    # Update bounds (refers to the original directionality from the SBML):
    glucose_index = reaction_names.index('D-glucose exchange')
    lbs[glucose_index] = -16.5
    ubs[glucose_index] = -16.5

    # Build matrix problem:
    lhs, equalities, rhs, variable_types, lbs, ubs, objective = \
        positive_constraint_based_modelling.build_super_daaave(
            s_matrix, lbs, ubs, rxn_mean_expression, rxn_stdev_expression)

    # Build solver problem
    lp = lp_solver.create_gurobi_lp(
        lhs, equalities, rhs, variable_types, lbs, ubs, objective)

    lp_solved = lp_solver.run_gurobi_lp(lp)  # Solve LP

    flux = lp_analyse.extract_flux_patterns_from_super_daaave(
        lp_solved, num_reactions)  # Extract flux pattern

    # Plotting some results
    plt.scatter(rxn_mean_expression, flux, s=rxn_stdev_expression)
    plt.title('Transcript Abundance vs Nearest Possible Flux')
    plt.xlabel('Abundance')
    plt.ylabel('Flux')
    plt.show()
