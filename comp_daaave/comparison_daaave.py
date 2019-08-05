'''
(c) University of Liverpool 2019

All rights reserved.

The ComparisonDaaave algorithm extends the principle of the original Daaaaave
algorithm that seeks flux patterns which correlate with measured
transcript/protein abundances.
Instead ComparisonDaaave algorithm seeks two flux patterns, the ratios of which
correlate with measured ratios of proteins (as from iTRAQ, for example).

The algorithm requires a model in SBML, experimental data in TSV format
(gene\tmean) - note stdevs are not yet accounted for within the objective,
and...

...the model must be forced into life by pushing flux through reactions.
The formulation of the problem is such that the 'best' solution is to select
flux1 as zero as this allows flux2 to also be zero and the ratio is ignored.
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=wrong-import-order
import math

from scipy.sparse import csr_matrix, eye, hstack, vstack

import lp_solver
import matplotlib.pyplot as plt
import numpy as np
import positive_constraint_based_modelling
import sbml_data_mapper


def build_comparison_daaave(s_matrix, lbs1_in, ubs1_in, lbs2_in, ubs2_in,
                            target, default_bound=10000):
    '''Build Comparision Daaaave.'''
    # Make bounds for positive formulation
    print('\tRecast bounds for positive fluxes')
    lbs1, ubs1 = positive_constraint_based_modelling.make_fluxes_positive(
        lbs1_in, ubs1_in)
    lbs2, ubs2 = positive_constraint_based_modelling.make_fluxes_positive(
        lbs2_in, ubs2_in)

    # The easiest way to write this formulation is to concatenate S and -S to
    # double the variables. The presolve deals with this very rapidly, so it's
    # not worth the hassle of the extra code
    num_species, num_reactions = s_matrix.get_shape()
    print('\tS matrix is', num_species, 'x', num_reactions)

    # Build lhs
    print('\tBuilding lhs')

    print('\tBuild Positive FBA block for S1')
    zeros_species_by_reactions = csr_matrix((num_species, num_reactions))
    zeros_reactions_by_reactions = csr_matrix((num_reactions, num_reactions))
    ones_reactions_by_reactions = eye(num_reactions).tocsr()
    neg_ones_reactions_by_reactions = -1 * ones_reactions_by_reactions

    # note that this should change if limits change
    thousands_reactions_by_reactions = default_bound * \
        ones_reactions_by_reactions

    p_stack1 = hstack([s_matrix,
                       -1 * s_matrix,
                       zeros_species_by_reactions,
                       zeros_species_by_reactions,
                       zeros_species_by_reactions])
    p_stack2 = hstack([neg_ones_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       thousands_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions])
    p_stack3 = hstack([zeros_reactions_by_reactions,
                       neg_ones_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       thousands_reactions_by_reactions,
                       zeros_reactions_by_reactions])
    p_stack4 = hstack([zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       ones_reactions_by_reactions,
                       ones_reactions_by_reactions,
                       zeros_reactions_by_reactions])
    p_stack5 = hstack([ones_reactions_by_reactions,
                       ones_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       neg_ones_reactions_by_reactions])

    pos_fba_block1 = vstack([p_stack1, p_stack2, p_stack3, p_stack4, p_stack5])

    print('\tBuild Positive FBA block for S2')
    # Seems pointless but we might need this later
    pos_fba_block2 = pos_fba_block1

    print('\tBuild MOMA block')
    target_reactions_by_reactions = csr_matrix(np.diagflat(target))
    neg_target_reactions_by_reactions = -1 * target_reactions_by_reactions

    # Build stacks
    print('\t\tBuilding stacks')
    m_stack1 = hstack([target_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       neg_ones_reactions_by_reactions,
                       neg_ones_reactions_by_reactions,
                       zeros_reactions_by_reactions])

    m_stack2 = hstack([neg_target_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       ones_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       neg_ones_reactions_by_reactions])

    # Vertical stack
    moma_block = vstack([m_stack1, m_stack2])

    print('\tCombine blocks')
    zeros_upper_right = csr_matrix(
        (num_species + (4 * num_reactions), 7 * num_reactions))
    zeros_lower_left = csr_matrix((2 * num_reactions, 4 * num_reactions))
    zeros_middle_left = csr_matrix(
        (num_species + (4 * num_reactions), 5 * num_reactions))
    zeros_middle_right = csr_matrix(
        (num_species + (4 * num_reactions), 2 * num_reactions))

    upper_stack = hstack([pos_fba_block1, zeros_upper_right])
    middle_stack = hstack(
        [zeros_middle_left, pos_fba_block2, zeros_middle_right])
    lower_stack = hstack([zeros_lower_left, moma_block])

    # Process lower stack to zero out rows associated with -1 flagged reactions
    print('\t\tEmptying lower stack of data-less constraints')
    lower_stack = lower_stack.tocsr()

    count = 0
    for row_idx, row in enumerate(target):
        if math.isnan(row):
            # Clear constraints on excesses
            lower_stack.data[lower_stack.indptr[row_idx]:
                             lower_stack.indptr[row_idx + 1]] = float(0)
            lower_stack.data[lower_stack.indptr[row_idx + num_reactions]:
                             lower_stack.indptr[row_idx + 1 + num_reactions]] \
                = float(0)  # Clear constraints on deficits

            count += 1

    print('\t\t\t', count * 2, 'data-less MoMA constraints removed')

    lhs = vstack([upper_stack, middle_stack, lower_stack])

    # Equalities
    print('\tBuilding equalities')
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_eqs = ['='] * num_species
    stack2_eqs = ['>'] * num_reactions
    stack3_eqs = ['>'] * num_reactions
    stack4_eqs = ['<'] * num_reactions
    stack5_eqs = ['='] * num_reactions
    # Stacks 6-10 are associated with the positive FBA formulation for S2
    stack6_eqs = ['='] * num_species
    stack7_eqs = ['>'] * num_reactions
    stack8_eqs = ['>'] * num_reactions
    stack9_eqs = ['<'] * num_reactions
    stack10_eqs = ['='] * num_reactions
    # Stacks 11-12 are associated with the MoMA formulation
    stack11_eqs = ['<'] * num_reactions
    stack12_eqs = ['<'] * num_reactions

    equalities = stack1_eqs + stack2_eqs + stack3_eqs + stack4_eqs + \
        stack5_eqs + stack6_eqs + stack7_eqs + stack8_eqs + stack9_eqs + \
        stack10_eqs + stack11_eqs + stack12_eqs

    # rhs
    print('\tBuild rhs')
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_rhs = [0.0] * num_species
    stack2_rhs = [0.0] * num_reactions
    stack3_rhs = [0.0] * num_reactions
    stack4_rhs = [1.0] * num_reactions
    stack5_rhs = [0.0] * num_reactions
    # Stacks 6-10 are associated with the positive FBA formulation for S2
    stack6_rhs = [0.0] * num_species
    stack7_rhs = [0.0] * num_reactions
    stack8_rhs = [0.0] * num_reactions
    stack9_rhs = [1.0] * num_reactions
    stack10_rhs = [0.0] * num_reactions
    # Stacks 11-12 are associated with the MoMA formulation
    stack11_rhs = [0.0] * num_reactions
    stack12_rhs = [0.0] * num_reactions

    rhs = stack1_rhs + stack2_rhs + stack3_rhs + stack4_rhs + stack5_rhs + \
        stack6_rhs + stack7_rhs + stack8_rhs + stack9_rhs + stack10_rhs + \
        stack11_rhs + stack12_rhs

    # _bounds
    print('\tBuild bounds')
    binl_bounds = [0] * (2 * num_reactions)
    binu_bounds = [1] * (2 * num_reactions)
    contl_bounds = [0.0] * num_reactions
    contu_bounds = [default_bound] * num_reactions

    # _bounds associated with Pos FBA blocks
    #             a+b    c+d          e
    fba1_lbs = lbs1 + binl_bounds + contl_bounds
    fba1_ubs = ubs1 + binu_bounds + contu_bounds
    fba2_lbs = lbs2 + binl_bounds + contl_bounds
    fba2_ubs = ubs2 + binu_bounds + contu_bounds

    # _bounds associated with 1-norm minimisation
    #           x             d
    moma_lbs = contl_bounds + contl_bounds
    moma_ubs = contu_bounds + contu_bounds

    lbs = fba1_lbs + fba2_lbs + moma_lbs
    ubs = fba1_ubs + fba2_ubs + moma_ubs

    # Type variables
    print('\tBuild variable types')
    a = ['continuous'] * num_reactions
    b = ['continuous'] * num_reactions
    c = ['binary'] * num_reactions
    # d = ['binary'] * num_reactions
    e = ['continuous'] * num_reactions
    x = ['continuous'] * num_reactions
    d = ['continuous'] * num_reactions
    variable_types = a + b + c + d + e + a + b + c + d + e + x + d

    # Objectives
    objective = {}

    for v in range(11 * num_reactions, 12 * num_reactions):
        objective[v] = 1
        objective[v + num_reactions] = 1

    print('\tBuilt...')
    print('\t\tlhs:', lhs.get_shape())
    print('\t\tVar:', len(variable_types))
    print('\t\tLBs:', len(lbs))
    print('\t\tUBs:', len(ubs))
    print('\t\tEqs:', len(equalities))
    print('\t\trhs:', len(rhs))

    return lhs, equalities, rhs, variable_types, lbs, ubs, objective


def _extract_flux_patterns(lp_soll, num_reactions):
    var_values = lp_soll.getVars()
    flux_vector1 = []

    for tic in range(4 * num_reactions, 5 * num_reactions):
        flux_obj = var_values[tic]
        # Tidy up some weird sign errors from Gurobi (minus zero?)
        f_squared = math.pow(flux_obj.X, 2)
        f_sqrt = round(math.sqrt(f_squared), 4)
        flux_vector1.append(f_sqrt)

    flux_vector2 = []
    for tic in range(9 * num_reactions, 10 * num_reactions):
        flux_obj = var_values[tic]
        # Tidy up some weird sign errors from Gurobi (minus zero?)
        f_squared = math.pow(flux_obj.X, 2)
        f_sqrt = round(math.sqrt(f_squared), 4)
        flux_vector2.append(f_sqrt)

    return flux_vector1, flux_vector2


def _plot(rxn_mean_ratios, flux1, flux2):
    '''Plot.'''
    x = []
    y = []

    for idx, rxn_mean_ratio in enumerate(rxn_mean_ratios):
        if rxn_mean_ratio > -1:
            if flux1[idx] and flux2[idx]:
                x.append(rxn_mean_ratio)
                y.append(flux2[idx] / flux1[idx])

    # Plotting some results
    plt.scatter(x, y)
    plt.title('Protein Ratio vs Nearest Possible Flux Ratio')
    plt.xlabel('log(Protein Ratio)')
    plt.ylabel('log(Flux Ratio)')
    plt.show()


def main():
    '''main method.'''
    model = 'data/yeast_5.21_MCISB.xml'
    data = 'data/yeastComparison.txt'

    print('- Fetch model and map in relative expression data')
    s_matrix, _, num_reactions, reaction_names, lbs, ubs, \
        rxn_mean_ratios \
        = sbml_data_mapper.read_and_map_sbml_and_relative_data(model, data)

    # Force flux through the system:
    lbs1 = lbs
    ubs1 = ubs
    lbs2 = lbs
    ubs2 = ubs

    glucose_index = reaction_names.index('D-glucose exchange')
    glucose_flux_value1 = -16.5
    glucose_flux_value2 = -5.00
    lbs1[glucose_index] = glucose_flux_value1
    ubs1[glucose_index] = glucose_flux_value1
    lbs2[glucose_index] = glucose_flux_value2
    ubs2[glucose_index] = glucose_flux_value2

    growth_index = reaction_names.index('growth')
    growth_flux_value1 = 1.00
    growth_flux_value2 = 0.50
    lbs1[growth_index] = growth_flux_value1
    ubs1[growth_index] = growth_flux_value1
    lbs2[growth_index] = growth_flux_value2
    ubs2[growth_index] = growth_flux_value2

    # Build ComparisonDaaave
    print('- Build ComparisonDaaave problem')

    # Build matrix problem
    lhs, equalities, rhs, variable_types, lbs, ubs, objective = \
        build_comparison_daaave(s_matrix, lbs1, ubs1, lbs2, ubs2,
                                rxn_mean_ratios)

    print('- Create Gurobi instance of problem')

    # Build solver problem
    lp = lp_solver.create_gurobi_lp(
        lhs, equalities, rhs, variable_types, lbs, ubs, objective)

    print('- Solve ComparisonDaaave')
    lp_solved = lp_solver.run_gurobi_lp(lp)  # Solve LP

    flux1, flux2 = _extract_flux_patterns(
        lp_solved, num_reactions)  # Extract flux patterns for models 1 and 2

    _plot(rxn_mean_ratios, flux1, flux2)
