'''
(c) University of Liverpool 2019

All rights reserved.
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
import math

from scipy import sparse
from scipy.sparse import csr_matrix, hstack, vstack


def build_positive_fba(s_matrix, lbs_in, ubs_in, objective,
                       default_bound=10000):
    '''
    S			is a CSR matrix containing the stoichiometries
    lbs			is a list of lower bounds for each reaction
    ubs			is a list of upper bounds for each reaction
    objective 	is a dictionary of reaction indices mapping to weights
                (all absent values are weighted=0)
    '''

    # Make bounds for positive formulation
    lbs, ubs = make_fluxes_positive(lbs_in, ubs_in)

    num_species, num_reactions = s_matrix.get_shape()

    print('\tBuild Positive FBA block for S1')
    zeros_species_by_reactions = csr_matrix((num_species, num_reactions))
    zeros_reactions_by_reactions = csr_matrix((num_reactions, num_reactions))
    ones_reactions_by_reactions = sparse.eye(num_reactions).tocsr()
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
    lhs = vstack([p_stack1, p_stack2, p_stack3, p_stack4, p_stack5])

    # Equalities
    print('\tBuilding equalities')
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_eqs = ['='] * num_species
    stack2_eqs = ['>'] * num_reactions
    stack3_eqs = ['>'] * num_reactions
    stack4_eqs = ['<'] * num_reactions
    stack5_eqs = ['='] * num_reactions
    equalities = stack1_eqs + stack2_eqs + stack3_eqs + stack4_eqs + stack5_eqs

    # rhs
    print('\tBuild rhs')
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_rhs = [0.0] * num_species
    stack2_rhs = [0.0] * num_reactions
    stack3_rhs = [0.0] * num_reactions
    stack4_rhs = [1.0] * num_reactions
    stack5_rhs = [0.0] * num_reactions
    rhs = stack1_rhs + stack2_rhs + stack3_rhs + stack4_rhs + stack5_rhs

    # _bounds
    print('\tBuild bounds')
    binl_bounds = [0] * (2 * num_reactions)
    binu_bounds = [1] * (2 * num_reactions)
    contl_bounds = [0.0] * num_reactions
    contu_bounds = [default_bound] * num_reactions

    # 	  a+b   c+d          e
    lbs = lbs + binl_bounds + contl_bounds
    ubs = ubs + binu_bounds + contu_bounds

    # Type variables
    print('\tBuild variable types')
    a = ['continuous'] * num_reactions
    b = ['continuous'] * num_reactions
    c = ['binary'] * num_reactions
    d = ['binary'] * num_reactions
    e = ['continuous'] * num_reactions
    variable_types = a + b + c + d + e

    # Objectives
    new_objective = {}

    for key in objective:
        val = objective[key]
        new_key = (4 * num_reactions) + key
        new_objective[new_key] = val

    return(lhs, equalities, rhs, variable_types, lbs, ubs, new_objective)


def build_super_daaave(s_matrix, lbs_in, ubs_in, target, stdevs,
                       default_bound=10000):
    '''
    SuperDaaave is an algorithm that seeks to correlate measured abundance and
    estimated flux.
    This is the Positive version of the algorithm that runs more quickly
    because all fluxes must be positive.
    '''

    # Make bounds for positive formulation
    print('\tReconfiguring bounds')
    lbs, ubs = make_fluxes_positive(lbs_in, ubs_in)

    # The easiest way to write this formulation is to concatenate S and -S to
    # double the variables. The presolve deals with this very rapidly, so it's
    # not worth the hassle of the extra code
    num_species, num_reactions = s_matrix.get_shape()
    print('\tS matrix is', num_species, 'x', num_reactions)

    # Build lhs
    print('\tBuilding lhs')

    print('\tBuild Positive FBA block')
    zeros_species_by_reactions = csr_matrix((num_species, num_reactions))
    zeros_reactions_by_reactions = csr_matrix((num_reactions, num_reactions))
    ones_reactions_by_reactions = sparse.eye(num_reactions).tocsr()
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
    pos_fba_block = vstack([p_stack1, p_stack2, p_stack3, p_stack4, p_stack5])

    print('\tBuild MOMA block')
    m_stack1 = hstack([ones_reactions_by_reactions,
                       neg_ones_reactions_by_reactions,
                       neg_ones_reactions_by_reactions,
                       zeros_reactions_by_reactions])
    m_stack2 = hstack([neg_ones_reactions_by_reactions,
                       ones_reactions_by_reactions,
                       zeros_reactions_by_reactions,
                       neg_ones_reactions_by_reactions])

    moma_block = vstack([m_stack1, m_stack2])

    print('\tCombine blocks')
    zeros_upper_right = csr_matrix(
        (num_species + (4 * num_reactions), 3 * num_reactions))
    zeros_lower_left = csr_matrix((2 * num_reactions, 4 * num_reactions))
    upper_stack = hstack([pos_fba_block, zeros_upper_right])
    lower_stack = hstack([zeros_lower_left, moma_block])
    lhs = vstack([upper_stack, lower_stack])

    # Equalities
    print('\tBuilding equalities')
    stack1_eqs = ['='] * num_species
    stack2_eqs = ['>'] * num_reactions
    stack3_eqs = ['>'] * num_reactions
    stack4_eqs = ['<'] * num_reactions
    stack5_eqs = ['='] * num_reactions
    stack6_eqs = ['<'] * num_reactions
    stack7_eqs = ['<'] * num_reactions
    equalities = stack1_eqs + stack2_eqs + stack3_eqs + \
        stack4_eqs + stack5_eqs + stack6_eqs + stack7_eqs

    # rhs
    print('\tBuild rhs')
    stack1_rhs = [0.0] * num_species
    stack2_rhs = [0.0] * num_reactions
    stack3_rhs = [0.0] * num_reactions
    stack4_rhs = [1.0] * num_reactions
    stack5_rhs = [0.0] * num_reactions
    stack6_rhs = [0.0] * num_reactions
    stack7_rhs = [0.0] * num_reactions
    rhs = stack1_rhs + stack2_rhs + stack3_rhs + \
        stack4_rhs + stack5_rhs + stack6_rhs + stack7_rhs

    # _bounds
    print('\tBuild bounds')
    binl_bounds = [0] * (2 * num_reactions)
    binu_bounds = [1] * (2 * num_reactions)
    contl_bounds = [0.0] * num_reactions
    contu_bounds = [default_bound] * num_reactions

    # Process target bounds to remove NaNs. Matrix block structure is
    # maintained, meaning NaN variables are completely free and some
    # constraints are inactive. These are left to the presolver to deal with.
    l_target = []
    u_target = []
    nan_count = 0

    for t in target:
        if math.isnan(t):  # Variable can take any value
            l_target.append(0)
            u_target.append(default_bound)
            nan_count += 1
        else:  # Variable is constrained to abundance data
            l_target.append(t)
            u_target.append(t)

    # 	  a+b   c+d          e             t         x             d
    lbs = lbs + binl_bounds + contl_bounds + l_target + contl_bounds + \
        contl_bounds
    ubs = ubs + binu_bounds + contu_bounds + u_target + contu_bounds + \
        contu_bounds

    # Type variables
    print('\tBuild variable types')
    a = ['continuous'] * num_reactions
    b = ['continuous'] * num_reactions
    c = ['binary'] * num_reactions
    # d = ['binary'] * num_reactions
    e = ['continuous'] * num_reactions
    t = ['continuous'] * num_reactions
    x = ['continuous'] * num_reactions
    d = ['continuous'] * num_reactions
    variable_types = a + b + c + d + e + t + x + d

    # Objectives
    objective = {}
    minimum_stdev = min(x for x in stdevs if x > 0)
    # Excesses and deficits (this is where 1/stdev is implemented)
    fetch = 0
    for v in range(6 * num_reactions, 7 * num_reactions):
        if not math.isnan(stdevs[fetch]):
            std = stdevs[fetch]
            # This accounts for STDEV and zeros
            obj = 1 / (std + (minimum_stdev / 2))
            objective[v] = obj
            objective[v + num_reactions] = obj
        fetch += 1

    print('\tBuilt...')
    print('\t\tlhs:', lhs.get_shape())
    print('\t\tVar:', len(variable_types), '=>', nan_count, 'NaNs opened')
    print('\t\tLBs:', len(lbs))
    print('\t\tUBs:', len(ubs))
    print('\t\tEqs:', len(equalities))
    print('\t\trhs:', len(rhs))
    print('\t\tObj:', len(objective))

    return(lhs, equalities, rhs, variable_types, lbs, ubs, objective)


def make_fluxes_positive(lbs, ubs):
    '''Converts lists of length num_reactions describing lower and upper bounds
    into lists of length 2*num_reactions describing lower and upper bounds
    associated with concatenated matrices S, -S.

    This preserves reaction directionality as written in the SBML
    (from lbs, ubs), but recasts all fluxes as zero to +ub
    '''
    pos_u_bounds = []
    pos_l_bounds = []
    neg_u_bounds = []
    neg_l_bounds = []
    typer = {}
    typer['Reversible'] = 0
    typer['Zero or Forwards'] = 0
    typer['Zero or Backwards'] = 0
    typer['Compulsory Backwards'] = 0
    typer['Compulsory Forwards'] = 0
    typer['Dead'] = 0
    typer['Other'] = 0

    for lb, ub in zip(lbs, ubs):
        if lb < 0 < ub:
            typer['Reversible'] += 1
            pos_u_bounds.append(ub)  # Alive
            pos_l_bounds.append(0.0)  # Alive
            neg_u_bounds.append(-1 * lb)  # Alive
            neg_l_bounds.append(0.0)  # Alive
        elif lb == 0 and ub > 0:
            typer['Zero or Forwards'] += 1
            pos_u_bounds.append(ub)  # Alive
            pos_l_bounds.append(0.0)  # Alive
            neg_u_bounds.append(0.0)  # Dead
            neg_l_bounds.append(0.0)  # Dead
        elif lb < 0 and ub == 0:
            typer['Zero or Backwards'] += 1
            pos_u_bounds.append(0.0)  # Dead
            pos_l_bounds.append(0.0)  # Dead
            neg_u_bounds.append(-1 * lb)  # Alive
            neg_l_bounds.append(0.0)  # Alive
        elif lb == 0 and ub == 0:
            typer['Dead'] += 1
            pos_u_bounds.append(0.0)  # Dead
            pos_l_bounds.append(0.0)  # Dead
            neg_u_bounds.append(0.0)  # Dead
            neg_l_bounds.append(0.0)  # Dead
        elif lb > 0 and ub > 0:
            typer['Compulsory Forwards'] += 1
            pos_u_bounds.append(ub)  # Alive
            pos_l_bounds.append(lb)  # Alive
            neg_u_bounds.append(0.0)  # Dead
            neg_l_bounds.append(0.0)  # Dead
        elif lb < 0 and ub < 0:
            typer['Compulsory Backwards'] += 1
            pos_u_bounds.append(0.0)  # Dead
            pos_l_bounds.append(0.0)  # Dead
            neg_u_bounds.append(-1 * lb)  # Alive
            neg_l_bounds.append(-1 * ub)  # Alive
        else:
            print('What!? THIS SHOULD NOT HAPPEN!')
            typer['Other'] += 1
            print('Other', lb, ub)

    sm = 0

    for key in typer:
        print('\t\t' + str(typer[key]) + '\t' + key)
        sm = sm + typer[key]

    u_bounds = pos_u_bounds + neg_u_bounds
    l_bounds = pos_l_bounds + neg_l_bounds

    return(l_bounds, u_bounds)
