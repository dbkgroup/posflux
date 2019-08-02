import math

from scipy import sparse
from scipy.sparse import dok_matrix, csr_matrix, hstack, vstack

import gurobipy
import numpy as np


def _buildPositiveFBA(S, lbsIn, ubsIn, objective, default_bound=1000000):
    '''
    S			is a CSR matrix containing the stoichiometries
    lbs			is a list of lower bounds for each reaction
    ubs			is a list of upper bounds for each reaction
    objective 	is a dictionary of reaction indices mapping to weights (all absent values are weighted=0)
    '''

    # Make bounds for positive formulation
    lbs, ubs = _makeFluxesPositive(lbsIn, ubsIn)

    sDimensions = S.get_shape()
    num_species, num_reactions = sDimensions

    print '\tBuild Positive FBA block for S1'
    zeros_species_by_reactions = csr_matrix((num_species, num_reactions))
    zeros_reactions_by_reactions = csr_matrix((num_reactions, num_reactions))
    ones_reactions_by_reactions = sparse.eye(num_reactions).tocsr()
    negOnes_reactions_by_reactions = -1 * ones_reactions_by_reactions
    thousands_reactions_by_reactions = default_bound * \
        ones_reactions_by_reactions  # note that this should change if limits change

    Pstack1 = hstack([S, -1 * S, zeros_species_by_reactions,
                      zeros_species_by_reactions, zeros_species_by_reactions])
    Pstack2 = hstack([negOnes_reactions_by_reactions, zeros_reactions_by_reactions,
                      thousands_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack3 = hstack([zeros_reactions_by_reactions, negOnes_reactions_by_reactions,
                      zeros_reactions_by_reactions, thousands_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack4 = hstack([zeros_reactions_by_reactions, zeros_reactions_by_reactions,
                      ones_reactions_by_reactions, ones_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack5 = hstack([ones_reactions_by_reactions, ones_reactions_by_reactions,
                      zeros_reactions_by_reactions, zeros_reactions_by_reactions, negOnes_reactions_by_reactions])
    LHS = vstack([Pstack1, Pstack2, Pstack3, Pstack4, Pstack5])

    # Equalities
    print '\tBuilding equalities'
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_eqs = ['='] * num_species
    stack2_eqs = ['>'] * num_reactions
    stack3_eqs = ['>'] * num_reactions
    stack4_eqs = ['<'] * num_reactions
    stack5_eqs = ['='] * num_reactions
    equalities = stack1_eqs + stack2_eqs + stack3_eqs + stack4_eqs + stack5_eqs

    # RHS
    print '\tBuild RHS'
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_RHS = [0.0] * num_species
    stack2_RHS = [0.0] * num_reactions
    stack3_RHS = [0.0] * num_reactions
    stack4_RHS = [1.0] * num_reactions
    stack5_RHS = [0.0] * num_reactions
    RHS = stack1_RHS + stack2_RHS + stack3_RHS + stack4_RHS + stack5_RHS

    # Bounds
    print '\tBuild bounds'
    binlBounds = [0] * (2 * num_reactions)
    binuBounds = [1] * (2 * num_reactions)
    contlBounds = [0.0] * num_reactions
    contuBounds = [default_bound] * num_reactions
    #	  a+b   c+d          e
    lbs = lbs + binlBounds + contlBounds
    ubs = ubs + binuBounds + contuBounds

    # Type variables
    print '\tBuild variable types'
    a = ['continuous'] * num_reactions
    b = ['continuous'] * num_reactions
    c = ['binary'] * num_reactions
    d = ['binary'] * num_reactions
    e = ['continuous'] * num_reactions
    variable_types = a + b + c + d + e

    # Objectives
    newObjective = {}
    for key in objective:
        val = objective[key]
        newKey = (4 * num_reactions) + key
        newObjective[newKey] = val

    return(LHS, equalities, RHS, variable_types, lbs, ubs, newObjective)


def build_comparison_daaave(S, lbs1In, ubs1In, lbs2In, ubs2In, target,
                            default_bound=1000000):

    # Make bounds for positive formulation
    print '\tRecast bounds for positive fluxes'
    lbs1, ubs1 = _makeFluxesPositive(lbs1In, ubs1In)
    lbs2, ubs2 = _makeFluxesPositive(lbs2In, ubs2In)

    # The easiest way to write this formulation is to concatenate S and -S to
    # double the variables. The presolve deals with this very rapidly, so it's
    # not worth the hassle of the extra code
    sDimensions = S.get_shape()
    num_species, num_reactions = sDimensions
    print '\tS matrix is', num_species, 'x', num_reactions

    # Build LHS
    print '\tBuilding LHS'

    print '\tBuild Positive FBA block for S1'
    zeros_species_by_reactions = csr_matrix((num_species, num_reactions))
    zeros_reactions_by_reactions = csr_matrix((num_reactions, num_reactions))
    ones_reactions_by_reactions = sparse.eye(num_reactions).tocsr()
    negOnes_reactions_by_reactions = -1 * ones_reactions_by_reactions
    thousands_reactions_by_reactions = default_bound * \
        ones_reactions_by_reactions  # note that this should change if limits change

    Pstack1 = hstack([S, -1 * S, zeros_species_by_reactions,
                      zeros_species_by_reactions, zeros_species_by_reactions])
    Pstack2 = hstack([negOnes_reactions_by_reactions, zeros_reactions_by_reactions,
                      thousands_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack3 = hstack([zeros_reactions_by_reactions, negOnes_reactions_by_reactions,
                      zeros_reactions_by_reactions, thousands_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack4 = hstack([zeros_reactions_by_reactions, zeros_reactions_by_reactions,
                      ones_reactions_by_reactions, ones_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack5 = hstack([ones_reactions_by_reactions, ones_reactions_by_reactions,
                      zeros_reactions_by_reactions, zeros_reactions_by_reactions, negOnes_reactions_by_reactions])
    PosFBA_block1 = vstack([Pstack1, Pstack2, Pstack3, Pstack4, Pstack5])

    print '\tBuild Positive FBA block for S2'
    PosFBA_block2 = PosFBA_block1  # Seems pointless but we might need this later

    print '\tBuild MOMA block'
    target_reactions_by_reactions = csr_matrix(np.diagflat(target))
    negTarget_reactions_by_reactions = -1 * target_reactions_by_reactions

    # Build stacks
    print '\t\tBuilding stacks'
    Mstack1 = hstack([target_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions,
                      zeros_reactions_by_reactions, negOnes_reactions_by_reactions, negOnes_reactions_by_reactions, zeros_reactions_by_reactions])
    Mstack2 = hstack([negTarget_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions,
                      zeros_reactions_by_reactions, ones_reactions_by_reactions, zeros_reactions_by_reactions, negOnes_reactions_by_reactions])

    # Vertical stack
    MOMA_block = vstack([Mstack1, Mstack2])

    print '\tCombine blocks'
    zeros_upper_right = csr_matrix(
        (num_species + (4 * num_reactions), 7 * num_reactions))
    zeros_lower_left = csr_matrix((2 * num_reactions, 4 * num_reactions))
    zeros_middle_left = csr_matrix(
        (num_species + (4 * num_reactions), 5 * num_reactions))
    zeros_middle_right = csr_matrix(
        (num_species + (4 * num_reactions), 2 * num_reactions))

    upper_stack = hstack([PosFBA_block1, zeros_upper_right])
    middle_stack = hstack(
        [zeros_middle_left, PosFBA_block2, zeros_middle_right])
    lower_stack = hstack([zeros_lower_left, MOMA_block])

    # Process lower stack to zero out rows associated with -1 flagged reactions
    print '\t\tEmptying lower stack of data-less constraints'
    lower_stack = lower_stack.tocsr()
    count = 0
    for row in range(len(target)):
        if math.isnan(target[row]):
            # Clear constraints on excesses
            lower_stack.data[lower_stack.indptr[row]                             :lower_stack.indptr[row + 1]] = float(0)
            lower_stack.data[lower_stack.indptr[row + num_reactions]:lower_stack.indptr[row +
                                                                                        1 + num_reactions]] = float(0)  # Clear constraints on deficits
            count += 1
    print '\t\t\t', count * 2, 'data-less MoMA constraints removed'

    LHS = vstack([upper_stack, middle_stack, lower_stack])

    # Equalities
    print '\tBuilding equalities'
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
    equalities = stack1_eqs + stack2_eqs + stack3_eqs + stack4_eqs + stack5_eqs + stack6_eqs + \
        stack7_eqs + stack8_eqs + stack9_eqs + stack10_eqs + stack11_eqs + stack12_eqs

    # RHS
    print '\tBuild RHS'
    # Stacks 1-5 are associated with the positive FBA formulation for S1
    stack1_RHS = [0.0] * num_species
    stack2_RHS = [0.0] * num_reactions
    stack3_RHS = [0.0] * num_reactions
    stack4_RHS = [1.0] * num_reactions
    stack5_RHS = [0.0] * num_reactions
    # Stacks 6-10 are associated with the positive FBA formulation for S2
    stack6_RHS = [0.0] * num_species
    stack7_RHS = [0.0] * num_reactions
    stack8_RHS = [0.0] * num_reactions
    stack9_RHS = [1.0] * num_reactions
    stack10_RHS = [0.0] * num_reactions
    # Stacks 11-12 are associated with the MoMA formulation
    stack11_RHS = [0.0] * num_reactions
    stack12_RHS = [0.0] * num_reactions
    RHS = stack1_RHS + stack2_RHS + stack3_RHS + stack4_RHS + stack5_RHS + stack6_RHS + \
        stack7_RHS + stack8_RHS + stack9_RHS + stack10_RHS + stack11_RHS + stack12_RHS

    # Bounds
    print '\tBuild bounds'
    binlBounds = [0] * (2 * num_reactions)
    binuBounds = [1] * (2 * num_reactions)
    contlBounds = [0.0] * num_reactions
    contuBounds = [default_bound] * num_reactions

    # Bounds associated with Pos FBA blocks
    #	  	   a+b    c+d          e
    fba1_lbs = lbs1 + binlBounds + contlBounds
    fba1_ubs = ubs1 + binuBounds + contuBounds
    fba2_lbs = lbs2 + binlBounds + contlBounds
    fba2_ubs = ubs2 + binuBounds + contuBounds

    # Bounds associated with 1-norm minimisation
    #		   x		     d
    moma_lbs = contlBounds + contlBounds
    moma_ubs = contuBounds + contuBounds

    lbs = fba1_lbs + fba2_lbs + moma_lbs
    ubs = fba1_ubs + fba2_ubs + moma_ubs

    # Type variables
    print '\tBuild variable types'
    a = ['continuous'] * num_reactions
    b = ['continuous'] * num_reactions
    c = ['binary'] * num_reactions
    d = ['binary'] * num_reactions
    e = ['continuous'] * num_reactions
    x = ['continuous'] * num_reactions
    d = ['continuous'] * num_reactions
    variable_types = a + b + c + d + e + a + b + c + d + e + x + d

    # Objectives
    objective = {}
    for v in range(11 * num_reactions, 12 * num_reactions):
        objective[v] = 1
        objective[v + num_reactions] = 1

    print '\tBuilt...'
    print '\t\tLHS:', LHS.get_shape()
    print '\t\tVar:', len(variable_types)
    print '\t\tLBs:', len(lbs)
    print '\t\tUBs:', len(ubs)
    print '\t\tEqs:', len(equalities)
    print '\t\tRHS:', len(RHS)

    return(LHS, equalities, RHS, variable_types, lbs, ubs, objective)


def _buildSuperDaaave(S, lbsIn, ubsIn, target, stdevs, default_bound=1000000):
    '''
    SuperDaaave is an algorithm that seeks to correlate measured abundance and estimated flux.
    This is the Positive version of the algorithm that runs more quickly because all fluxes must be positive.
    '''

    # Make bounds for positive formulation
    print '\tReconfiguring bounds'
    lbs, ubs = _makeFluxesPositive(lbsIn, ubsIn)

    # The easiest way to write this formulation is to concatenate S and -S to
    # double the variables. The presolve deals with this very rapidly, so it's
    # not worth the hassle of the extra code
    sDimensions = S.get_shape()
    num_species, num_reactions = sDimensions
    print '\tS matrix is', num_species, 'x', num_reactions

    # Build LHS
    print '\tBuilding LHS'

    print '\tBuild Positive FBA block'
    zeros_species_by_reactions = csr_matrix((num_species, num_reactions))
    zeros_reactions_by_reactions = csr_matrix((num_reactions, num_reactions))
    ones_reactions_by_reactions = sparse.eye(num_reactions).tocsr()
    negOnes_reactions_by_reactions = -1 * ones_reactions_by_reactions
    thousands_reactions_by_reactions = default_bound * \
        ones_reactions_by_reactions  # note that this should change if limits change

    Pstack1 = hstack([S, -1 * S, zeros_species_by_reactions,
                      zeros_species_by_reactions, zeros_species_by_reactions])
    Pstack2 = hstack([negOnes_reactions_by_reactions, zeros_reactions_by_reactions,
                      thousands_reactions_by_reactions, zeros_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack3 = hstack([zeros_reactions_by_reactions, negOnes_reactions_by_reactions,
                      zeros_reactions_by_reactions, thousands_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack4 = hstack([zeros_reactions_by_reactions, zeros_reactions_by_reactions,
                      ones_reactions_by_reactions, ones_reactions_by_reactions, zeros_reactions_by_reactions])
    Pstack5 = hstack([ones_reactions_by_reactions, ones_reactions_by_reactions,
                      zeros_reactions_by_reactions, zeros_reactions_by_reactions, negOnes_reactions_by_reactions])
    PosFBA_block = vstack([Pstack1, Pstack2, Pstack3, Pstack4, Pstack5])

    print '\tBuild MOMA block'
    Mstack1 = hstack([ones_reactions_by_reactions, negOnes_reactions_by_reactions,
                      negOnes_reactions_by_reactions, zeros_reactions_by_reactions])
    Mstack2 = hstack([negOnes_reactions_by_reactions, ones_reactions_by_reactions,
                      zeros_reactions_by_reactions, negOnes_reactions_by_reactions])
    MOMA_block = vstack([Mstack1, Mstack2])

    print '\tCombine blocks'
    zeros_upper_right = csr_matrix(
        (num_species + (4 * num_reactions), 3 * num_reactions))
    zeros_lower_left = csr_matrix((2 * num_reactions, 4 * num_reactions))
    upper_stack = hstack([PosFBA_block, zeros_upper_right])
    lower_stack = hstack([zeros_lower_left, MOMA_block])
    LHS = vstack([upper_stack, lower_stack])

    # Equalities
    print '\tBuilding equalities'
    stack1_eqs = ['='] * num_species
    stack2_eqs = ['>'] * num_reactions
    stack3_eqs = ['>'] * num_reactions
    stack4_eqs = ['<'] * num_reactions
    stack5_eqs = ['='] * num_reactions
    stack6_eqs = ['<'] * num_reactions
    stack7_eqs = ['<'] * num_reactions
    equalities = stack1_eqs + stack2_eqs + stack3_eqs + \
        stack4_eqs + stack5_eqs + stack6_eqs + stack7_eqs

    # RHS
    print '\tBuild RHS'
    stack1_RHS = [0.0] * num_species
    stack2_RHS = [0.0] * num_reactions
    stack3_RHS = [0.0] * num_reactions
    stack4_RHS = [1.0] * num_reactions
    stack5_RHS = [0.0] * num_reactions
    stack6_RHS = [0.0] * num_reactions
    stack7_RHS = [0.0] * num_reactions
    RHS = stack1_RHS + stack2_RHS + stack3_RHS + \
        stack4_RHS + stack5_RHS + stack6_RHS + stack7_RHS

    # Bounds
    print '\tBuild bounds'
    binlBounds = [0] * (2 * num_reactions)
    binuBounds = [1] * (2 * num_reactions)
    contlBounds = [0.0] * num_reactions
    contuBounds = [default_bound] * num_reactions

    # Process target bounds to remove NaNs. Matrix block structure is
    # maintained, meaning NaN variables are completely free and some
    # constraints are inactive. These are left to the presolver to deal with.
    lTarget = []
    uTarget = []
    nanCount = 0
    for t in target:
        if math.isnan(t):  # Variable can take any value
            lTarget.append(0)
            uTarget.append(default_bound)
            nanCount += 1
        else:  # Variable is constrained to abundance data
            lTarget.append(t)
            uTarget.append(t)

    #	  a+b   c+d          e             t         x             d
    lbs = lbs + binlBounds + contlBounds + lTarget + contlBounds + contlBounds
    ubs = ubs + binuBounds + contuBounds + uTarget + contuBounds + contuBounds

    # Type variables
    print '\tBuild variable types'
    a = ['continuous'] * num_reactions
    b = ['continuous'] * num_reactions
    c = ['binary'] * num_reactions
    d = ['binary'] * num_reactions
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
        if math.isnan(stdevs[fetch]) == False:
            std = stdevs[fetch]
            # This accounts for STDEV and zeros
            obj = 1 / (std + (minimum_stdev / 2))
            objective[v] = obj
            objective[v + num_reactions] = obj
        fetch += 1

    print '\tBuilt...'
    print '\t\tLHS:', LHS.get_shape()
    print '\t\tVar:', len(variable_types), '=>', nanCount, 'NaNs opened'
    print '\t\tLBs:', len(lbs)
    print '\t\tUBs:', len(ubs)
    print '\t\tEqs:', len(equalities)
    print '\t\tRHS:', len(RHS)
    print '\t\tObj:', len(objective)

    return(LHS, equalities, RHS, variable_types, lbs, ubs, objective)


def _makeFluxesPositive(lbs, ubs):
    '''Converts lists of length num_reactions describing lower and upper bounds into lists
    of length 2*num_reactions describing lower and upper bounds associated with concatenated
    matrices S, -S.

    This preserves reaction directionality as written in the SBML (from lbs, ubs), but
    recasts all fluxes as zero to +ub
    '''
    posUbounds = []
    posLbounds = []
    negUbounds = []
    negLbounds = []
    typer = {}
    typer['Reversible'] = 0
    typer['Zero or Forwards'] = 0
    typer['Zero or Backwards'] = 0
    typer['Compulsory Backwards'] = 0
    typer['Compulsory Forwards'] = 0
    typer['Dead'] = 0
    typer['Other'] = 0

    i = 0
    for lb in lbs:
        ub = ubs[i]
        if lb < 0 and ub > 0:
            typer['Reversible'] += 1
            posUbounds.append(ub)  # Alive
            posLbounds.append(0.0)  # Alive
            negUbounds.append(-1 * lb)  # Alive
            negLbounds.append(0.0)  # Alive
        elif lb == 0 and ub > 0:
            typer['Zero or Forwards'] += 1
            posUbounds.append(ub)  # Alive
            posLbounds.append(0.0)  # Alive
            negUbounds.append(0.0)  # Dead
            negLbounds.append(0.0)  # Dead
        elif lb < 0 and ub == 0:
            typer['Zero or Backwards'] += 1
            posUbounds.append(0.0)  # Dead
            posLbounds.append(0.0)  # Dead
            negUbounds.append(-1 * lb)  # Alive
            negLbounds.append(0.0)  # Alive
        elif lb == 0 and ub == 0:
            typer['Dead'] += 1
            posUbounds.append(0.0)  # Dead
            posLbounds.append(0.0)  # Dead
            negUbounds.append(0.0)  # Dead
            negLbounds.append(0.0)  # Dead
        elif lb > 0 and ub > 0:
            typer['Compulsory Forwards'] += 1
            posUbounds.append(ub)  # Alive
            posLbounds.append(lb)  # Alive
            negUbounds.append(0.0)  # Dead
            negLbounds.append(0.0)  # Dead
        elif lb < 0 and ub < 0:
            typer['Compulsory Backwards'] += 1
            posUbounds.append(0.0)  # Dead
            posLbounds.append(0.0)  # Dead
            negUbounds.append(-1 * lb)  # Alive
            negLbounds.append(-1 * ub)  # Alive
        else:
            print 'What!? THIS SHOULD NOT HAPPEN!'
            typer['Other'] += 1
            print 'Other', reaction_names[i], lbs[i], ubs[i]
        i += 1

    sum = 0
    for key in typer:
        print '\t\t' + str(typer[key]) + '\t' + key
        sum = sum + typer[key]

    uBounds = posUbounds + negUbounds
    lBounds = posLbounds + negLbounds

    return(lBounds, uBounds)
