'''
(c) University of Liverpool 2019

All rights reserved.
'''
import math


def extract_flux_patterns_from_positive_fba(lp_sol, num_reactions):
    '''Extract flux patterns.'''
    var_values = lp_sol.getVars()
    flux_vector = []

    for tic in range(4 * num_reactions, 5 * num_reactions):
        component = var_values[tic]
        flux_vector.append(component.X)

    return flux_vector


def extract_flux_patterns_from_super_daaave(lp_sol, num_reactions):
    '''Extract flux patterns.'''
    # Equivalent to fetching variables e
    fluxes = lp_sol.getVars()
    flux_vector = []

    for tic in range(num_reactions):
        flux_obj = fluxes[tic + (4 * num_reactions)]
        # Tidy up some weird sign errors from Gurobi (minus zero?)
        f_squared = math.pow(flux_obj.X, 2)
        f_sqrt = round(math.sqrt(f_squared), 4)
        flux_vector.append(f_sqrt)

    return flux_vector
