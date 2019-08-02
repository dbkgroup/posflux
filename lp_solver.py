'''
(c) University of Liverpool 2019

All rights reserved.
'''
from gurobipy import Model, GRB, LinExpr


def _run_gurobi_lp(lp):
    lp.optimize()
    return lp


def _create_gurobi_lp(lhs, equalities, rhs, variable_types, lbs, ubs, objs):
    '''
    lhs				is a CSR matrix of size c * v (constraints * variables)
                    with nonzero coefficients
    equalities		is a list of size c containing equalities, inequalities
    RHS				is a list of size c of values (for now no functions are included)
    variable_types	is a list of variables types
    lbs,ubs			are lists of size v containing lower and upper bounds
    objective		is the dictionary of non-zero elements and their weights
                    (virtually of size v, though in practice often smaller)
    '''
    # Create model
    lp = Model("LP")

    # Create variables, note all variables here are continuous
    print('\tCreate Gurobi variables')

    for var_index, variable_type in enumerate(variable_types):
        # Identify variable type
        if variable_type == 'continuous':
            var_type = GRB.CONTINUOUS
        elif variable_type == 'binary':
            var_type = GRB.BINARY
        else:
            var_type = 'unknown'

        # Identify whether objective
        if var_index in objs:
            objective_value = objs.get(var_index)
        else:
            objective_value = 0.0

        # Create variable
        lp.addVar(lbs[var_index], ubs[var_index], objective_value,
                  var_type, 'v' + str(var_index))

    # Integrate new variables
    lp.update()

    lpvars = lp.getVars()  # List of variable objects

    # Create constraints
    print('\tCreate Gurobi constraints')

    for cons_index, _ in enumerate(equalities):
        lin = []  # Coefficients
        refs = []  # Variables
        row = lhs.getrow(cons_index)
        _, variables = row.nonzero()

        lin = row.data

        for r in variables:
            # coeff = row.getcol(r).toarray()[0][0]
            refs.append(lpvars[r])

        if equalities[cons_index] == '>':
            equality = GRB.GREATER_EQUAL
        elif equalities[cons_index] == '<':
            equality = GRB.LESS_EQUAL
        elif equalities[cons_index] == '=':
            equality = GRB.EQUAL
        else:
            equality = ''

        lp.addConstr(LinExpr(lin, refs), equality,
                     rhs[cons_index], 'c' + str(cons_index))

    # Integrate new constraints
    lp.update()

    return lp
