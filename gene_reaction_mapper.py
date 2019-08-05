'''
(c) University of Liverpool 2019

All rights reserved.
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import itertools


def map_genes_to_reactions(gprs, gene_to_expression, gene_to_expression_sd):
    '''TODO'''
    rxn_exp = []
    rxn_exp_sds = []

    for gpr in gprs:
        exp = []
        exp_sds = []

        for enzyme_complex in _parse_gpr(gpr):
            exp.append(min([(gene_to_expression[gene]
                             if gene in gene_to_expression
                             else 1e-6)
                            for gene in enzyme_complex]))

            exp_sds.append(min(
                [(gene_to_expression_sd[gene]
                  if gene in gene_to_expression_sd
                  else 1e-6)
                 for gene in enzyme_complex]))

        rxn_exp.append(max(exp) if exp else float('NaN'))
        rxn_exp_sds.append(max(exp_sds) if exp_sds else float('NaN'))

    return rxn_exp, rxn_exp_sds


def _add_gpr(gpr,
             gene_to_expression,
             gene_to_expression_sd,
             gurobipy_model,
             reactions,
             all_gene_instances,
             gene_to_gene_instances,
             rxn_to_gene_usage_to_gene_exp_sd):
    '''TODO'''
    enzyme_complexes = []

    # Add reaction instances:
    reaction_name = 'r' + str(len(reactions))
    r_var = gurobipy_model.addVar(name=reaction_name)
    reactions.append(r_var)

    reaction_gene_to_gene_instances = {}

    for i, enzyme_complex in enumerate(_parse_gpr(gpr)):
        # Add complex instances:
        _add_complex(i, enzyme_complex, reaction_name, enzyme_complexes,
                     gurobipy_model, all_gene_instances,
                     gene_to_gene_instances,
                     reaction_gene_to_gene_instances)

    # Add reaction-gene constraints (for standard deviation calculation):
    gene_usage_to_gene_expression_sd = []
    rxn_to_gene_usage_to_gene_exp_sd.append(gene_usage_to_gene_expression_sd)

    for gene, gene_instances in reaction_gene_to_gene_instances.items():
        r_g = gurobipy_model.addVar(name=gene + "_" + reaction_name)
        gene_expression = gene_to_expression[gene] \
            if gene in gene_to_expression and gene_to_expression[gene] > 0 \
            else 1e-6

        lin_expr = gurobipy.LinExpr(
            [1.0 / gene_expression] * len(gene_instances), gene_instances)

        gurobipy_model.update()
        gurobipy_model.addConstr(lin_expr, gurobipy.GRB.EQUAL, r_g)
        gene_usage_to_gene_expression_sd.append(
            (r_g, gene_to_expression_sd[gene]
             if gene in gene_to_expression_sd
             and gene_to_expression_sd[gene] > 0
             else 1e-6))

    # Add reaction constraints:
    lin_expr = gurobipy.LinExpr(
        [1.0] * len(enzyme_complexes), enzyme_complexes)
    gurobipy_model.update()
    gurobipy_model.addConstr(lin_expr, gurobipy.GRB.EQUAL, r_var)


def _add_complex(i, enzyme_complex, reaction_name, enzyme_complexes,
                 gurobipy_model, all_gene_instances, gene_to_gene_instances,
                 reaction_gene_to_gene_instances):
    '''Add complex.'''
    x_name = 'x' + str(i) + '_' + reaction_name
    x_var = gurobipy_model.addVar(name=x_name)
    enzyme_complexes.append(x_var)

    # Add gene instances:
    if isinstance(enzyme_complex, str):
        _add_gene(gurobipy_model,
                  enzyme_complex,
                  all_gene_instances,
                  gene_to_gene_instances,
                  reaction_gene_to_gene_instances,
                  x_var,
                  x_name)
    else:
        for gene in enzyme_complex:
            _add_gene(gurobipy_model,
                      gene,
                      all_gene_instances,
                      gene_to_gene_instances,
                      reaction_gene_to_gene_instances,
                      x_var,
                      x_name)


def _add_gene(m, gene, all_gene_instances, gene_to_gene_instances,
              reaction_gene_to_gene_instances, x, x_name):
    '''TODO'''
    if gene not in gene_to_gene_instances:
        gene_to_gene_instances[gene] = []

    if gene not in reaction_gene_to_gene_instances:
        reaction_gene_to_gene_instances[gene] = []

    g_var = m.addVar(name=gene + '_' + x_name)

    gene_to_gene_instances[gene].append(g_var)
    reaction_gene_to_gene_instances[gene].append(g_var)
    all_gene_instances.append(g_var)

    # Integrate new variables:
    m.update()

    # Add complex constraints:
    m.addConstr(x, gurobipy.GRB.LESS_EQUAL, g_var)

    return g_var


def _parse_gpr(gpr, consider_splice_variants=False):
    '''TODO'''
    gpr = gpr.replace('(', ' ( ')
    gpr = gpr.replace(')', ' ) ')
    gpr = gpr.replace(' AND ', ' and ')
    gpr = gpr.replace(' OR ', ' or ')

    isoenzymes = []
    _evaluate_statements(gpr.split(), isoenzymes, consider_splice_variants)

    # Replace strings with [strings]
    for i, term in enumerate(isoenzymes):
        if isinstance(term, str):
            isoenzymes[i] = [term]

    # Remove duplicates:
    isoenzymes.sort()
    return list(i for i, _ in itertools.groupby(isoenzymes))


def _evaluate_statements(tokens, isoenzymes, consider_splice_variants):
    '''TODO'''
    if not tokens:
        return

    if len(tokens) == 1:
        _add_isoenzyme(tokens[0], isoenzymes)
        return

    has_parentheses, l_paren, r_paren = _has_parentheses(tokens)

    if not has_parentheses:
        isoenzyme = _evaluate_statement(tokens, consider_splice_variants)
        _add_isoenzyme(isoenzyme, isoenzymes)
        return

    if l_paren + 1 == r_paren:  # Empty parenthesis
        tokens[l_paren:r_paren + 1] = []
    else:
        tokens[l_paren:r_paren + 1] = \
            [_evaluate_statement(tokens[l_paren + 1:r_paren],
                                 consider_splice_variants)]

    _evaluate_statements(tokens, isoenzymes, consider_splice_variants)


def _has_parentheses(token_lst):
    '''TODO'''
    left_lst = _find(token_lst, '(')

    if not left_lst:
        return False, -1, -1

    left = left_lst[-1]
    right = _find(token_lst, ')', left)[0]

    return True, left, right


def _find(lst, obj, start=0):
    '''TODO'''
    return [i for i, elem in enumerate(lst) if elem == obj and i >= start]


def _evaluate_statement(tokens, consider_splice_variants):
    '''TODO'''
    if len(tokens) == 1:
        return _consider_splice_variants(tokens[0], consider_splice_variants)

    if 'or' not in tokens:
        tokens = filter(lambda tokens: tokens != 'and', tokens)
        values = []

        for token in tokens:
            if isinstance(token, str):
                values.append(token)
            else:
                values.extend(token)

        return [_consider_splice_variants(v, consider_splice_variants)
                for v in values]

    or_index = tokens.index('or')
    lhs = _evaluate_statement(tokens[:or_index], consider_splice_variants)
    rhs = _evaluate_statement(tokens[or_index + 1:], consider_splice_variants)

    if isinstance(lhs, list) or isinstance(rhs, list):
        return lhs, rhs

    return lhs, rhs


def _consider_splice_variants(token, consider_splice_variants):
    '''TODO'''
    return token if consider_splice_variants or '.' not in token \
        else token[:token.find('.')]


def _add_isoenzyme(isoenzyme, isoenzymes):
    '''TODO'''
    if isinstance(isoenzyme, tuple):
        for term in isoenzyme:
            _add_isoenzyme(term, isoenzymes)
    elif isinstance(isoenzyme, list) and any(isinstance(term, tuple)
                                             for term in isoenzyme):
        # Special case for parenthesised OR terms, e.g. 'AAA and (BBB or CCC)'
        for i, term in enumerate(isoenzyme):
            if isinstance(term, str):
                isoenzyme[i] = [term]
            else:
                isoenzyme[i] = list(term)
        isoenzymes.extend([list(elem)
                           for elem in list(itertools.product(*isoenzyme))])
    else:
        isoenzymes.append(isoenzyme)
