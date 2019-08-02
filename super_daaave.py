import LPSolver as solver
import LPanalyse as analyse
import Positive_Constraint_based_Modelling as pcbm
import SBMLDataMapper as sbml_data_mapper
import matplotlib.pyplot as plt


model = 'data/yeast_5.21_MCISB.xml'
data = 'data/genedata_75.txt'

print '- Fetch model and map in scaled expression data'
Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, rxn_mean_expression, rxn_stdev_expression = sbml_data_mapper._read_and_map_sbml_and_expression_data(
    model, data, 'D-glucose exchange', 16.5)

# Update bounds (this refers to the original directionality from the SBML file)
glucoseIndex = reaction_names.index('D-glucose exchange')
lbs[glucoseIndex] = -16.5
ubs[glucoseIndex] = -16.5

print '- Build SuperDaaave problem'
LHS, equalities, RHS, variable_types, lbs, ubs, objective = pcbm._buildSuperDaaave(
    Smatrix, lbs, ubs, rxn_mean_expression, rxn_stdev_expression)  # Build matrix problem

print '- Create Gurobi instance of problem'
LP = solver._create_GurobiLP(
    LHS, equalities, RHS, variable_types, lbs, ubs, objective)  # Build solver problem

print '- Solve SuperDaaave'
LPsolved = solver._run_GurobiLP(LP)  # Solve LP
flux = analyse._extract_flux_patterns_from_SuperDaaave(
    LPsolved, num_reactions)  # Extract flux pattern

# Plotting some results
plt.scatter(rxn_mean_expression, flux, s=rxn_stdev_expression)
plt.title('Transcript Abundance vs Nearest Possible Flux')
plt.xlabel('Abundance')
plt.ylabel('Flux')
plt.show()
